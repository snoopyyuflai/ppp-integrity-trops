function [measId,extraInputs] = measUpdate(obj,epoch,obs,corrData,measRemovedSlip,varargin)


p = inputParser;

p.addParameter('measExclude',[]);
p.addParameter('extraInputs',[]);

% parse the results
parse(p, varargin{:});
res        = p.Results;
measExclude = res.measExclude;
extraInputs = res.extraInputs;

%% Pull a few things out of the filter object

% Run parameters
PARAMS = obj.PARAMS;

% State- has been propagated in the time update
statePropagated = obj.state;

% Receiver position
pos = obj.pos;

% Receiver velocity
vel = obj.vel;

%
covPropagated  = obj.cov;
nState = size(covPropagated,1);
epoch0 = epoch;

% Initialize
predMeas = [];
H = zeros(0,nState);
R = [];
measId  = [];       % MeasID
meas    = [];       % actual measurements
prnConstInds = [];
el = [];
az = [];

obj.resids      = [];
obj.measRemoved = [];
gnssMeas = [];

%% Loop through each measurement type and add it to the list!
for idx = 1:length(obs)
    obsi = obs{idx};
    
    switch obsi.type
        case navsu.internal.MeasEnum.GNSS
            [predMeasi,Hi,Ri,el,az,prnConstInds,measIdi,measi,...
                measIdRemovedLow,extraInputs,dop] = handleGnssMeas(obj,epoch0,obsi,...
                corrData,'extraInputs',extraInputs);
            gnssMeas = obsi;
        case navsu.internal.MeasEnum.Position
            [predMeasi,Hi,Ri,measIdi,measi] = handlePositionMeas(obj,obsi);
            
        case navsu.internal.MeasEnum.Velocity
            [predMeasi,Hi,Ri,measIdi,measi] = handleVelocityMeas(obj,obsi);
        otherwise
            continue;
    end
    [predMeas,H,R,measId,meas] = catMeas(predMeas,predMeasi,H,Hi,R,Ri,measId,measIdi,meas,measi);
end

%% Pseudomeasurements
% Apply constraint by add/remove measurements
if PARAMS.measUse.noVertVel 
    [predMeasi,Hi,Ri,measIdi,measi] = handleVehicleConstraintPseudomeas(obj);
    
    [predMeas,H,R,measId,meas] = catMeas(predMeas,predMeasi,H,Hi,R,Ri,measId,measIdi,meas,measi);
end

nMeas = size(H,1);

%% Measurement exclusion may have been forced from the outside- DO SOMETHING ABOUT IT!
if ~isempty(measExclude)
    measMask = false(size(H,1),1);
    for mdx = 1:length(measExclude)
        measMask = measMask | matches(measExclude(mdx),measId);
    end
    
    % Pull these values out from stuff :)
    H(measMask,:) = [];
    predMeas(measMask) = [];
    meas(measMask) = [];
    R(measMask,:) = [];
    R(:,measMask) = [];
    measId(measMask) = [];
end

%% Do the measurement update
if nMeas > 0
    % Measurement were available- do the update.
    % the update step of EKF, compute the current linearization of
    % jocobian matrices
    [H,delta_z,residsPost,K,~,~,~,measIdRemoved,measId] = ...
        navsu.ppp.measUpdateExclude(H,covPropagated,R,predMeas,PARAMS,meas,measId);
    
    % 9. Update state estimates
    stateNew = statePropagated + K * delta_z;
    
    % Update covariance
    cov = (eye(nState) - K * H) * covPropagated;
    
    obj.epochLastGnssUpdate = epoch;
    
    
else
    % No measurement update
    stateNew = statePropagated;
    
    cov = covPropagated;
    
    residsPost = [];
    measIdRemoved= [];
    
end

%% Update the position and velocity values
% obj.R_b_e = (eye(3) - navsu.geo.crossProdMatrix(stateNew(obj.INDS_STATE.ATTITUDE))) * obj.R_b_e;

% 'state' is the difference of vellocity and position (PPP), so need to
% subtract or add back
vel = vel - stateNew(obj.INDS_STATE.VEL); % velocity indices
pos = pos - stateNew(obj.INDS_STATE.POS); % position indices

obj.clockBias  = stateNew(obj.INDS_STATE.CLOCK_BIAS);
obj.clockDrift = stateNew(obj.INDS_STATE.CLOCK_DRIFT);

% put updated values into object
% obj.R_b_e = R_b_e;
obj.vel   = vel;
obj.pos   = pos;
obj.cov  = cov;
obj.state = stateNew;
obj.posPrevTc = pos;

% obj.imuBiasStates = obj.imuBiasStates + stateNew([obj.INDS_STATE.ACC_BIAS obj.INDS_STATE.W_BIAS]);

measIdRemovedFull = [measExclude(:); measIdRemoved(:);];

% Deal with resets if any of the removed measurements were carrier phases
if ~isempty(measIdRemovedFull)
    for jdx = 1:size(measIdRemovedFull,1)
        if ~isempty(measIdRemovedFull(jdx).TypeID) ...
                    && measIdRemovedFull(jdx).TypeID == navsu.internal.MeasEnum.GNSS ...
                    && measIdRemovedFull(jdx).subtype == navsu.internal.MeasEnum.Carrier
            % carrier phase- reset ambiguity by just removing the state
            obj.removeFlexState([measIdRemovedFull(jdx).prn ...
                                 measIdRemovedFull(jdx).const ...
                                 1 ...
                                 measIdRemovedFull(jdx).freq] );
        end
    end
end

%% Save things for output
if ~isempty(gnssMeas) && nMeas > 0
    measGnss = measId([measId.TypeID] == navsu.internal.MeasEnum.GNSS);
    if ~isempty(measGnss)
        obj.allSatsSeen = sortrows(unique([[[measGnss.prn]' [measGnss.const]']; ...
                                           obj.allSatsSeen],'rows'),2);
    end
    % Save residuals
    obj.resids.measId = measId;
    obj.resids.resids = residsPost;
    obj.resids.epochs = epoch0*ones(size(measId));
    
    obj.resids.el = el;
    obj.resids.az = az;
    obj.resids.prnConstInds = prnConstInds;
    obj.resids.epochsElAz = epoch0*ones(size(el));
    
    obj.resids.dop = dop;
    obj.resids.dopEpoch = epoch0;
    
    % Only actually keeping one of the low measurements per satelite
    if ~isempty(measIdRemovedLow)
        prnConstLow = [[measIdRemovedLow.prn]' [measIdRemovedLow.const]'];
        [~,indsLowUn] = unique(prnConstLow,'rows');
        measIdRemovedLow = measIdRemovedLow(indsLowUn);
    else
        measIdRemovedLow = [];
    end
    
    measRemovedIdAny = [measIdRemovedLow; measIdRemoved; measRemovedSlip];
    measRemovedReason = [obj.RemovedEl*ones(size(measIdRemovedLow)); ...
                         obj.RemovedResid*ones(size(measIdRemoved));...
                         obj.RemovedSlip*ones(size(measRemovedSlip))];
    
    obj.measRemoved.id = measRemovedIdAny;
    obj.measRemoved.reason = measRemovedReason;
    obj.measRemoved.epoch  = epoch0*ones(size(measRemovedReason));
end


end


function [predMeas,H,R,measId,meas] = catMeas(predMeas,predMeasi,H,Hi,R,Ri,measId,measIdi,meas,measi)


% concatenate the measurement information!
% Add everything
predMeas = [predMeas; predMeasi];
H         = [H; Hi];
measId    = [measId; measIdi];
meas      = [meas; measi];

R2 = zeros(size(R,1)+size(Ri,1),size(R,1)+size(Ri,1));
R2(1:size(R,1),1:size(R,1)) = R;
R2((end-(size(Ri,1)-1)):end,(end-(size(Ri,1)-1)):end) = Ri;
R = R2;


end