function measRemovedSlip = manageStatesMulti(obj,epoch,obs)

multiFilter = length(obj) > 1;

if multiFilter
    objList = obj;
    obj = objList(1);
    nFilter = length(objList);
else
    objList = [];
end

PARAMS = obj.PARAMS;

gnssMeas = navsu.ppp.pullMeasFromList(obs, navsu.internal.MeasEnum.GNSS);

gnssMeas = navsu.ppp.measMask(gnssMeas, PARAMS.measMask);


if isempty(gnssMeas)
    % Currently, we need a GNSS measurement in order to proceed.
    measRemovedSlip = [];
    return;
end

% Check for cycle slips so that these can be removed and reset
% cycle slips: Discontinuity in phase lock loop
measRemovedSlip = obj.checkCycleSlips(epoch, gnssMeas, PARAMS);

if multiFilter && ~isempty(measRemovedSlip)
    for fdx = 2:nFilter
        for mdx = 1:size(measRemovedSlip,1)
            objList(fdx).removeFlexState(measRemovedSlip(mdx,:));
        end
    end
end

%% Include new states in the filter, if there is new satellite
% State types:
% 1 = carrier phase ambiguity- one per CP
% 2 = TEC state- one per LOS
% 3 = Satellite specific DCB
% 4 = Code phase multipath
% 5 = Carrier phase multipath
% 6 = Ephemeris error

stateTypes = {'cp'; ... % 1  
    'L1DELAYSTATE'; ... % 2
    'TECSTATE';     ... % 2
    'RX_DCB_GLO';...    % 3
    'RX_DCB_GPS';...    % 3
    'MP_CODE';...       % 4
    'MP_CARR';...       % 5
    'EPH'};             % 6
  
stateUse   = [true; ...
    PARAMS.states.iono && strcmp(PARAMS.states.ionoMode, 'L1DELAYSTATE'); ...
    PARAMS.states.iono && strcmp(PARAMS.states.ionoMode, 'TECSTATE'); ...
    PARAMS.states.RX_DCB_GLO;...
    PARAMS.states.RX_DCB_GPS;...
    PARAMS.states.MP_CODE; ...
    PARAMS.states.MP_CARR; ...
    PARAMS.states.EPH];

for sdx = 1:length(stateTypes)
    if ~stateUse(sdx)
        continue;
    end
    
    % Given the state that is being examined, the currently available
    % measurements, and the current states being tracked, determine what,
    % if any, new states need to be added and what can be removed. 
    [measInfoAvail,stateNotNeeded] = navsu.ppp.stateAddInfo(stateTypes{sdx},gnssMeas,obj);
    
    if ~isempty(stateNotNeeded)
        % Remove from the state estimate, covariance, FLEX_STATES, and
        % FLEX_STATES_INFO
        obj.state(obj.INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
        obj.cov(obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.cov(:,obj.INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
        obj.INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
        obj.INDS_STATE.FLEX_STATES = obj.INDS_STATE.FLEX_STATE_MIN-1+...
            (1:size(obj.INDS_STATE.FLEX_STATES_INFO,1))';
        
        if multiFilter
            for fdx = 2:nFilter
                objList(fdx).state(objList(fdx).INDS_STATE.FLEX_STATES(stateNotNeeded)) = [];
                objList(fdx).cov(objList(fdx).INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
                objList(fdx).cov(:,objList(fdx).INDS_STATE.FLEX_STATES(stateNotNeeded),:) = [];
                objList(fdx).INDS_STATE.FLEX_STATES_INFO(stateNotNeeded,:) = [];
                objList(fdx).INDS_STATE.FLEX_STATES = objList(fdx).INDS_STATE.FLEX_STATE_MIN-1+...
                    (1:size(objList(fdx).INDS_STATE.FLEX_STATES_INFO,1))';
            end
        end
    end    
    
    for idx = 1:size(measInfoAvail,1)
        % Add this new state
        infoAdd = measInfoAvail(idx,:);

        % Pull initial values for the state and covariance
        [stateAdd,covAdd] = navsu.ppp.initStateCov(stateTypes{sdx}, infoAdd, PARAMS, gnssMeas);
        
        indStateAdd = obj.INDS_STATE.FLEX_STATE_MIN+length(obj.INDS_STATE.FLEX_STATES);
        
        % state
        obj.state = [obj.state; stateAdd];
        
        % covariance
        covNew = zeros(indStateAdd);
        covNew(1:size(obj.cov,1),1:size(obj.cov,1)) = obj.cov;
        covNew(indStateAdd,indStateAdd) = covAdd;
        obj.cov = covNew;
        
        % FLEX_STATE_INFO
        obj.INDS_STATE.FLEX_STATES_INFO = [obj.INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
        
        % FLEX_STATE index
        obj.INDS_STATE.FLEX_STATES = [obj.INDS_STATE.FLEX_STATES; indStateAdd];
        
        if multiFilter
            for fdx = 2:nFilter
                objList(fdx).state = [objList(fdx).state; stateAdd];
                
                % covariance
                covNew = zeros(indStateAdd);
                covNew(1:size(objList(fdx).cov,1),1:size(objList(fdx).cov,1)) = objList(fdx).cov;
                covNew(indStateAdd,indStateAdd) = covAdd;
                objList(fdx).cov = covNew;
                
                % FLEX_STATE_INFO
                objList(fdx).INDS_STATE.FLEX_STATES_INFO = [objList(fdx).INDS_STATE.FLEX_STATES_INFO(:,1:4); infoAdd];
                
                % FLEX_STATE index
                objList(fdx).INDS_STATE.FLEX_STATES = [objList(fdx).INDS_STATE.FLEX_STATES; indStateAdd];
            end
        end
    end
end


end