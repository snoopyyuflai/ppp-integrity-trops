function computePl(obj)

PARAMS = obj.PARAMS;
% Pull out position delta and covariances for all necessary subsets
% aivFlag = [obj.ssList.aivFlag];
ssType  = obj.subsetType;
pSats   = obj.subsetPsat;

indsSsPrimary = find(ssType == 1);  % Index of subset
nSsPrimary = length(indsSsPrimary); % Number of subset

posSs  = nan(3,nSsPrimary);
covSs  = nan(3,3,nSsPrimary);
pSatSs = pSats(indsSsPrimary);

posAiv = obj.aivFilter.pos;
covAiv    = obj.aivFilter.cov(obj.aivFilter.INDS_STATE.POS,obj.aivFilter.INDS_STATE.POS);

% pull current estimate of positions and covariances out from subsetfilter
for idx = 1:nSsPrimary
    posSs(:,idx)   =  obj.subsetFilters(indsSsPrimary(idx)).pos;
    
    % INDS_STATE = indices of position in the whole covariance matrix
    covSs(:,:,idx) = obj.subsetFilters(indsSsPrimary(idx)).cov(...
        obj.subsetFilters(indsSsPrimary(idx)).INDS_STATE.POS,...
        obj.subsetFilters(indsSsPrimary(idx)).INDS_STATE.POS);
    
%     obj.ssList(indsSsPrimary(idx)).ekf.cov(obj.ssList(indsSsPrimary(idx)).ekf.INDS_STATE.POS, ...
%         obj.ssList(indsSsPrimary(idx)).ekf.INDS_STATE.POS);
end

dPos = posAiv-posSs;

%% Rotate everything to ENU
llh = navsu.geo.xyz2llh(posAiv');

[~,R] = navsu.geo.xyz2enu([0 0 0],llh(1)*pi/180,llh(2)*pi/180);
% R = Rotational Matrix
% Convert the position deltas and the covariances to ENU
dPosEnu = R*dPos;
covAivEnu = R*covAiv*R';
covSsEnu  = nan(size(covSs));
for idx = 1:size(covSsEnu,3)
    covSsEnu(:,:,idx) = R*squeeze(covSs(:,:,idx))*R';
end

%% Compute protection levels!
PFA    = obj.PARAMS.pfa; %false alert
PHMI   = obj.PARAMS.phmi; %maximum integrity risk
PL_TOL = obj.PARAMS.pl_prec; %protection level precision

covAll = cat(3,covAivEnu,covSsEnu); % concatenate aiv and all subset cov
dxi     = [zeros(3,1) dPosEnu]'; % the 'dy'
p_subset = [1; pSatSs];
sigma = [squeeze(covAll(1,1,1:end)) squeeze(covAll(2,2,1:end)) squeeze(covAll(3,3,1:end))].^(1/2);

% pull out identical subsets
ssIdentical = find(sqrt(sigma(2:end,1).^2 - sigma(1,1).^2) == 0);
sigma(ssIdentical+1,:) = [];
dxi(ssIdentical+1,:) = [];
p_subset(ssIdentical+1)= [];

nSubsetsi = size(sigma,1);

Tmult = norminv(PFA/3/nSubsetsi);

% Initialize outputs
anyFail = [0 0 0];
pl_ss_exact_optfa = nan(3,1);
pl_ss_approx1_optfa = nan(3,1);
pl_ss_approx2_optfa = nan(3,1);
pl_ss_exact = nan(3,1);
pl_ss_approx1 = nan(3,1);
pl_ss_approx2 = nan(3,1);
pl_chi2_exact = nan(3,1);
pl_chi2_approx1 = nan(3,1);
pl_chi2_approx2 = nan(3,1);

% Saving some things in case FDE is desired
dpos = dxi;
thresh = nan(size(dpos));

computeExact = 1; % ?
for ddx = 1:3
    %All in view sigma
    sigma_v_0 = sigma(1,ddx);
    
    %Solution separation sigma
    sigma_v_ss = sqrt(sigma(:,ddx).^2 - sigma(1,ddx).^2);
    
    % if there are any imaginary things... just kick it out.
    if any(imag(sigma_v_ss)) % imaginary partj
        sigma_v_ss = abs(sigma_v_ss);
        %         error('Subset covariance smaller than AIV')
    end
    
    %%%%%%%% PLs based on solution separation test statistic %%%%%%%%
    PFA_1coord  = PFA/3;
    PHMI_1coord = PHMI/3;
    
    % Optimized allocation of false alert %
    Topt = ones(nSubsetsi,1)*solSep.pl.compute_protection_level_v2(sigma_v_ss(2:end),0*sigma_v_ss(2:end), 0*sigma_v_ss(2:end)+2, PFA_1coord, PL_TOL);
    Topt(1) = 0;
    thresh(:,ddx) = Topt;
    if ~any(abs(dxi(2:end,ddx)) > Topt(2:end))
        if computeExact && PARAMS.pl_ss_exact_optfa
            pl_ss_exact_optfa(ddx) = solSep.pl.ss_exact_vpl(sigma_v_0, sigma_v_ss, p_subset, Topt, PHMI_1coord,PL_TOL);
        end
        if PARAMS.pl_ss_approx1_optfa
            vpli = solSep.pl.ss_approx1_vpl(sigma_v_0, sigma_v_ss, p_subset, Topt, PHMI_1coord,PL_TOL);
            
            pl_ss_approx1_optfa(ddx) = vpli;
        end
        
        if PARAMS.pl_ss_approx2_optfa
            pl_ss_approx2_optfa(ddx) = solSep.pl.ss_approx2_vpl(sigma_v_0, sigma_v_ss, p_subset, Topt, PHMI_1coord,PL_TOL);
        end
    else
        anyFail(ddx) = 1;
    end
    
    % Equal allocation of false alert %
    T   = -norminv(.5*PFA_1coord/(nSubsetsi-1))*sigma_v_ss;
    if ~any(abs(dxi(2:end,ddx)) > Topt(2:end))
        multi = 1;
    else
        multi = -1;
    end
    
    if computeExact && PARAMS.pl_ss_exact
        pl_ss_exact(ddx) = multi*solSep.pl.ss_exact_vpl(sigma_v_0, sigma_v_ss, p_subset, T, PHMI_1coord,PL_TOL);
    end
    
    if PARAMS.pl_ss_approx1
        pl_ss_approx1(ddx) = multi*solSep.pl.ss_approx1_vpl(sigma_v_0, sigma_v_ss, p_subset, T, PHMI_1coord,PL_TOL);
    end
    if PARAMS.pl_ss_approx2
        pl_ss_approx2(ddx) = multi*solSep.pl.ss_approx2_vpl(sigma_v_0, sigma_v_ss, p_subset, T, PHMI_1coord,PL_TOL);
    end
    
end
if PARAMS.pl_ss_exact_optfa
    obj.pl = pl_ss_exact_optfa;
    
elseif PARAMS.pl_ss_approx1_optfa
	obj.pl = pl_ss_approx1_optfa; 
    
elseif PARAMS.pl_ss_approx2_optfa
    obj.pl = pl_ss_approx2_optfa;
    
elseif PARAMS.pl_ss_exact
    obj.pl = pl_ss_exact;
    
elseif PARAMS.pl_ss_approx1
    obj.pl = pl_ss_approx1;
    
elseif PARAMS.pl_ss_approx2
    obj.pl = pl_ss_approx2;
end


%% Rotate everything to local frame
if 0
llh = xyz2llh(posAiv');

% [~,R] = XYZ2ENU([0 0 0],llh(1)*pi/180,llh(2)*pi/180);
if ~isempty(obj.ssList(1,1).ekf.lastAccMeas)
    R = obj.ssList(1,1).ekf.R_b_e';
else
    vel0 = obj.ssList(1,1).ekf.vel;
    pos0 = obj.ssList(1,1).ekf.pos;
    xi  = -vel0./norm(vel0);
    z0i = pos0./norm(pos0);
    yi  = -cross(xi,z0i);
    zi = cross(xi,yi);
    
    R  = [xi yi zi]';
end

% Convert the position deltas and the covariances to ENU
dPosEnu = R*dPos;
covAivEnu = R*covAiv*R';
covSsEnu  = nan(size(covSs));
for idx = 1:size(covSsEnu,3)
    covSsEnu(:,:,idx) = R*squeeze(covSs(:,:,idx))*R';
end

end
%% Compute protection levels in local frame
if 0
    PFA    = 1e-6; %false alert
    PHMI   = 1e-7; %maximum integrity risk
    PL_TOL = 0.1; %protection level precision
    
    covAll = cat(3,covAivEnu,covSsEnu);
    dxi     = [zeros(3,1) dPosEnu]';
    p_subset = [1 pSatSs]';
    sigma = [squeeze(covAll(1,1,1:end)) squeeze(covAll(2,2,1:end)) squeeze(covAll(3,3,1:end))].^(1/2);
    
    % pull out identical subsets
    ssIdentical = find(sqrt(sigma(2:end,1).^2 - sigma(1,1).^2) == 0);
    sigma(ssIdentical+1,:) = [];
    dxi(ssIdentical+1,:) = [];
    p_subset(ssIdentical+1)= [];
    
    nSubsetsi = size(sigma,1);
    
    Tmult = norminv(PFA/3/nSubsetsi);
    
    % Initialize outputs
    anyFailLoc = [0 0 0];
    pl_ss_approx1_optfa = nan(3,1);
    
    for ddx = 1:3
        %All in view sigma
        sigma_v_0 = sigma(1,ddx);
        
        %Solution separation sigma
        sigma_v_ss = sqrt(sigma(:,ddx).^2 - sigma(1,ddx).^2);
        
        % if there are any imaginary things... just kick it out.
        if any(imag(sigma_v_ss))
            sigma_v_ss = abs(sigma_v_ss);
        end
        
        %%%%%%%% PLs based on solution separation test statistic %%%%%%%%
        PFA_1coord  = PFA/3;
        PHMI_1coord = PHMI/3;
        
        % Optimized allocation of false alert %
        Topt = ones(nSubsetsi,1)*compute_protection_level_v2(sigma_v_ss(2:end),0*sigma_v_ss(2:end), 0*sigma_v_ss(2:end)+2, PFA_1coord, PL_TOL);
        Topt(1) = 0;
        if ~any(abs(dxi(2:end,ddx)) > Topt(2:end))
            
            if PARAMS.solSep.pl_ss_approx1_optfa
                vpli = ss_approx1_vpl(sigma_v_0, sigma_v_ss, p_subset, Topt, PHMI_1coord,PL_TOL);
                
                pl_ss_approx1_optfa(ddx) = vpli;
            end
        else
            anyFailLoc(ddx) = 1;
        end
    end
    
    plLoc = pl_ss_approx1_optfa;
end
%% If there was a failure, do something about it!
if any(anyFail) && 0
    switch PARAMS.solSep.fdeMode
        case 'reinit'
            % Determine which was the faulty subset
            [~,indWorst] =  max(max((abs(dpos)./thresh),[],2));
            
            measExclude = cat(1,obj.measOut,obj.ssList(indWorst).measOut);
            obj = solSepWrap(PARAMS);
            obj.measOut = measExclude;
            
        case 'background'
            'fdsafa';
            % Create a new ss object
            % Determine which was the faulty subset
            [~,indWorst] =  max(max((abs(dpos)./thresh),[],2));
            
            ss2 = solSepWrap(PARAMS);
            
            % Keep the same groups
            ss2.ssGroups = obj.ssGroups;
            measExclude = obj.ssList(indWorst).measOut;
            %             measExclude = cat(1,obj.measOut,obj.ssList(indWorst).measOut);
            ss2.measOut = cat(1,obj.measOut,obj.ssList(indWorst).measOut);
            
            
            measExcludeRow0 = reshape(measExclude,size(measExclude,1)*size(measExclude,2),size(measExclude,3));
            measExcludeRow0(isnan(measExcludeRow0(:,1)),:) = [];
            % The first subset is the all in view and is the the one that
            % excluded the fault
            ss2.ssList(1) = obj.ssList(indWorst);
            
            % Go through and find each of the background subset
            ssTypes = [obj.ssList.subsetType]';
            indsBg = find(ssTypes == 2);
            measOutBg = obj.ssList(indsBg).measOutRows;
            indsNewSs = find(cellfun(@(x) any(ismember(x,measExcludeRow0,'rows')),measOutBg));
            
            ss2.ssList(2:(length(indsNewSs)+1)) = obj.ssList(indsBg(indsNewSs));
            
            
            ss2.initialized = true;
            ss2.lastID = obj.lastID;
            % Some things need to be cleaned up-
            ss2.ssList(1).aivFlag = true;
            ss2.ssList(1).subsetType = 0;
            ss2.ssList(1).measOut = [];
            ss2.ssList(1).pSat = 1;
            ss2.ssList(1).whichGroup = [];
            ss2.ssList(1).pSatList = [];
            %             ss2.ssList(1).ekf.posPrevTc = [];
            
            % Change the subset types of each subset
            for idx = 2:length(ss2.ssList)
                ss2.ssList(idx).subsetType = 1;
                ss2.ssList(idx).aivFlag = false;
                % actually need to remove the faulted measurement group from all
                % of the measurement exclusion info. Probably can just
                % remove the group as well?
                
                measOuti    = ss2.ssList(idx).measOut;
                whichGroupi = ss2.ssList(idx).whichGroup;
                pSatListi   = ss2.ssList(idx).pSatList;
                
                % bump everything up
                for rdx = 1:size(measExcludeRow0)
                    indRemove = find(ismember(permute(measOuti,[1 3 2]),measExcludeRow0(rdx,:),'rows'));
                    
                    % Bump things up for the measurement exclusion
                    measOuti(indRemove:(end-1),:,:) = measOuti((indRemove+1):end,:,:);
                    measOuti(end,:,:) = nan;
                    
                    % Bump elements up for the groups
                    whichGroupi(indRemove:(end-1)) = whichGroupi((indRemove+1):end);
                    whichGroupi(end) = nan;
                    
                    % Bump things up for the psat
                    pSatListi(indRemove:(end-1)) = pSatListi((indRemove+1):end);
                    pSatListi(end) = 0;
                    
                    % Put it all back in the solution separation object!
                    ss2.ssList(idx).measOut    = measOuti;
                    ss2.ssList(idx).whichGroup = whichGroupi;
                    ss2.ssList(idx).pSatList   = pSatListi;
                end
                
                % compute psat here, please
                groupsi = unique(ss2.ssList(idx).whichGroup);
                groupsi(isnan(groupsi)) = [];
                pSatGroupsi = zeros(size(groupsi));
                for gdx = 1:length(pSatGroupsi)
                    groupsi(gdx) = 1-prod(1-ss2.ssList(idx).pSatList(ss2.ssList(idx).whichGroup == groupsi(gdx)));
                end
                ss2.ssList(idx).pSat = prod(groupsi);
            end
            
            obj = ss2;
            
    end
end













end