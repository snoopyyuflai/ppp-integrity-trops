function initialize(obj,obs,corrData)

% Try to initialize the all in view filter

aivFilter = copy(obj.filter);

[measUsedAiv] = aivFilter.initialize(obs,corrData);


% Remove all subsets
obj.subsetFilters = repelem(navsu.estimators.pppFilter,0,0);
obj.subsetMeasOut = repelem(navsu.internal.MeasID,0,0);
obj.subsetPsat    = [];
obj.subsetId      = [];
obj.subsetWhichGroup  = [];
obj.subsetPsatList = [];
obj.subsetType = [];
obj.groups = repelem(navsu.internal.MeasID,0,0);
obj.groupsPsat = [];

% Need a bunch of GNSS satellites. please only proceed if we have enough
% satllites :)
if ~isempty(measUsedAiv)
    indsGnss = find(cat(1,measUsedAiv.TypeID) == navsu.internal.MeasEnum.GNSS);
    measUsedGnss = measUsedAiv(indsGnss);
    measUsedCode = measUsedGnss(find(cat(1,measUsedGnss.subtype) == navsu.internal.MeasEnum.Code));
    
    % cat(1,measUsedGnss.subtype) is used to called the whole subtype into
    % an array of measUsedCode
    
    if isempty(indsGnss)
        nSv = 0;
    else
        prnConst = [cat(1,measUsedCode.prn) cat(1,measUsedCode.const)];
        nSv =length(unique(prnConst,'rows'));
    end
    if nSv < 8
        % not enough number of satellite available!!!
        % differential ppp so needs twice of 4?
        aivFilter.initialized = 0;
    end
end

if aivFilter.initialized == 1  
    obj.aivFilter = aivFilter;
    
    % do the subset management
    obj.manageSubsets(measUsedAiv,obs,corrData,aivFilter)
    
    % Update the protection levels
    computePl(obj);
    
    % If the initialization was successful, then we can say so and keep going.
    obj.initialized = true;
end





end