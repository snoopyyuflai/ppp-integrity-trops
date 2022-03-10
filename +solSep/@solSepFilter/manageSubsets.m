function manageSubsets(obj,measId,obs,corrData,filterSpawn)

if isempty(measId)
    return;
end

faultResetSubsets =  strcmp(obj.PARAMS.fdeMode,'background');

%% Build necessary subsets from the currently available measurements
ssBasis = obj.PARAMS.ssBasis; % basics of ss, Psat + GNSS

ssNeedList = [];   % actually needed in ss
psatNeedList = [];
% loop thorught all measurement type, e.g. gnss, positions, velocity...
for idx = 1:length(ssBasis.measId)
    ssBasisIdi = ssBasis.measId(idx);
    psatBasisi = ssBasis.psat(idx);
    
    % pull measurement type out
    measIdType = measId(cat(1,measId.TypeID) == ssBasisIdi.TypeID);
    
    if isempty(measIdType)
        continue;
    end
    
    % This is fairly clumsy- should find a better way to handle this in the
    % future :)
    switch ssBasisIdi.TypeID
        case navsu.internal.MeasEnum.GNSS
            %           %  prn  |  const |  freq  | subtype
            measInfo = [cat(1,measIdType.prn) cat(1,measIdType.const)  cat(1,measIdType.freq) double(cat(1,measIdType.subtype))];
            basisInfo = repmat([ssBasisIdi.prn ssBasisIdi.const ssBasisIdi.freq double(ssBasisIdi.subtype)],size(measInfo,1),1);
            
            % Filter out unwanted and repeated measurements
            
            % For inf values, do nothing
            
            % Zero out when we don't care about the specific value here
            measInfo(basisInfo == 0) = 0;
            
            % For specific values, pick out measurements with that value.
            measInfo((basisInfo ~= 0) & (~isinf(basisInfo)) & measInfo ~= basisInfo) = nan;
            % remove rows with nans
            measInfo(any(isnan(measInfo),2),:) = [];
            
            % remove repeated measurements
            measInfo = unique(measInfo,'rows');
            
            % Make measurement IDs and add them to the list.
            ssNeedList = cat(1,ssNeedList,navsu.internal.MeasIdGnss(measInfo(:,1),measInfo(:,2),measInfo(:,3),measInfo(:,4)));
            psatNeedList = cat(1,psatNeedList, repmat(psatBasisi,size(measInfo,1),1));
            
        case navsu.internal.MeasEnum.Position
            %           %  id  |  xyz
            measInfo = [double(cat(1,measIdType.id)) double(cat(1,measIdType.xyz)) ];
            basisInfo = repmat([ssBasisIdi.id double(ssBasisIdi.xyz)],size(measInfo,1),1);
            
            % For inf values, do nothing
            
            % Zero out when we don't care about the specific value here
            measInfo(basisInfo == 0) = 0;
            
            % For specific values, pick out measurements with that value.
            measInfo((basisInfo ~= 0) & (~isinf(basisInfo)) & measInfo ~= basisInfo) = nan;
            % remove rows with nans
            measInfo(any(isnan(measInfo),2),:) = [];
            
            measInfo = unique(measInfo,'rows');
            
            ssNeedList = cat(1,ssNeedList,navsu.internal.MeasIdPos(measInfo(:,1),measInfo(:,2)));
            psatNeedList = cat(1,psatNeedList, repmat(psatBasisi,size(measInfo,1),1));
        case navsu.internal.MeasEnum.Velocity
            %           %  id  |  xyz
            measInfo = [double(cat(1,measIdType.id)) double(cat(1,measIdType.xyz)) ];
            basisInfo = repmat([ssBasisIdi.id double(ssBasisIdi.xyz)],size(measInfo,1),1);
            
            % For inf values, do nothing
            
            % Zero out when we don't care about the specific value here
            measInfo(basisInfo == 0) = 0;
            
            % For specific values, pick out measurements with that value.
            measInfo((basisInfo ~= 0) & (~isinf(basisInfo)) & measInfo ~= basisInfo) = nan;
            % remove rows with nans
            measInfo(any(isnan(measInfo),2),:) = [];
            
            measInfo = unique(measInfo,'rows');
            
            ssNeedList = cat(1,ssNeedList,navsu.internal.MeasIdVel(measInfo(:,1),measInfo(:,2)));
            psatNeedList = cat(1,psatNeedList, repmat(psatBasisi,size(measInfo,1),1));
    end
end

%% Groups lol
% Construct groups, not subset
% Group: Put fault together, or it can say it convert faults into another
% type by grouping number of them to reduce total fault number, and as a
% result reduce the associated number of subset

% iterate through all satellites/measurements, and put them into group, the
% size of group is specified previously
for idx = 1:length(ssNeedList)
    % if there is any satellite/fault being not already added to the group?
    % more of a sanity check
    if ~any(ssNeedList(idx) == obj.groups(:)) 
        
        % check if there's space to add it to a group
        % Find if there is any empty satellite/measurement id in the group,
        % if there is, meaning there is space to add current
        % satellite/measurements in
        indAdd  = min(find(obj.groups == navsu.internal.MeasID));
        
        if isempty(indAdd)
            % need to add a new group/ need to construct new group
            groupAdd = repelem(navsu.internal.MeasID,1,obj.PARAMS.subsetGroupSize);
            groupAdd(1) = ssNeedList(idx);
            
            % also include the probability of 'satellite' failure
            psatAdd = repmat(0,1,obj.PARAMS.subsetGroupSize);
            psatAdd(1) = psatNeedList(idx);
            
            obj.groups = cat(1,obj.groups,groupAdd);
            obj.groupsPsat = cat(1,obj.groupsPsat,psatAdd);
            
        else
            % add (satellite/fault) to an existing group
            
            
            obj.groups(indAdd) = ssNeedList(idx);
            obj.groupsPsat(indAdd) = psatNeedList(idx);
            
            
            %%%%%%%%%% THIS SHOULD ACTUALLY BE ADDED TO EXISTING SUBSETS AS
            %%%%%%%%%% WELL- LOL PLEASE FIX ME EVENTUALLY HAHA :)
        end
    end
end

% Probility of failure of any member of each group
% Probability of falures occur in each group, given the probability of one
% satellite/fault faliure, and number of fault in each group
% prod(1-obj.groupsPsat,2) := Probability of all satellites/faults in each
% group are not fault
pSatGroups = 1-prod(1-obj.groupsPsat,2);


%% Build the list of subsets to create, based on content of groups just created
nGroupsOut = obj.PARAMS.nOutSubset+faultResetSubsets;
% one fault might have multiple measurements, e.g. subsetGroupSize = 2
nMeasOutMax = nGroupsOut.*obj.PARAMS.subsetGroupSize;

ssNeeded     = [];
pSatNeeded   = [];
ssTypeNeeded = [];
ssWhichGroup = [];
ssPsatList   = [];

for idx = 1:nGroupsOut
    indsi = nchoosek(1:size(obj.groups,1),idx); % idex of faults to be excluded
    
    ssGroupsi = []; % Satellite Idex of Group in ss
    pSatSsi = ones(size(indsi,1),1);
    ssPsatListi = [];
    
    for jdx = 1:size(indsi,2)
        ssGroupsi = cat(2,ssGroupsi,obj.groups(indsi(:,jdx),:));
        ssPsatListi = cat(2,ssPsatListi,obj.groupsPsat(indsi(:,jdx),:));
        pSatSsi = pSatSsi.*pSatGroups(indsi(:,jdx));
    end
    
    % Put these in a nan-padded full size matrix to make concatenation
    % easier
    ssGroupsPad = repelem(navsu.internal.MeasID,size(ssGroupsi,1),nMeasOutMax);
    ssGroupsPad(:,1:size(ssGroupsi,2)) = ssGroupsi;
    
    ssPsatListPad = nan(size(ssPsatListi,1),nMeasOutMax);
    ssPsatListPad(:,1:size(ssGroupsi,2)) = ssPsatListi;
    
    % group 1 + group 2...so need to identify this measurement is in which
    % group
    ssWhichGroupi = repmat(kron(1:obj.PARAMS.subsetGroupSize,ones(1,nGroupsOut)),size(ssGroupsi,1),1);
    ssWhichGroupi(ssGroupsPad == navsu.internal.MeasID) = nan;
    
    if idx == nGroupsOut && faultResetSubsets
        typei = 2; % on-deck subset % next in the queue
    else
        typei = 1;
    end
    
    ssNeeded = cat(1,ssNeeded,ssGroupsPad);
    ssPsatList = cat(1,ssPsatList,ssPsatListPad);
    pSatNeeded = cat(1,pSatNeeded,pSatSsi);
    ssTypeNeeded = cat(1,ssTypeNeeded,typei*ones(size(ssGroupsi,1),1));
    ssWhichGroup = cat(1,ssWhichGroup,ssWhichGroupi);
end

%% Construct subsetfilter based on groups
for idx = 1:size(ssNeeded,1)
    
    ssNeededi = ssNeeded(idx,:);
    % check if this subset is already included
    if any(ssNeededi == obj.subsetMeasOut)
        % This measurement subset was found- keep it moving
        continue;
    end
    
    indNewSs = length(obj.subsetFilters)+1;
    
    % Subset type (primary or background)
    obj.subsetType(indNewSs,1)      = ssTypeNeeded(idx);
    % Measurement exclusion list
    obj.subsetMeasOut(indNewSs,:) = ssNeeded(idx,:);
    % Subset aggregate psat
    obj.subsetPsat(indNewSs,:)    = pSatNeeded(idx);
    % Subset ID
    obj.subsetId(indNewSs,1)        = obj.lastId+1;
    % Which group each thing comes from
    obj.subsetWhichGroup(indNewSs,:) = ssWhichGroup(idx,:);
    % Individual subset psats
    obj.subsetPsatList(indNewSs,:) = ssPsatList(idx,:);
    
    if obj.aivFilter.initialized == 1
        % the AIV was JUST initialized
        
        % Just put in an empty filter
        if indNewSs == 1
            % Need to just start this over for some reason with the
            % superclass pppfilter, maybe should fix at some point
            obj.subsetFilters = copy(obj.filter);
        else
            obj.subsetFilters(indNewSs,1) = copy(obj.filter);
        end
        
        % The measurements to be excluded from this initialization
        measExcludei = ssNeeded(idx,:);
        
        % Start the filter estimates, with the measurements being excluded
        obj.subsetFilters(indNewSs,1).initialize(obs,corrData...
                                                ,'measExclude',measExcludei...
                                                ,'tropsPreloadData',obj.aivFilter.tropsPreloadData);
        
    elseif obj.aivFilter.initialized == 2
        % the AIV has been running, so just use that one
        obj.subsetFilters(indNewSs,1) = copy(filterSpawn);
        
        
        %%%%% TO DO:
    end
    
    % fields: convert struct into cell
    excludeTypes = fields(obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThresh);
    
    % clean up excluded measurments in the corresponding subset
    for fdx = 1:length(excludeTypes)
        if strcmp(excludeTypes{fdx},'GNSS')
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThreshLarge.(excludeTypes{fdx}).Code = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThresh.(excludeTypes{fdx}).Code = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThreshLarge.(excludeTypes{fdx}).Carrier = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThresh.(excludeTypes{fdx}).Carrier = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThreshLarge.(excludeTypes{fdx}).Doppler = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThresh.(excludeTypes{fdx}).Doppler = Inf;
        else
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThreshLarge.(excludeTypes{fdx}) = Inf;
            obj.subsetFilters(indNewSs,1).PARAMS.measUse.excludeThresh.(excludeTypes{fdx}) = Inf;
        end
    end
    obj.lastId = obj.lastId+1;
    
end





end