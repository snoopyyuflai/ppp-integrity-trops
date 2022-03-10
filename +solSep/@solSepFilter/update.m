function  extraInputs = update(obj,epoch,obs,corrData)

aivFilter0 = copy(obj.aivFilter);

% Update the all in view filter
[measUsedAiv,extraInputs] = obj.aivFilter.update(epoch,obs,corrData);

% Subset come from AIV in each time stamp, so in each estimation,
% the subset need to be reinitialized every time

measRemovedAiv = [];
if ~isempty(measUsedAiv)
    % only including measurements that were rejected- no cycle slips or low
    % elevation satellites get passed to the subset filters.  those will be
    % handled by the subset filters. 
    if ~isempty(obj.aivFilter.measRemoved)
        measRemovedAiv = obj.aivFilter.measRemoved.id(obj.aivFilter.measRemoved.reason == 2);
    end
    manageSubsets(obj,measUsedAiv,obs,corrData,aivFilter0); 
end
% Update the subset filters
for idx = 1:length(obj.subsetFilters) 
    
    % If subsets were just initialized, they may not need to be updated. 
    
    % If the time of the last measurement update (this could be also during
    % a least squares initialization), then there is no need to update the
    % filter at this time step. 
    
  [measUsedSsi] = obj.subsetFilters(idx).update(epoch,obs,corrData,'measExclude',[obj.subsetMeasOut(idx,:)'; ...
       measRemovedAiv],'extraInputs',extraInputs);
       
end


if obj.aivFilter.initialized == 0
    % the aill in view filter was de-initialized. need to start the whole
    % thing over. 
    obj.initialized = 0;
end


if ~isempty(measUsedAiv)
    % Update the protection levels
    computePl(obj);
else
   obj.pl = nan(3,1);
end

end