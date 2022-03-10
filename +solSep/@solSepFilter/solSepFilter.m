classdef solSepFilter < matlab.mixin.Copyable
    
    properties %(SetAccess = immutable)
        filter
        
        %     end
        %
        %     properties  %(SetAccess = public)
        % The all in view filter
        aivFilter %navsu.estimators.AbstractNavFilter
        
        % The primary subset filters
        subsetFilters  navsu.estimators.pppFilter
        
        % Excluded measurements for each subset. This is a cell array of
        % lists of navsu.internal.MeasID
        subsetMeasOut navsu.internal.MeasID
        
        % Probability of failure for each aggregate subset
        subsetPsat
        
        % ID number just to keep track of each subset
        subsetId
        
        % The last used ID number
        lastId = 0
        
        % Indicator of which group these subsets are from
        subsetWhichGroup
        
        % psat for each measurement in each subset
        subsetPsatList
        
        % type of subset, 1 = primary, 2 = background
        subsetType
        
        % Groupings to construct each subset
        groups navsu.internal.MeasID
        
        % Psat for each element in the subset groupings
        groupsPsat
        
        % Parameters related to running the subsets
        PARAMS
        
        % Protection level information
        pl = nan(3,1);
        
        % Flag indicating if this has been initialized
        initialized = false
        
    end
    
    methods
        
        function obj = solSepFilter(filter)
            % Constructor- you MUST supply an initial copy of the filter
            obj.filter = filter;
            
            obj.PARAMS = initParams(obj);
            
        end
    end
    
    
    methods
        
        PARAMS = initParams(obj);
        
        % Initialize the aiv and subsets
        initialize(obj,obs,corrData)
        
        % Update the filters
        update(obj,epoch,obs,corrData)
        
        % Manage subsets given a list of measurement IDs-
        % Can require spawning new subsets
        manageSubsets(obj,measId,obs,corrData,filterSpawn)
        
        % Compute the protection level
        computePl(obj);
        
        % Save some data
        outData = saveState(obj,outData,epoch,obs);
        
        % Check the set of measurements given
        obs = checkMeas(obj,obs0)
    end
    
    methods 
        % Plot the data
        plotOutput(obj,outputs,varargin)
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end