function [measId,extraInputs] = update(obj, epoch, obs, corrData, varargin)

%%
% Manage the states in the filter :), with multi filters
% New satellites available, then need to add them to the state (along with
% new covariance etc...), and if the satellites are no longer in view, it
% has to be removed from the state. Also check and exclude cycle slip
measRemovedSlip = obj.manageStatesMulti(epoch, obs);

% Time update, EKF-Model Predict
obj.timeUpdate(epoch)

% Measurement update, Update GNSS data in next time step,measurement and
% EKF-Update
[measId, extraInputs] = obj.measUpdate(epoch, obs, corrData, measRemovedSlip, ...
                                       varargin{:});

% Make sure that the filter knows that it is running.
if (epoch - obj.epochLastGnssUpdate) < 10
    % if it's been too long, assume we're not initialized anymore
    obj.initialized = 2;
else
    obj.initialized = 0;
end

end
