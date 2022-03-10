function outData = runPppSs(ss,obs,corrData,varargin)

% Check what measurements are in here and sort out what's usable
obs = ss.checkMeas(obs);

% Sync all measurements 
[obsMap,epochs] = navsu.ppp.syncMeas(obs);

% Number of time steps
nEpochs = length(epochs);

% Initialize the progress bar
wb = navsu.internal.loadingBar(nEpochs);

% Data to save
outData = [];

%% Run the loop
for tdx = 1:nEpochs
    % Pull the measurement from the full list
    obsi = navsu.ppp.stripMeas(tdx,obs,obsMap);
    
    epochi = epochs(tdx);
    
    % Main filters
    if ss.initialized
        % do the update
        ss.update(epochi,obsi,corrData);
    else
       % need to initialize
        ss.initialize(obsi,corrData);
    end
    
    % Save some things for output
    outData = ss.saveState(outData,epochi,obsi);
    
    % Update the progress bar    
    wb.update(tdx);
end

% Close the progress bar
wb.close;


end

