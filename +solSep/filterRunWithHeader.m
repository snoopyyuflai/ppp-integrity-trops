% Function to run filters given a code injection header input 


function [runFilter, outData, truthFile, gnssMeas, wheelMeas,imuMeas] =  filterRunWithHeader(filenameHeader)

persistent obsGnssRaw corrData imuMeasRaw epochs wheelsRaw posRaw velRaw

run(filenameHeader)


% Check the 'inputs'
if ~exist('keepDataPersistent','var')
   keepDataPersistent = false; 
end

if ~exist('tStartPush','var')
    % Optional input from the header regarding pushing the start of the run
    tStartPush = 0;
end

if ~exist('tDuration','var')
    % Optional input from the header regarding the duration of the run
    tDuration = Inf; 
end

if ~exist('downsampleFactor','var')
   % Optional input from the header regarding downsampling the data
   downsampleFactor = 1;
end

if ~exist('orbMode','var')
    orbMode = 'PRECISE';
end
if ~exist('clkMode','var')
    clkMode = 'PRECISE';
end

if ~exist('plotOutputPos','var')
    plotOutputPos = 'REF';
end
%% Read observation file
if isempty(obsGnssRaw) %exist('obsGnssRaw','var')
    disp('Reading observation file')
    [~,~,extObs] = fileparts(filenameObs);
    
    switch extObs
        case {'.ASC','.GPS'}
            % Novatel ASCII log
            constellations = navsu.readfiles.initConstellation(constUse(1),constUse(2),constUse(3),constUse(4),constUse(5));
            gnssImuFull = navsu.internal.readSPANdata(filenameObs,'RANGEA',true,...
                'RAWIMUA',true,'constellations',constellations,...
                'correctHalfCycle',true);
            
            
            obsGnssRaw.meas      = gnssImuFull.RANGEA.obsData;
            obsGnssRaw.PRN       = constellations.PRN;
            obsGnssRaw.constInds = constellations.constInds;
            obsGnssRaw.epochs   = gnssImuFull.RANGEA.epochs;
            epochs = obsGnssRaw.epochs;
            
            tLock = gnssImuFull.RANGEA.tLock;
            tLock.epochs = obsGnssRaw.epochs;
            
            obsGnssRaw.tLock = tLock;
            
            if ~isempty(gnssImuFull.RAWIMUA.headerWeek)
                imuMeasRaw = gnssImuFull.RAWIMUA;
            else
                imuMeasRaw = gnssImuFull.RAWIMUSXA;
            end
            
            
            imuOrder = [-2 1 3];
            
            imuMeasRaw.acc  = imuMeasRaw.acc*accScale;
            imuMeasRaw.gyro = imuMeasRaw.gyro*gyroScale;
            
            imuMeasRaw.acc = [sign(imuOrder(1))*imuMeasRaw.acc(:,abs(imuOrder(1))) ...
                sign(imuOrder(2))*imuMeasRaw.acc(:,abs(imuOrder(2))) ...
                sign(imuOrder(3))*imuMeasRaw.acc(:,abs(imuOrder(3)))];
            
            imuMeasRaw.gyro = [sign(imuOrder(1))*imuMeasRaw.gyro(:,abs(imuOrder(1))) ...
                sign(imuOrder(2))*imuMeasRaw.gyro(:,abs(imuOrder(2))) ...
                sign(imuOrder(3))*imuMeasRaw.gyro(:,abs(imuOrder(3)))];
            
            imuMeasRaw.type = navsu.internal.MeasEnum.IMU;
            
            antMod = [];
            antOffset = [0 0 0]';
            
            % wheel speed information
            if isfield(gnssImuFull,'VEHICLEDATASA') && ~isempty(gnssImuFull.VEHICLEDATASA.headerWeek)
                wheelsRaw = gnssImuFull.VEHICLEDATASA;
            else
               wheelsRaw = []; 
            end
            
            

            'fdasf';
        case {'.19o','.18o', '.rnx','.18O','.19O','.obs','.17o','.17O'}
            % typical RINEX file
            
            [obsStruc, constellations, epochs, date, pos, interval, antoff, antmod,...
                rxmod] = navsu.readfiles.loadRinexObs(filenameObs,'constellations',...
                navsu.readfiles.initConstellation(constUse(1),constUse(2),constUse(3),constUse(4),constUse(5)));
            
            obsGnssRaw.meas      = obsStruc;
            obsGnssRaw.PRN       = constellations.PRN;
            obsGnssRaw.constInds = constellations.constInds;
            obsGnssRaw.epochs    = epochs;
            obsGnssRaw.tLock     = [];
            
            imuMeasRaw = [];
            
            wheelsRaw = [];
    end
end

%%
if isempty(corrData) % ~exist('corrData','var')
    disp('Loading corrections')
    % looking for most recent products
    jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) navsu.time.epochs2jd(max(epochs))]+0.5)-0.5;
    jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);
    
    [doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
    doyProd = floor(doyProd);
    
    corrData = navsu.svOrbitClock('configFile',configFile,'constUse',constUse);
    
    corrData.orbMode = orbMode;
    corrData.clkMode = clkMode;
    
    if strcmp(corrData.orbMode,'PRECISE')
        if exist('igsAc','var')
            % if a specific IGS analysis center was specified in the header
            % file, set it now
            corrData.settings.gpsEphCenter = igsAc;
            corrData.settings.gpsClkCenter = igsAc;
            corrData.settings.gloEphCenter = igsAc;
            corrData.settings.gloClkCenter = igsAc;
            corrData.settings.galEphCenter = igsAc;
            corrData.settings.galClkCenter = igsAc;
        end
        
        % load MGEX ephemeris data into our correction object
        corrData.initOrbitData(yearProd,doyProd);
        
        % load MGEX clock data into our correction object
        corrData.initClockData(yearProd,doyProd);
        
        % load antenna phase center data for satellites
        filenameAtx = 'igs14_sats_only.atx';
        
        corrData.initAtxData(filenameAtx);
        
    else
        corrData.initBroadcastData(yearProd,doyProd);
    end
    
    % ionospheric data
    corrData.initIonoData(yearProd,doyProd);
    
    % dcb products
    corrData.settings.dcbSource = 6;
    corrData.initDcb(yearProd,doyProd);
    
end

%% preprocess observations
[gnssMeas, dcbCorr0] = navsu.ppp.preprocessGnssObs(obsGnssRaw, corrData,...
    'epochStart',epochs(1)+tStartPush,'epochEnd',epochs(1)+tStartPush+tDuration,...
    'downsampleFactor',downsampleFactor);
if ~isempty(gnssMeas)
    gnssMeas = {gnssMeas};
end

imuMeas = navsu.ppp.preprocessImuMeas(imuMeasRaw,'epochStart',...
    epochs(1)+tStartPush,'epochEnd',epochs(1)+tStartPush+tDuration);
if ~isempty(imuMeas)
    imuMeas = {imuMeas};
end


wheelMeas = navsu.ppp.preprocessWheels(wheelsRaw,'epochStart',...
    epochs(1)+tStartPush,'epochEnd',epochs(1)+tStartPush+tDuration);
if ~isempty(wheelMeas)
    wheelMeas = {wheelMeas};
end

posMeas = navsu.ppp.preprocessPosMeas(posRaw,'epochStart',...
    epochs(1)+tStartPush,'epochEnd',epochs(1)+tStartPush+tDuration);
if ~isempty(posMeas)
    posMeas = {posMeas};
end

velMeas = navsu.ppp.preprocessVelMeas(velRaw,'epochStart',...
    epochs(1)+tStartPush,'epochEnd',epochs(1)+tStartPush+tDuration);
if ~isempty(velMeas)
    velMeas = {velMeas};
end



%% do the ppp lol
outData = navsu.ppp.runPpp(runFilter,[gnssMeas imuMeas wheelMeas posMeas velMeas],corrData);

%%
close all force;

% runFilter.plotOutput(outData,'truthFile',truthFile,'outputPos',plotOutputPos);

if exist('outputFile')
    save(outputFile,'outData','runFilter','truthFile')
end


%%
if ~keepDataPersistent
   clear corrData obsGnssRaw
end


end