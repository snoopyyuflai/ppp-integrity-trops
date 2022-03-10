% Sctipt to run a simple PPP algorithm using the navsu repo

%% inputs
% RINEX v3 observation file
filenameGnss = 'D:\PNT Data\Roof logs\swift-gnss-20200312-093212.sbp.obs';
% filenameGnss = 'C:\Users\kazgu\Desktop\Stanford\swift-gnss-20200312-093212.sbp.obs';

% need a configutation file to set where to put downloaded products.  The
% default included is called default.ini
configFile = 'config.ini';

% Truth position of the Xona roof antenna
truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];

% Three letter code indicating desired IGS analysis center
igsAc = 'GRG';

%
% this should be multi-constellation!
constUse = [1 1 1 0 0];  % GPS | GLO | GAL | BDS | QZSS

filter = navsu.estimators.pppFilter;
% filter = navsu.estimators.leastSq;

filter.PARAMS.states.RX_DCB_GLO = false;
filter.PARAMS.Q.POS = 1;
filter.PARAMS.Q.VEL = 1;
filter.PARAMS.Q.AMB = 0.002;

filter.PARAMS.measMask.f1 = [0 0 1]';
filter.PARAMS.measMask.f2 = [0 0 0]';
filter.PARAMS.measMask.f3 = [0 0 0]';

filter.PARAMS.measUse.noVertVel = 1;

% Initialize the filter
solutionSep = solSep.solSepFilter(filter);

solutionSep.PARAMS.nOutSubset= 1;

%% Read observation file
if ~exist('obsStruc','var')   
    disp('Reading observation file')
    [obsStruc, constellations, epochs, date, pos, interval, antoff, antmod,...
        rxmod] = navsu.readfiles.loadRinexObs(filenameGnss,'constellations',...
        navsu.readfiles.initConstellation(constUse(1),constUse(2),constUse(3),constUse(4),constUse(5)));
    
    obsGnssRaw.meas      = obsStruc;
    obsGnssRaw.PRN       = constellations.PRN;
    obsGnssRaw.constInds = constellations.constInds;
    obsGnssRaw.epochs    = epochs;
    obsGnssRaw.tLock     = [];
end

epochStart = min(epochs)+200*60+60*10;
downsampleFac = 30;

%%
if ~exist('corrData','var')
    disp('Loading corrections')
    % looking for most recent products
    jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) navsu.time.epochs2jd(max(epochs))]+0.5)-0.5;
    jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);
    
    [doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
    doyProd = floor(doyProd);
    
    corrData = navsu.svOrbitClock('configFile',configFile,'constUse',constUse);
    
    corrData.settings.gpsEphCenter = igsAc;
    corrData.settings.gpsClkCenter = igsAc;
    corrData.settings.gloEphCenter = igsAc;
    corrData.settings.gloClkCenter = igsAc;
    corrData.settings.galEphCenter = igsAc;
    corrData.settings.galClkCenter = igsAc;
    
    % load MGEX ephemeris data into our correction object
    corrData.initOrbitData(yearProd,doyProd);
    
    % load MGEX clock data into our correction object
    corrData.initClockData(yearProd,doyProd);
    
    % Load broadcast data
%     corrData.initBroadcastData(yearProd,doyProd);
    
    % load antenna phase center data for satellites
    filenameAtx = 'igs14_sats_only.atx';
    
    corrData.initAtxData(filenameAtx);
    
    % dcb products
    corrData.settings.dcbSource = 6;
    corrData.initDcb(yearProd,doyProd);
    
    % ionospheric data
    corrData.initIonoData(yearProd,doyProd);
end

%% preprocess observations
[gnssMeas, dcbCorr0] = navsu.ppp.preprocessGnssObs(obsGnssRaw,...
     corrData,'downsampleFac',downsampleFac,'epochStart',epochStart,...
     'epochEnd',epochStart+60*10);
 
 % Build position measurement
posMeas = navsu.ppp.buildPosMeas(gnssMeas.epochs,truePosEcef',100,1);

% velMeas = navsu.ppp.buildVelMeas(gnssMeas.epochs,[0 0 0]',0.1,1);

 
%% Estimate!
corrData.orbMode = 'BROADCAST';
corrData.clkMode = 'BROADCAST';

corrData.orbMode = 'PRECISE';
corrData.clkMode = 'PRECISE';

whatToRun = solutionSep;

outData = solSep.runPppSs(whatToRun,{gnssMeas posMeas},corrData);

%%
close all;

whatToRun.plotOutput(outData,'truthFile',truePosEcef);




