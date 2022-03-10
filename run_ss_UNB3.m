% Sctipt to run a simple PPP algorithm using the navsu repo
% clear variables;
close all;
%% inputs

% RINEX v3 observation file
% filenameGnss = 'D:\PNT Data\Roof logs\swift-gnss-20200312-093212.sbp.obs';
filenameGnss = 'D:\PHD_Research\data\GNSS_OBS\swift-gnss-20200312-093212.obs';
% need a configutation file to set where to put downloaded products.  The
% default included is called default.ini
configFile = 'config.ini';
% Username and password file for NASA data/products download. See: 
% [1] https://cddis.nasa.gov/Data_and_Derived_Products/CDDIS_Archive_Access.html
% [2] https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html
netrcFile = 'D:\PHD_Research\Access\cddislogin.netrc'; % change to your directory
% cddis also requires a cookie file when using cURL, see [1]
cookieFile = 'D:\PHD_Research\Access\nasa_cddis_cookies.txt'; % change to your directory

% Truth position of the Xona roof antenna
truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];

% Three letter code indicating desired IGS analysis center
igsAc = 'GRG';

%
% this should be multi-constellation!
constUse = [1 1 1 0 0];  % GPS | GLO | GAL | BDS | QZSS



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

% epochStart = min(epochs)+200*60+60*10;
% downsampleFac = 30;

epochStart = min(epochs);
downsampleFac = 10; % e.g. total: 300 secs -> 300/10 = 30 time steps, but 
                    % sill 300 secs total
%% load navigation data
if ~exist('corrData','var')
    disp('Loading corrections')
    % looking for most recent products
    jdRange0 = floor([navsu.time.epochs2jd(min(epochs)) navsu.time.epochs2jd(max(epochs))]+0.5)-0.5;
    jdRangeProd = (min(jdRange0)-1):(max(jdRange0)+1);
    
    [doyProd,yearProd] = navsu.time.jd2doy(jdRangeProd);
    doyProd = floor(doyProd);
    
    corrData = navsu.svOrbitClock('configFile', configFile, ...
                                'constUse', constUse, ... % defaults to GPS only
                                'netrcFile', netrcFile, ... % not required on Mac!
                                'cookieFile', cookieFile); % not required on Mac!
    
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
    filenameAtx = 'D:\PHD_Research\data\ATX\igs14_sats_only.atx';
    
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
     'epochEnd',epochStart+60*30);
 
 % Build position measurement
% posMeas = navsu.ppp.buildPosMeas(gnssMeas.epochs,truePosEcef',100,1);
% 
% velMeas = navsu.ppp.buildVelMeas(gnssMeas.epochs,[0 0 0]',0.1,1);

 
%% Initialize the filter
filter = navsu.estimators.pppFilter;
% filter = navsu.estimators.leastSq;

filter.PARAMS.states.RX_DCB_GLO = false;
% filter.PARAMS.Q.POS = 1;
% filter.PARAMS.Q.VEL = 1;
filter.PARAMS.Q.VEL = 0;
% filter.PARAMS.Q.AMB = 0.002;

filter.PARAMS.measMask.f1 = [0 0 1]';
filter.PARAMS.measMask.f2 = [0 0 0]';
filter.PARAMS.measMask.f3 = [0 0 0]';

% filter.PARAMS.measUse.noVertVel = 1;
filter.PARAMS.measUse.noVertVel = 0;

% select tropspheric model
% filter.PARAMS.tropModel = 'GPT2W_1';
% filter.PARAMS.tropModel_FileName = 'D:\PHD_Research\data\GPT_GRID\gpt2_1w.grd.txt';

% Initialize the filter
solutionSep = solSep.solSepFilter(filter);

solutionSep.PARAMS.subsetGroupSize = 2;
solutionSep.PARAMS.nOutSubset = 2;



%% Estimate!
% Run the filter!
% corrData.orbMode = 'BROADCAST';
% corrData.clkMode = 'BROADCAST';

corrData.orbMode = 'PRECISE';
corrData.clkMode = 'PRECISE';

whatToRun = solutionSep;

% outData = solSep.runPppSs(whatToRun,{gnssMeas posMeas},corrData);
outData = solSep.runPppSs(whatToRun,{gnssMeas},corrData);

%%
close all;

whatToRun.plotOutput(outData,'truthFile',truePosEcef);

%%
% process trops
solutionSep.aivFilter.tropsAvg = solutionSep.aivFilter.tropsAvg(1:end-3)+cat(1,outData.tropo);
% save UNB3 outdata
UNB3_pos_filename = 'D:\PHD_Research\output\UNB3_pos.mat';
UNB3_pl_filename = 'D:\PHD_Research\output\UNB3_pl.mat';
UNB3_trops_filename = 'D:\PHD_Research\output\UNB3_trops.mat';
UNB3_pos_all = cat(2, outData.pos);
UNB3_pl_all  = cat(2, outData.pl);
UNB3_trops_all  = solutionSep.aivFilter.tropsAvg;
save(UNB3_pos_filename,'UNB3_pos_all');
save(UNB3_pl_filename,'UNB3_pl_all');
save(UNB3_trops_filename,'UNB3_trops_all');

