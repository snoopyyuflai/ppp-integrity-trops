% Header example

% To run solSep.filterRunWithHeader, you need to pass in a header, which is
% really just a code injection.  Everything below should just be straight
% up MATLAB code that sets a bunch of settings for your filter run.

% Required outputs
%   filenameObs - the filename and path of the observation data
%   constUse    - [1x5] constellation selection boolean [1 1 0 0 0] selects
%                 GPS and GLONASS
%   filter      - a filter to run- navsu.estimators.AbstractNavFilter

filenameObs = 'E:\Box Sync\Novatel Data\7500_OpenSky\Nov_1500_01Mar18_OpenSky_IMU.GPS';

constUse = [1 1 0 0 0];  % GPS | GLO | GAL | BDS | QZSS

filter = navsu.estimators.inertialPppFilter;

% Optional outputs - output an empty variable if you don't want to specify.
%   configFile  - should probably include one if you don't want your
%                 precise data to be downloaded to the current directory   
%   truthFile   - truth file for comparison (if you want)
%   igsAc       - three letter IGS code for the analysis center to download
%                 the precise data from
%   outputFile  - name for an output file.  Default will just have a
%                 similar name to the filenameObs and be in the same
%                 folder. 
%   tStartPush  - exclude this many seconds of measurements at the
%                 beginning of the run
%   tDuration   - run this many seconds 
%   keepDataPersistent - don't delete the observation and correction data
%                 from the persistent variable- useful for consecutive runs

configFile = 'config.ini';

truthFile = 'E:\Box Sync\Novatel Data\7500_OpenSky\Truth_IGS08_OpenSky.txt';

igsAc = 'grm';


% Really optional inputs
% You can just make whatever changes to filter parameters.
% Example: 
%   filter = navsu.estimators.inertialPppFilter;
%   filter.PARAMS.states.RX_DCB_GLO = false;
%   filter.PARAMS.measMask.f1 = [1 1 1]';
%   filter.PARAMS.measMask.f2 = [0 0 0]';

filter.PARAMS.IMU_ARM = [0.328 0.327 1.310]';

