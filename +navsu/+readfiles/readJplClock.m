function [epochs, clk, clk_sig] = readJplClock(filename, dataLines)
% readJplClock
% DESCRIPTION:
%   Parser for JPL real-time clock products. 
%
% INPUT:
%   filename    - the name of the file (please include the path as well)
%
% OUTPUT:
%   epochs      - Nx1 vector of GPS epochs of clock outputs
%   clk         - Nx32 matrix of GPS clock biases [s]
%   clk_sig     - Nx32 matrix of GPS clock bias sigmas [s]
%
% See also: navsu.readfiles.readRinexClock

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["jpl_epoch", "val1", "val2", "val3", "satAndType"];
opts.VariableTypes = ["double", "double", "double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, "satAndType", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "satAndType", "EmptyFieldRule", "auto");

% Import the data
data = readtable(filename, opts);

%% Make it into something useful
indsClkBias = navsu.strFindCell(data.satAndType,'Clk');

% 
epochs = navsu.time.jpl2epochs(data.jpl_epoch(indsClkBias));

constellation = ones(size(epochs));


c = navsu.constants('c');
bias = data.val2(indsClkBias)/c;

rate = data.val3(indsClkBias)/c;

% need to pull PRN
svnStr = cell2mat(cellfun(@(x) x(15:16),data.satAndType(indsClkBias),'UniformOutput',false));
svns = str2num(svnStr);

prns = navsu.svprn.svn2prn(svns,epochs,constellation(1));

epochsUn = unique(epochs);
prnUn = 1:32;

% Make it into a nice 32xN matrix
clk     = nan(length(prnUn),length(epochsUn));
clk_sig = nan(length(prnUn),length(epochsUn));

for idx = 1:length(prnUn)
    indsPrni = find(prns == prnUn(idx));
    
    epochsi = epochs(indsPrni);
    biasi = bias(indsPrni);
    ratei = rate(indsPrni);
    
    [~,ixb] = ismember(epochsi,epochsUn);
    
    clk(idx,ixb) = biasi;
    clk_sig(idx,ixb) = ratei;
end

epochs = epochsUn;

end


























