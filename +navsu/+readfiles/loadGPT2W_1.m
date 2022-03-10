function tropsGrid = loadGPT2W_1(filename)
fid = fopen(filename,'r');
% read first comment line
line = fgetl(fid);

% initialization
pgrid = zeros([64800, 5]);
Tgrid = zeros([64800, 5]);
Qgrid = zeros([64800, 5]);
dTgrid = zeros([64800, 5]);
u = zeros([64800, 1]);
Hs = zeros([64800, 1]);
ahgrid = zeros([64800, 5]);
awgrid = zeros([64800, 5]);
lagrid = zeros([64800, 5]);
Tmgrid = zeros([64800, 5]);

% loop over grid points
% 64800 for the 1 degree grid (GP)
for n = 1:64800
    
    % read data line
    line = fgetl(fid);
    vec = sscanf(line,'%f');
        
    % read mean values and amplitudes
    pgrid(n,1:5)  = vec(3:7);          % pressure in Pascal
    Tgrid(n,1:5)  = vec(8:12);         % temperature in Kelvin
    Qgrid(n,1:5)  = vec(13:17)./1000;  % specific humidity in kg/kg
    dTgrid(n,1:5) = vec(18:22)./1000;  % temperature lapse rate in Kelvin/m
    u(n)          = vec(23);           % geoid undulation in m
    Hs(n)         = vec(24);           % orthometric grid height in m
    ahgrid(n,1:5) = vec(25:29)./1000;  % hydrostatic mapping function coefficient, dimensionless
    awgrid(n,1:5) = vec(30:34)./1000;  % wet mapping function coefficient, dimensionless
	lagrid(n,1:5) = vec(35:39);    	   % water vapor decrease factor, dimensionless
	Tmgrid(n,1:5) = vec(40:44);    % mean temperature in Kelvin
    
end
fclose (fid);
tropsGrid.pgrid = pgrid;
tropsGrid.Tgrid = Tgrid;
tropsGrid.Qgrid = Qgrid;
tropsGrid.dTgrid = dTgrid;
tropsGrid.u = u;
tropsGrid.Hs = Hs;
tropsGrid.ahgrid = ahgrid;
tropsGrid.awgrid = awgrid;
tropsGrid.lagrid = lagrid;
tropsGrid.Tmgrid = Tmgrid;
end

