function [latpp, lonpp, fpp] = ionoPiercePoint(latR, lonR, azS, elS)
% ionoPiercePoint
% DESCRIPTION:
%   Computation of the ionosphere piercing point (IPP).
% INPUT:
%   latR = receiver position latitude  [rad]
%   lonR = receiver position longitude [rad]
%   azS  = satellite azimuth [rad]
%   elS  = satellite elevation [rad]
% OUTPUT:
%   latpp = ionosphere piercing point latitude  [rad]
%   lonpp = ionosphere piercing point longitude [rad]
%   fpp   = slant factor (mapping function)
%
% See also: 
%
%
%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta 3
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Antonio Herrera Olmo, 2012
%  Contributors:     Eugenio Realini, 2013
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------
% 
%   Modified10/11/2021 by Fabian Rothmaier: vectorized the function

R  = 6378.1363; %Earth radius [km]
hI = 350;       %ionosphere thin shell height [km]

k = (R/(R+hI))*cos(elS);
phipp = (pi/2) - elS - asin(k);

%latitude of the ionosphere piercing point
latpp = asin(sin(latR).*cos(phipp) + cos(latR).*sin(phipp).*cos(azS));

%longitude of the ionosphere piercing point
case1 = ((latpp >  70*pi/180) & (tan(phipp).*cos(azS)      > tan((pi/2) - latR))) | ...
        ((latpp < -70*pi/180) & (tan(phipp).*cos(azS + pi) > tan((pi/2) + latR)));

asinTerm = asin(sin(phipp).*sin(azS./cos(latpp)));

lonpp = lonR + asinTerm;
lonpp(case1) = lonpp(case1) + pi - 2*asinTerm(case1);


%slant (obliquity) factor
fpp = (1 - k.^2).^(-1/2);

end
