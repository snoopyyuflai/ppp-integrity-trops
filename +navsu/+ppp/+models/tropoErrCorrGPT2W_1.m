function [tropo_corr, Mw, tropData] = tropoErrCorrGPT2W_1(el,h,lat,lon,mjd,epochi,tropsPreloadData)
%%GPT2W model for tropospheric delay, with 1x1 degree accuracy
dmjd = mjd;
nstat = size(el,1);
lat = lat*pi/180;
lon = lon*pi/180;
dlat = repelem(lat, nstat)';
dlon = repelem(lon, nstat)';
hell = repelem(h, nstat)';
zd = 90 - el;
zd = zd*pi/180; % zenith distance
it = 1;
[p,T,dT,Tm,e,ah,aw,la,undu] = navsu.ppp.models.gpt2_1w(dmjd,dlat,dlon,hell,nstat,it,tropsPreloadData);
p = p';
Tm = Tm';
e = e';
ah = ah';
aw = aw';
la = la';


% compute dry delay
zhd = navsu.ppp.models.saasthyd(p,dlat,hell);
% compute wet delay
zwd = navsu.ppp.models.asknewet(e,Tm,la);
% compute mapping functions
[vmf1h,vmf1w] = navsu.ppp.models.vmf1(ah,aw,dmjd,dlat,zd);
% compute total delay
ztd = vmf1h.*zhd + vmf1w.*zwd;
% save data
tropo_corr = ztd;
Mw = vmf1w;
tropData.gmfwSave   = Mw;
tropData.dall       = zhd + zwd;
tropData.ddry       = zhd;
tropData.dwet       = zwd;
% remove bad data...
tropo_corr(tropo_corr > 1e3) = 0;
if ~isreal(tropo_corr)
    tropo_corr = zeros(size(tropo_corr));
end

end