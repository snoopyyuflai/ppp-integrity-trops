function [trop0,m,tropDataSave] = tropDelay(el, az, h, lat, lon, doy, ...
                                            PARAMS, tropData, INDS_STATIONS, ...
                                            epochi, tropsPreloadData)

% Elevation and azimuth are in degrees
% latitude and lontitude are in degrees
% lat = ellipsoidal latitude
% lon = ellipsoidal lontitude

switch PARAMS.tropModel
    case 'UNB3'
        [trop0, m, tropDataSave] = navsu.ppp.models.tropoErrCorrUnb3(el, h, lat, doy);
%     otherwise
%         error(['"%s" Tropo error model currently not implemented.\n', ...
%                'Can only do "UNB3".'], PARAMS.tropModel);
%     case 'IGS'
%         [trop0,m,tropDataSave] = tropo_error_correction_IGS(el,az,h,lat,lon,doy,PARAMS,tropData,...
%             INDS_STATIONS,epochi);
%         
%     case 'SAAS'
%         % saastomoinen
%         [trop0] = saastamoinen_model_SU(lat, lon, h, el);
%         m = zeros(size(trop0));
%         tropDataSave.ddry = zeros(size(trop0));
%         tropDataSave.gmfwSave  = zeros(size(trop0));
%         
%     case 'GO'
%         [trop0,m,tropDataSave] = tropo_error_correction_go(el,h,lat,doy);
%         
%         
%     case 'GMF' 
%         mjd = jd2mjd(epochs2jd(epochi));
%         
%         [trop0,m,tropDataSave] = tropo_error_correction_gmf(el,h,lat,lon,mjd,epochi);

    case 'GPT2W_1' 
        mjd = navsu.time.jd2mjd(navsu.time.epochs2jd(epochi));
        [trop0,m,tropDataSave] = navsu.ppp.models.tropoErrCorrGPT2W_1(el,h,lat,lon,mjd,epochi,tropsPreloadData);

    case 'NoDelay' 
        trop0 = zeros(size(el));
        m     = zeros(size(el));
        tropDataSave.gmfwSave   = zeros(size(el));
        tropDataSave.dall       = zeros(size(el));
        tropDataSave.ddry       = zeros(size(el));
        tropDataSave.dwet       = zeros(size(el));
  
    otherwise
        error(['"%s" Tropo error model currently not implemented.\n', ...
               'Can only do "UNB3".'], PARAMS.tropModel);


end

end