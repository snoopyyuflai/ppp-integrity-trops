function outData = saveState(obj,outData,epoch,obs)
% Create a structure with all the information that you want to save :)


outState = [];

% If this an IMU update, don't save this value
if length(obs) == 1 && ~isempty(navsu.ppp.pullMeasFromList(obs,navsu.internal.MeasEnum.IMU))
    return;
end

% outState.epoch = epoch;
% outState.aivPos  = obj.aivFilter.pos;
% outState.pl      = obj.pl;
% outState.aivCovPos = obj.aivFilter.cov(obj.aivFilter.INDS_STATE.POS,obj.aivFilter.INDS_STATE.POS);
% outState.resids = obj.aivFilter.resids;
% outState.residsInfo = [];
% outState.measRemoved = obj.aivFilter.measRemoved;

outState = obj.aivFilter.saveOutStatePpp(outData,epoch,obs);
outState.pl      = obj.pl;

outState.fullSet = copy(obj);


% Save a bunch of stuff for the subset filters
ssCov = cat(3,obj.subsetFilters.cov);
covPosSs = ssCov(obj.aivFilter.INDS_STATE.POS,obj.aivFilter.INDS_STATE.POS,:);
outState.ssCovPos = covPosSs;

posSs = cat(2,obj.subsetFilters.pos);

outState.ssPos = posSs;
outState.ssEpochs = repmat(epoch,1,size(posSs,2));
outState.ssId  = obj.subsetId';
outState.ssMeasId = obj.subsetMeasOut';


% if isempty(outData) || isempty([outData(:).residsInfo])
%     gnssMeas = [];
%     for idx = 1:length(obs)
%         obsi = obs{idx};
%         switch obsi.type
%             case navsu.internal.MeasEnum.GNSS
%                 gnssMeas = obsi;
%         end
%     end
%     
%     if ~isempty(gnssMeas)
%         
%         % this is the first one- include some more information
%         outState.residsInfo.rangeInfo = gnssMeas.range;
%         outState.residsInfo.rangeInfo.obs = [];
%         outState.residsInfo.rangeInfo.lockTime = [];
%         
%         outState.residsInfo.dopplerInfo = gnssMeas.doppler;
%         outState.residsInfo.dopplerInfo.obs = [];
%     end
% end

outData = [outData; outState];

end