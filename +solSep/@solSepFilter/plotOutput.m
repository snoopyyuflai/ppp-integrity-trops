function plotOutput(obj,outputs,varargin)
%
% No all in view solution for now
% obj.aivFilter.plotOutputPpp(outputs,varargin{:})
% 'fdasf';

%%
p = inputParser;

p.addParameter('truePosEcef',[]);
p.addParameter('truthFile',[]);
p.addParameter('outputPos','REF');  % reference position (IMU) or antenna phase center

% parse the results
parse(p, varargin{:});
res = p.Results;
truePosEcef          = res.truePosEcef;          % truth position to compare to
truthFile            = res.truthFile;      %
outputPos = res.outputPos;

% Plot the position and clock bias in ENU
%% Pull out saved values
xyz = cat(2,outputs.pos)';
xyzCov = cat(3,outputs.covPos);
epochs = cat(2,outputs.epoch)';
b = zeros(size(epochs));
pl = cat(2,outputs.pl);

%% Subset positions and covariances
xyzSsList = cat(2,outputs.ssPos);
xyzCovSsList = cat(3,outputs.ssCovPos);
idSsList = cat(2,outputs.ssId);
epochSsList = cat(2,outputs.ssEpochs);

measOutList = cat(2,outputs.ssMeasId);

ids = unique(idSsList);
measOut = [];
nEpochs = length(epochs);

xyzSs = nan(nEpochs,3,length(ids));
xyzCovSs = nan(3,3,nEpochs,length(ids));

for idx = 1:length(ids)
    indsId = find(idSsList == ids(idx));
    
    epochsi = epochSsList(indsId)';
    
    [~,ib] = ismember(epochsi,epochs);
    
    xyzSs(ib,:,idx) = xyzSsList(:,indsId)';
    xyzCovSs(:,:,ib,idx) = xyzCovSsList(:,:,indsId);
    
    measOut = cat(2,measOut,measOutList(:,indsId(end)));
    
end


%%
llh0 = navsu.geo.xyz2llh(xyz(1,:));


llh = navsu.geo.xyz2llh(xyz(:,:));

utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));

% Convert estimated position to enu
enu    = nan(size(xyz));
enuStd = nan(size(xyz));

for idx = 1:size(xyz,1)
    [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
        navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
    
    [~,RxyzEnu] = navsu.geo.xyz2enu([0 0 0],llh(idx,1)*pi/180,llh(idx,2)*pi/180);
    
    covEnui = RxyzEnu'*squeeze(xyzCov(:,:,idx))*RxyzEnu;
    enuStd(idx,:) = sqrt(diag(covEnui));
end
estim.pos = xyz;
estim.epochs = epochs;
estim.enuPos = enu;
estim.enuStd = enuStd;

% Convert subsets to ENU
enuSs    = nan(size(xyzSs));
enuStdSs = nan(size(xyzSs));

for idx = 1:size(xyzSs,1)
    for jdx = 1:size(xyzSs,3)
        [enuSs(idx,1,jdx),enuSs(idx,2,jdx),enuSs(idx,3,jdx)] = ...
            navsu.thirdparty.cart2utm(xyzSs(idx,1,jdx),xyzSs(idx,2,jdx),xyzSs(idx,3,jdx),utmZone);
        
        [~,RxyzEnu] = navsu.geo.xyz2enu([0 0 0],llh(idx,1)*pi/180,llh(idx,2)*pi/180);
        
        covEnui = RxyzEnu'*squeeze(xyzCovSs(:,:,idx,jdx))*RxyzEnu;
        enuStdSs(idx,:,jdx) = sqrt(diag(covEnui));
    end
end


if ~isempty(truthFile)
    % Parse the truth data
    [~,posTruth, epochsTruth,cov,~,~,velEnu,stdEnu] = navsu.internal.parsePosSolFile(truthFile);
    
    truth.pos = posTruth;
    truth.epochs = epochsTruth;
    truth.enuPos = nan(size(posTruth));
    truth.enuStd = stdEnu;
    
    llh0 = navsu.geo.xyz2llh(posTruth(1,:));
    % Convert truth data to ENU
    utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));
    for idx = 1:size(posTruth,1)
        [truth.enuPos(idx,1),truth.enuPos(idx,2),truth.enuPos(idx,3)] = ...
            navsu.thirdparty.cart2utm(posTruth(idx,1),posTruth(idx,2),posTruth(idx,3),utmZone);
    end
    
    % Interpolate truth data to match the estimated data :)
    truthEnuInterp = nan(size(estim.enuPos));
    truthXyzInterp = nan(size(estim.enuPos));
    if isempty(truth.epochs)
        truthEnuInterp = repmat(truth.enuPos,size(estim.enuPos,1),1);
        truthXyzInterp = repmat(truth.pos,size(estim.enuPos,1),1);
    else
        % Pad truth gaps with NaNs
        indsGap = find(diff(truth.epochs) > 1);
        truth.enuPos(indsGap,:) = nan(length(indsGap),3);
        truth.pos(indsGap,:) = nan(length(indsGap),3);
        
        interpMethod = 'linear';
        truthEnuInterp(:,1) = interp1(truth.epochs,truth.enuPos(:,1),estim.epochs,interpMethod);
        truthEnuInterp(:,2) = interp1(truth.epochs,truth.enuPos(:,2),estim.epochs,interpMethod);
        truthEnuInterp(:,3) = interp1(truth.epochs,truth.enuPos(:,3),estim.epochs,interpMethod);
        
        truthXyzInterp(:,1) = interp1(truth.epochs,truth.pos(:,1),estim.epochs,interpMethod);
        truthXyzInterp(:,2) = interp1(truth.epochs,truth.pos(:,2),estim.epochs,interpMethod);
        truthXyzInterp(:,3) = interp1(truth.epochs,truth.pos(:,3),estim.epochs,interpMethod);
    end
    
    truth.enuPosInterp = truthEnuInterp;
    truth.xyzPosInterp = truthXyzInterp;
else
    % There is no truth.  But this you must learn for yourself.
    truth = [];
end


%% plot

if ~isempty(truth)
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]);
    
    yplot = [estim.enuPos-truth.enuPosInterp]';
    yplotStd = pl;
    tplot = (epochs-epochs(1))'/60; % minutes
    ylabels = {'East [m]' 'North [m]' 'Up [m]' 'Clock [m]'};
    for idx = 1:3
        axes(ha(idx))
        
        plot(tplot,abs(yplot(idx,:)))
        hold on;
        plot(tplot,yplotStd(idx,:),'k','linewidth',2)
        %         plot(tplot,yplotStd(idx,:),'k')
        
        ylabel(ylabels{idx})
        grid on
        legend('Absolute Error (m)', 'Protection Level (m)')
        if idx == 1
            title('AIV Error and Protection Level')
        end
        if idx < 3
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
    
end

%% plot all in view and subset covariances


figure;
ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]);

prnsOut = nan(size(measOut));
constOut = nan(size(measOut));
for rdx = 1:size(measOut,1)
    for cdx = 1:size(measOut,2)
        if isa(measOut(rdx,cdx),'navsu.internal.MeasIdGnss')
           prnsOut(rdx,cdx) = measOut(rdx,cdx).prn;
           constOut(rdx,cdx) = measOut(rdx,cdx).const;
        end
    end
end


yplot1 = estim.enuStd;
%     yplotStd = [estim.enuStd b*navsu.constants.c]';
yplot2 = enuStdSs;%./estim.enuStd;
tplot = (epochs-epochs(1))/60; % minutes
ylabels = {'East [m]' 'North [m]' 'Up [m]' };
for idx = 1:3
    axes(ha(idx))
    
%     plot(tplot,squeeze(yplot2(:,idx,:))','linewidth',0.5)
%     hold on;
%     plot(tplot,yplot1(:,idx),'k','linewidth',2)
    
    h(1) = plot(tplot,yplot1(:,idx),'k','linewidth',2);
    hold on;
    plot(tplot,squeeze(yplot2(:,idx,:))','linewidth',0.5)
    ylabel(ylabels{idx})
    grid on
    legend('All-In-View (m)')
    uistack(h(1),'top');
    
    if idx == 1
        title('Subset and AIV covariances')
    end
    
    if idx < 3
        xticklabels('')
    else
        xlabel('Time [min]')
    end
end

%% plot all in view and subset position errors
if ~isempty(truth)
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]);
    yplot1 = [estim.enuPos-truth.enuPosInterp];
    
    % yplot1 = estim.enuStd;
    yplot2 = enuSs-truth.enuPosInterp;
    tplot = (epochs-epochs(1))/60; % minutes
    ylabels = {'East [m]' 'North [m]' 'Up [m]' };
    
    
    for idx = 1:3
        axes(ha(idx))
        
        h(1) = plot(tplot,yplot1(:,idx),'k','linewidth',2);
        hold on;
        
        s = plot(tplot,squeeze(yplot2(:,idx,:))','linewidth',0.5);
        
        for jdx = 1:length(s)
            s(jdx).DataTipTemplate.DataTipRows(1).Label = 't';
            s(jdx).DataTipTemplate.DataTipRows(2).Label = 'resid';
            
            
            for mdx = 1:size(prnsOut,1)
                row = dataTipTextRow('PRN',prnsOut(mdx,jdx)*ones(size(yplot2,1),1));
                s(jdx).DataTipTemplate.DataTipRows(2+2*(mdx-1)+1) = row;
                 row = dataTipTextRow('const',constOut(mdx,jdx)*ones(size(yplot2,1),1));
                s(jdx).DataTipTemplate.DataTipRows(2+2*(mdx-1)+2) = row;
            end
%             row = dataTipTextRow('PRN',double(indsUnPh(idx,1)).*ones(nEpochs,1));
%             s(idx).DataTipTemplate.DataTipRows(3) = row;
%             row = dataTipTextRow('const',double(indsUnPh(idx,2)).*ones(nEpochs,1));
%             s(idx).DataTipTemplate.DataTipRows(4) = row;
%             row = dataTipTextRow('sig',double(indsUnPh(idx,3)).*ones(nEpochs,1));
%             s(idx).DataTipTemplate.DataTipRows(5) = row;
        end

        legend('All-In-View (m)')
        uistack(h(1),'top'); % bring AIV to the front
        ylabel(ylabels{idx})
        grid on
        
        if idx == 1
            title('AIV and subset position errors')
        end
        
        if idx < 3
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
end

end



























