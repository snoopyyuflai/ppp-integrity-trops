% clear variables;
close all;
%% load result
GPT2W_pos_filename = 'D:\PHD_Research\output\GPT2W_pos.mat';
GPT2W_pl_filename = 'D:\PHD_Research\output\GPT2W_pl.mat';
GPT2W_trops_filename = 'D:\PHD_Research\output\GPT2W_trops.mat';
GPT2W_epoch_filename = 'D:\PHD_Research\output\GPT2W_epoch.mat';
UNB3_pos_filename = 'D:\PHD_Research\output\UNB3_pos.mat';
UNB3_pl_filename = 'D:\PHD_Research\output\UNB3_pl.mat';
UNB3_trops_filename = 'D:\PHD_Research\output\UNB3_trops.mat';
NoDelay_pos_filename = 'D:\PHD_Research\output\NoDelay_pos.mat';
NoDelay_pl_filename = 'D:\PHD_Research\output\NoDelay_pl.mat';
NoDelay_trops_filename = 'D:\PHD_Research\output\NoDelay_trops.mat';

GPT2W_pos_all = load(GPT2W_pos_filename);
GPT2W_pl_all = load(GPT2W_pl_filename);
GPT2W_epoch_all = load(GPT2W_epoch_filename);
GPT2W_trops_all = load(GPT2W_trops_filename);
UNB3_pos_all = load(UNB3_pos_filename);
UNB3_pl_all = load(UNB3_pl_filename);
UNB3_trops_all = load(UNB3_trops_filename);
NoDelay_pos_all = load(NoDelay_pos_filename);
NoDelay_pl_all = load(NoDelay_pl_filename);
NoDelay_trops_all = load(NoDelay_trops_filename);

GPT2W_pos_all = GPT2W_pos_all.GPT2W_pos_all';
GPT2W_pl_all = GPT2W_pl_all.GPT2W_pl_all';
GPT2W_trops_all = GPT2W_trops_all.GPT2W_trops_all;
GPT2W_epoch_all = GPT2W_epoch_all.GPT2W_epoch_all';
UNB3_pos_all = UNB3_pos_all.UNB3_pos_all';
UNB3_pl_all = UNB3_pl_all.UNB3_pl_all';
UNB3_trops_all = UNB3_trops_all.UNB3_trops_all;
NoDelay_pos_all = NoDelay_pos_all.NoDelay_pos_all';
NoDelay_pl_all = NoDelay_pl_all.NoDelay_pl_all';
NoDelay_trops_all = NoDelay_trops_all.NoDelay_trops_all;

epochs = GPT2W_epoch_all;

truePosEcef = [-2706115.1823 -4278731.1983 3866392.5504];
truthFile = truePosEcef;

estim_size = size(GPT2W_pos_all);
%% convert xyz to enu
xyzpos_all = zeros(size(GPT2W_pos_all,1),3*size(GPT2W_pos_all,2));
enupos_all = zeros(size(GPT2W_pos_all,1),3*size(GPT2W_pos_all,2));
xyzpos_all(:,1:3) = GPT2W_pos_all;
xyzpos_all(:,4:6) = UNB3_pos_all;
xyzpos_all(:,7:9) = NoDelay_pos_all;

for j = 1:3
    xyz = xyzpos_all(:,3*(j-1)+1:3*j);
    llh0 = navsu.geo.xyz2llh(xyz(1,:));
    utmZone = navsu.thirdparty.findUtmZone(llh0(1),llh0(2));
    for idx = 1:size(xyz,1)
        [enu(idx,1),enu(idx,2),enu(idx,3)] = ...
            navsu.thirdparty.cart2utm(xyz(idx,1),xyz(idx,2),xyz(idx,3),utmZone);
    end
    enupos_all(:,3*(j-1)+1:3*j) = enu;
end
GPT2W_enupos_all   = enupos_all(:,1:3);
UNB3_enupos_all    = enupos_all(:,4:6);
NoDelay_enupos_all = enupos_all(:,7:9);
%% parse the truth
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
truthEnuInterp = nan(estim_size);
truthXyzInterp = nan(estim_size);
if isempty(truth.epochs)
    truthEnuInterp = repmat(truth.enuPos,estim_size(1),1);
    truthXyzInterp = repmat(truth.pos,estim_size(1),1);
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

%% plot absolute position error and protection level
close all
ylabels = {'NoDelay' 'UNB3' 'GPT2W'};
titles = {'All-in-View Absolute Error and Protection Level in ECEF-X [m]'...
          'All-in-View Absolute Error and Protection Level in ECEF-Y [m]'...
          'All-in-View Absolute Error and Protection Level in ECEF-Z [m]'};
GPT2W_pos_error = GPT2W_pos_all - truthXyzInterp;
UNB3_pos_error = UNB3_pos_all - truthXyzInterp;
NoDelay_pos_error = NoDelay_pos_all - truthXyzInterp;
yplot = zeros(size(NoDelay_pos_error,1), size(NoDelay_pos_error,2) ,3);
pl_yplot = zeros(size(GPT2W_pl_all,1), size(GPT2W_pl_all,2) ,3);
yplot(:,:,1) = NoDelay_pos_error;
yplot(:,:,2) = UNB3_pos_error;
yplot(:,:,3) = GPT2W_pos_error;
pl_yplot(:,:,1) = NoDelay_pl_all;
pl_yplot(:,:,2) = UNB3_pl_all;
pl_yplot(:,:,3) = GPT2W_pl_all;

tplot = (epochs-epochs(1))'/60; % minutes

for idx = 1:1:3 % xyz
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
    for idk = 1:1:3 % mdoels
        axes(ha(idk))
        plot(tplot,abs(yplot(:,idx,idk)))
        hold
        plot(tplot,pl_yplot(:,idx,idk),'k','linewidth',2)
        legend('Error', 'Protection Level')
        grid on
        ylabel(ylabels{idk})
        xlabel('Time [min]')
        if idk == 1
            title(titles{idx})
        end
        if idk < 3
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
end
%% plot total error
close all
yplot = zeros(size(NoDelay_pos_error,1), size(NoDelay_pos_error,2));
yplot(:,1) = vecnorm(NoDelay_pos_error,2,2);
yplot(:,2) = vecnorm(UNB3_pos_error,2,2);
yplot(:,3) = vecnorm(GPT2W_pos_error,2,2);
figure;
ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
for idx = 1:1:3 % xyz
    axes(ha(idx))
    plot(tplot,abs(yplot(:,idx)))
    grid on
    ylabel(ylabels{idx})
    xlabel('Time [min]')
    if idx == 1
        title('L^2 norm of Position Error')
    end
    if idx < 3
        xticklabels('')
    else
        xlabel('Time [min]')
    end
end
%% plot tropospheric delay contrast
yplot = zeros(size(NoDelay_trops_all,1), size(NoDelay_trops_all,2));
yplot(:,1) = NoDelay_trops_all(1:end);
yplot(:,2) = UNB3_trops_all(1:end);
yplot(:,3) = GPT2W_trops_all(1:end);
ylabels = {'NoDelay' 'UNB3' 'GPT2W'};
figure;
ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
for idx = 1:1:3 % xyz
    axes(ha(idx))
    plot(tplot,abs(yplot(:,idx)))
    grid on
    ylabel(ylabels{idx})
    xlabel('Time [min]')
    if idx == 1
        title('Tropospheric Delay [m]')
    end
    if idx < 3
        xticklabels('')
    else
        xlabel('Time [min]')
    end
end
%% plot error in ENU
ylabels = {'NoDelay' 'UNB3' 'GPT2W'};
titles = {'All-in-View Absolute Error in ECEF-E [m]'...
          'All-in-View Absolute Error in ECEF-N [m]'...
          'All-in-View Absolute Error in ECEF-U [m]'};
GPT2W_enupos_error = GPT2W_enupos_all - truthEnuInterp;
UNB3_enupos_error = UNB3_enupos_all - truthEnuInterp;
NoDelay_enupos_error = NoDelay_enupos_all - truthEnuInterp;
yplot = zeros(size(NoDelay_enupos_error,1), size(NoDelay_enupos_error,2) ,3);
pl_yplot = zeros(size(GPT2W_pl_all,1), size(GPT2W_pl_all,2) ,3);
yplot(:,:,1) = NoDelay_enupos_error;
yplot(:,:,2) = UNB3_enupos_error;
yplot(:,:,3) = GPT2W_enupos_error;

trops_yplot = zeros(size(NoDelay_trops_all,1), size(NoDelay_trops_all,2));
trops_yplot(:,1) = NoDelay_trops_all(1:end);
trops_yplot(:,2) = UNB3_trops_all(1:end);
trops_yplot(:,3) = GPT2W_trops_all(1:end);

tplot = (epochs-epochs(1))'/60; % minutes

for idx = 1:1:3 % enu
    figure;
    ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
    for idk = 1:1:3 % mdoels
        axes(ha(idk))
        plot(tplot,abs(yplot(:,idx,idk)))
        hold
        plot(tplot,abs(trops_yplot(:,idk)))
        legend('ENU Error', 'Trops delay')
        grid on
        ylabel(ylabels{idk})
        xlabel('Time [min]')
        if idk == 1
            title(titles{idx})
        end
        if idk < 3
            xticklabels('')
        else
            xlabel('Time [min]')
        end
    end
end
%% overlap result-with Nodelay
close all
yplot = zeros(size(NoDelay_pos_error,1), size(NoDelay_pos_error,2));
yplot(:,1) = vecnorm(NoDelay_pos_error,2,2);
yplot(:,2) = vecnorm(UNB3_pos_error,2,2);
yplot(:,3) = vecnorm(GPT2W_pos_error,2,2);
figure;
hold
% ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
for idx = 1:1:3 % xyz
%     axes(ha(idx))
    plot(tplot,abs(yplot(:,idx)))
    grid on
    xlabel('Time [min]')
    if idx == 1
        title('L^2 norm of Position Error')
    end
    if idx < 3
%         xticklabels('')
    else
        xlabel('Time [min]')
    end
end
legend('NoDelay', 'UNB3', 'GPT2W')
ylabel('Error [m]')
%% overlap result-with Nodelay
close all
yplot = zeros(size(NoDelay_pos_error,1), size(NoDelay_pos_error,2));
yplot(:,1) = vecnorm(NoDelay_pos_error,2,2);
yplot(:,2) = vecnorm(UNB3_pos_error,2,2);
yplot(:,3) = vecnorm(GPT2W_pos_error,2,2);
figure;
hold
% ha = navsu.thirdparty.tightSubplot(3,1,0.02,[0.1 0.1],[0.07 0.05]); 
for idx = 2:1:3 % xyz
%     axes(ha(idx))
    plot(tplot,abs(yplot(:,idx)))
    grid on
    xlabel('Time [min]')
    if idx == 2
        title('L^2 norm of Position Error')
    end
    if idx < 3
%         xticklabels('')
    else
        xlabel('Time [min]')
    end
end
legend('UNB3', 'GPT2W')
ylabel('Error [m]')

