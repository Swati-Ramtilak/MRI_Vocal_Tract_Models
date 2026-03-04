function [f1, err] = fit_zte_hist_model(img, SI_thres, Nminpts, Nbins, binRange, N_bkwidth_left, N_bkwidth_right, gaussEqn, noisePeak, gaussC, plotfit,savepath)
% Default values for optional arguments
if nargin < 8 || isempty(gaussEqn)
    gaussEqn = 'a*exp(-((x-b)/c)^2)';
end
if nargin < 9 || isempty(noisePeak)
    noisePeak = 0;
end   
if nargin < 10 || isempty(gaussC)
    gaussC = 5;
end
if nargin < 11 || isempty(plotfit)
    plotfit = 0;
end
if isempty(binRange) || length(binRange) ~= 2 || any(~isfinite(binRange))
    binRange = [min(img(:)), max(img(:))];
end
% Histogram and smoothing
[N,edges] = histcounts(img, Nbins,'BinLimits', binRange);
N(1:10) = 0;
Nidx = find(N > 0);
Ntmp = smooth(smooth(smooth(N(Nidx))));
edgestmp = edges(Nidx);
sg_pk_idx = 1;
iter = 0;
scale = 0.8;                                                               
% Main loop to adjust SI_thres if necessary
while sg_pk_idx == 1
    if iter > 0
        SI_thres = SI_thres * scale;                                       %To decrease threshold to search for peaks at lower intensity levels
    end
    if noisePeak == 0
        idx = find(edgestmp(1:end-1) > SI_thres);
    else
        idx = find(edgestmp(1:end-1) < SI_thres);
    end
    Ntmp1 = Ntmp(idx);
    edges1 = edgestmp(idx);   
    sg_pk_N = max(Ntmp1);
    sg_pk_idx = find(Ntmp1 == sg_pk_N, 1, 'last');
    iter = iter + 1;
end
% Define the region of interest around the peak
idx_vct = max([2 sg_pk_idx - N_bkwidth_left]):min(sg_pk_idx + N_bkwidth_right, length(Ntmp1));
N_sg = Ntmp1(idx_vct);
edges_sg = edges1(idx_vct);
sg_X_idx = find(N_sg > 0);
sg_X = edges_sg(sg_X_idx);
sg_Y = N_sg(sg_X_idx)';                                                    %Non linear least square fitting on gaussian curve to histogram region identified
if length(sg_X) > Nminpts
    [atmp, btmp] = max(sg_Y);
    startPoints0 = [atmp, sg_X(btmp), gaussC];
    [f0, err0] = fit(sg_X', sg_Y', gaussEqn, 'Start', startPoints0);
    startPoints1 = [atmp, sg_X(btmp), 0.1];
    [f1, err1] = fit(sg_X', sg_Y', gaussEqn, 'Start', startPoints1);
    startPoints2 = [atmp, sg_X(btmp), 100];
    [f2, err2] = fit(sg_X', sg_Y', gaussEqn, 'Start', startPoints2);
    startPoints3 = [atmp, sg_X(btmp), 1000];
    [f3, err3] = fit(sg_X', sg_Y', gaussEqn, 'Start', startPoints3);
    startPoints4 = [atmp, sg_X(btmp), 0.01];
    [f4, err4] = fit(sg_X', sg_Y', gaussEqn, 'Start', startPoints4);
    % Choose the fit with the smallest error
    [err, idx] = min([err0.rmse, err1.rmse, err2.rmse, err3.rmse, err4.rmse]);
    switch idx
        case 1
            f1 = f0;
        case 2
            f1 = f1;
        case 3
            f1 = f2;
        case 4
            f1 = f3;
        case 5
            f1 = f4;
    end    
    f1.c = abs(f1.c);
    if plotfit == 1
        figure;
        plot(f1, sg_X, sg_Y); 
        timestamp = datestr(now, 'yyyymmdd_HHMMSS'); 
        imageFileName = ['fit_plot_' timestamp '.png'];
        fullImagePath = fullfile(savepath, imageFileName);
        saveas(gcf, fullImagePath);
        close(gcf);
    end
else
    f1 = [];
    err = [];
    disp('Fitting failed due to insufficient data points.');
end
