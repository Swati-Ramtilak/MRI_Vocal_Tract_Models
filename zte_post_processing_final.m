%% Section1 :Load Images
loading = 1;
if loading == 1
    clear; 
    close all;  
    %User selects directory of files
    imgdir = uigetdir('Please select the file folder') ;
    if imgdir(end) ~= '/'
        imgdir = [imgdir '/'];
    end
    fnext = '*.dcm';
    dicomm = 1;
    savepath = uigetdir('Please select the file folder') ;
    %Load all DICOM data into imagefnames
    imgfnames = dir(fullfile(imgdir, fnext));
    slices = length(imgfnames);
    %reads the metadata from DICOM      
    fn = [imgdir imgfnames(1).name];                                       
    dinfo_tmp = dicominfo(fn);
    %accessing metadata from DICOM
    rbw = dinfo_tmp.PixelBandwidth * double(dinfo_tmp.Width);
    rows = double(dinfo_tmp.Height); 
    cols = double(dinfo_tmp.Width);
    imsize = rows;
    nslice = slices;
    mag = zeros(nslice,rows,cols);
    % assume 16:9 screen display ratio
    imx = sqrt(nslice/(16*9)); 
    imrows = floor(9*imx+0.5);
    imcols = floor(16*imx+0.5);
    imst = 1;              
    implane = 1; 
    if imrows * imcols > nslice
        imcols = imcols - 1;
    else
        imst = max(1, floor(nslice - imrows*imcols)/2);
    end
    rfmt = [imrows imcols imst implane];    
    dinfo = [];
    dcminfo = [];                                                          
    prof = zeros(1,nslice);                                                
    for slidx = 1:nslice   
        fn = [imgdir imgfnames(slidx).name];
        dinfo{slidx} = dicominfo(fn);
        dcminfo = [dcminfo dinfo{slidx}.SliceLocation];                    
        mag(slidx,:,:) = dicomread(fn);
    end    
    [slord , slidx] = sort(dcminfo, 'descend');
    mag = mag(slidx,:,:);
    dcminfo = dcminfo(slidx);
    dinfo = dinfo(slidx);
    mag = double(mag);
end
% display original images
figure;
imagesc(reframeslice(mag, rfmt(1),rfmt(2), rfmt(3))); colormap('gray');axis off;
title('original images');
imageFileName = 'original_dicom.png';  
fullImagePath = fullfile(savepath, imageFileName);
saveas(gcf, fullImagePath);
close gcf;

%% Section 2: Define Imaging parameters
bone_water_fraction = 0.2;                                                 %fraction of water content in bone tissue, set to 20%.
T1csf = 3000;
T1tissue = 1000;
T1bone = 230;
if double(dinfo{1}.AcquisitionMatrix(4)) == 0
    TR = dinfo{1}.RepetitionTime/double(dinfo{1}.AcquisitionMatrix(3));    %Calculate TR as k space is partially filled with every repetition 
else
    TR = dinfo{1}.RepetitionTime/double(dinfo{1}.AcquisitionMatrix(4));
end
theta = dinfo{1}.FlipAngle*pi/180;         
E1csf = exp(-TR/T1csf);
E1tissue = exp(-TR/T1tissue);
E1bone = exp(-TR/T1bone);
Mcsf = (1-E1csf)*sin(theta)/(1-E1csf*cos(theta));
Mtissue = (1-E1tissue)*sin(theta)/(1-E1tissue*cos(theta));
Mbone = (1-E1bone)*sin(theta)/(1-E1bone*cos(theta));
sig_ratio = Mcsf/Mtissue;                                                  % csf vs tissue ratio
bone_sig_ratio = Mbone/Mtissue*bone_water_fraction;

%% Section 3: Setting a Threshold First Guess
gaussEqn = 'a*exp(-((x-b)/c)^2)';
Nbins = 500; Nminpts = 8; noisePeak = 1; plotfit = 1; gaussC = sqrt(5);
mag_kw = mag;
mag_kw_org = mag_kw;
N_bkwidth_left = 10;
N_bkwidth_right = 10;
[avg_org,SI_thres,binRange,N,binC] = find_avg_zte_signal(mag_kw, Nbins);
%bar(binC, N, 'hist');

%% Section 4: Fitting Histogram Model
f1 = fit_zte_hist_model(mag_kw, SI_thres, Nminpts, Nbins, binRange, N_bkwidth_left, N_bkwidth_right, gaussEqn, savepath);
if ~isempty(f1)
    avg_org = f1.b;                                                        %Mean of curve
else
    disp('error: f1 empty');
end
[f1_noise] = fit_zte_hist_model(mag_kw, SI_thres, Nminpts, Nbins, binRange, 20, 5, gaussEqn, 1, [], 1, savepath);
if ~isempty(f1_noise) && ~isempty(f1)
    SI_thres = max([f1.b - abs(f1.c) * 3, f1_noise.b+3*abs(f1_noise.c)]);
else
    disp('error: f1_noise empty');
end

%% Section 5: Local Thresholding and Smoothing
SI_thres1 = (1-sig_ratio)*avg_org;
SI_thres0 = sig_ratio*avg_org;
winZ = ceil(slices./2.^[2 3 4]);
winY = ceil(rows./2.^[2 3 4]);
winX = ceil(cols./2.^[2 3 4]);
stepZ = 16;                                                                %4 16;
stepY = 16;                                                                %4 16;
stepX = 16;                                                                %4 16;
[N,edges] = histcounts(mag_kw, Nbins);
B = zeros(size(mag_kw));
B(:,:,:) = avg_org;
Bpre = B;
for mres = 1:length(winX)
    W = [winZ(mres) winY(mres) winX(mres)];                                % window size
    sz = W(1)*W(2)*W(3);
    for i = 1:stepZ:nslice-W(1)
        for j = 1:stepY:imsize-W(2)
            for k = 1:stepX:imsize-W(3)
                zvct = i : i + W(1);
                yvct = j : j + W(2);
                xvct = k : k + W(3);
                %blk = mag_kw(zvct,yvct,xvct);
                [N, X] = histcounts(mag_kw(zvct,yvct,xvct), edges);
                N = smooth(N);
                maxN = max(N(find(X(1:end-1)>SI_thres)));
                if maxN > 0.02*sz                                          %To check if significat signal is present in the histogram
                    idx = find(N == maxN, 1, 'last');
                    avgX = X(idx);              
                    zpreidx = max([zvct(1) 2]);
                    ypreidx = max([yvct(1) 2]);
                    xpreidx = max([xvct(1) 2]);
                    avg_pre = sum(sum(sum(Bpre(zvct, yvct, xvct))))/sz;
                    nb_vox = [B(zpreidx, ypreidx, xpreidx-1),  B(zpreidx, ypreidx-1, xpreidx),   B(zpreidx-1, ypreidx, xpreidx), ...
                              B(zpreidx,ypreidx-1, xpreidx-1), B(zpreidx-1, ypreidx, xpreidx-1), B(zpreidx-1, ypreidx-1, xpreidx), ...
                              B(zpreidx-1, ypreidx-1, xpreidx-1)];
                    Bnb_avgX = mean(nb_vox);
                    maxdiff = max(abs(nb_vox - avgX));
                    if ((maxdiff <= SI_thres1) || (avgX > Bnb_avgX))
                        B(zvct,yvct,xvct) = avgX;       
                    else
                        B(zvct,yvct,xvct) = max([avgX Bnb_avgX-SI_thres1  avg_pre - SI_thres1]);
                    end
                end
            end
        end
    end
    Bpre = B;
end

%% Section 6: Local thresholding and smoothing
Br = zeros(size(mag_kw));
Br(:,:,:) = avg_org;
Bpre = Br;
for mres = 1:length(winX)
    W = [winZ(mres) winY(mres) winX(mres)];                                % window size
    sz = W(1)*W(2)*W(3);
    for i = nslice:-stepZ:W(1)+1
        for j = imsize:-stepY:W(2)+1
            for k =imsize:-stepX:W(3)+1
                zvct = i :-1: i - W(1);
                yvct = j :-1: j - W(2);
                xvct = k :-1: k - W(3);
                [N, X] = histcounts(mag_kw(zvct, yvct, xvct), edges);
                N = smooth(N);
                maxN = max(N(find(X(1:end-1)>SI_thres)));
                if maxN > 0.02*sz
                    idx = find(N == maxN, 1, 'last');
                    avgX = X(idx);              
                    zpreidx = min([zvct(1) nslice-1]);
                    ypreidx = min([yvct(1) imsize-1]);
                    xpreidx = min([xvct(1) imsize-1]);
                    avg_pre = sum(sum(sum(Bpre(zvct, yvct, xvct))))/sz;
                    nb_vox = [Br(zpreidx, ypreidx, xpreidx+1),  Br(zpreidx, ypreidx+1, xpreidx),   Br(zpreidx+1, ypreidx, xpreidx), ...
                              Br(zpreidx,ypreidx+1, xpreidx+1), Br(zpreidx+1, ypreidx, xpreidx+1), Br(zpreidx+1, ypreidx+1, xpreidx), ...
                              Br(zpreidx+1, ypreidx+1, xpreidx+1)];
                    Bnb_avgX = mean(nb_vox);
                    maxdiff = max(abs(nb_vox - avgX));
                    if ((maxdiff <= SI_thres1) || (avgX > Bnb_avgX))
                        Br(zvct,yvct,xvct) = avgX;       
                    else
                        Br(zvct,yvct,xvct) = max([avgX Bnb_avgX-SI_thres1  avg_pre - SI_thres1]);
                    end
                end
            end
        end
    end
    Bpre = Br;
end

%% Section 7: Setting a Threshold After smoothing and local thresholding
Bdiff = abs(B - avg_org);
Brdiff = abs(Br - avg_org);
Bf = B.*(Bdiff>=Brdiff) + Br.*(Bdiff<Brdiff);
idxmtx = (mag_kw>SI_thres0).*(mag_kw>Bf);
idx = find(idxmtx > 0);
Bf(idx) = mag_kw(idx);
Bgauss = imgaussfilt3(Bf, 15);
mag_kw_B = mag_kw./Bgauss; 
plotfit = 1;
gauss_thres = 1e-1; 
N_bkwidth_left = 10;
N_bkwidth_right = 10;
noisePeak = 0;
mag_kw_B_filt = imgaussfilt3(mag_kw_B,1); 
[avg_org,SI_thres,binRange,N,binC ] = find_avg_zte_signal(mag_kw_B, Nbins);
%Xbar(binC, N, 'hist');
f1 = fit_zte_hist_model(mag_kw_B, SI_thres, Nminpts, Nbins, binRange, N_bkwidth_left, N_bkwidth_right, gaussEqn, noisePeak, [], plotfit, savepath);
SI_thres_low = f1.b - 3*abs(f1.c); 
SI_thres_high = f1.b + 3*abs(f1.c);

%% Section 8: Generate Masks
mask = mag_kw_B_filt>SI_thres_low*1.1;
mask1 = mask;
for m = 1:nslice 
     slice = int16(squeeze(mask(m,:,:)));
     slice = imfill(slice);
     mask1(m,:,:) = slice;
end
mask = mask1;
for m = 1:cols 
     slice = int16(squeeze(mask(:,m,:)));
     slice = imfill(slice);
     mask1(:,m,:) = slice;
end
mask = mask1;
for m = 1:rows  
     slice = int16(squeeze(mask(:,:,m)));
     slice = imfill(slice);
     mask1(:,:,m) = slice;
end
mask1 = imgaussfilt(double(mask1),3) > 0;
se = strel('sphere', 8);                                                   %MOD - cube, sphere
mask1er = imerode(mask1, se);
mask1gauss = imgaussfilt3(double(mask1er),3); 
mask_gauss = mag_kw_B_filt > 0.3;                                          %MOD
mask_gauss = bwareaopen(1-mask_gauss, 2000);                               %MOD 3000
mask_gauss = imgaussfilt3(double(1-mask_gauss), 2);

%% Section 9: Creating a Mask for Air Regions
remove_air = 0 ;                                                            % Assign Value 1 to use this section of code
if remove_air == 1
    nseeds = 80;
    rfmt1 = [5 10 75 5];
    region_tmp = zeros(size(mag_kw_B));
    region_final = zeros(size(mag_kw_B));
    %thres = 0.03;                                                          % MOD - %0.1 0.03;
    thres = 0.1; 
    dist_sz = {50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
        50 50 50 50 50 50 50 50 50 50 ...
              };              
    decrstr = ['Pick ' num2str(nseeds) ' locations for air'];
    fig1 = figure('Name', decrstr, 'NumberTitle', 'off');
    imagesc(reframeslice(mag_kw_B,rfmt1(1),rfmt1(2), rfmt1(3), 0, rfmt1(4)));axis off; colormap('gray');
    set(fig1, 'Units', 'pixels');
    [posx, posy] = ginput(nseeds);
    if length(posx) < nseeds
        disp('not enough points!');
        return;
    end
    posx = round(posx);
    posy = round(posy);
    loc = {};
    for m = 1:nseeds
         if m == 14
             m= m+1;
         end
        loc{m} = location3d(posx(m), posy(m), rfmt1, imsize, imsize);
        [P, region_tmp] = myregionGrowing(mag_kw_B, loc{m}, thres*avg_org, dist_sz{m}, false, true, false,[],0.7*avg_org);
        region_final = region_final + double(region_tmp);
    end
    mask_air = region_final > 0;
    maska_gauss = double(imfill(mask_air,'holes'));
    maska_gauss = imgaussfilt3(1 - maska_gauss, 2);
end

%% Section 10: Combine masks and remove air
[avg_org, SI_thres, binRange] = find_avg_zte_signal(mag_kw_B);
mag_kw_B = mag_kw_B/avg_org;
maginv = abs(mag_kw_B-max(max(max(mag_kw_B))));
logMag = -log(mag_kw_B);
logMag = logMag + abs(min(min(min(logMag))));
npow = 1;
final_mask = (mask1gauss.*mask_gauss).^npow;
logMag_mask = logMag.*final_mask;
if remove_air == 1
    logMag_mask = logMag_mask.*maska_gauss;
end
x = find(isnan(logMag_mask));
logMag_mask(x) = 0;
x = find(isinf(logMag_mask));
logMag_mask(x) = 0;
mask0 = imgaussfilt3(mag_kw, 25);                                          %MOD - 19, 25
mask0 = mask0>0.32*max(mask0(:)); 
mask0 = imgaussfilt3(double(mask0), 1);
mag_med = imgaussfilt3(mag, 3);
mask2 = mag_med > 0.1 * max(mag_med(:));                                   %MOD - 0.06 0.18; 0.1
mask2 = bwareaopen(1 - mask2, 7000);                                       %MOD - 2000, 3000
mask2 = imgaussfilt3(double(1 - mask2), 2);
logMag_mask = logMag_mask .* (mask2.^4) .* mask0;

%% Section 11:Saving Data
imagesc(reframeslice(mag_kw_B, rfmt(1), rfmt(2), rfmt(3)));axis off;colormap('gray');
imageFileName = 'processed_dicom.png';  
fullImagePath = fullfile(savepath, imageFileName);
saveas(gcf, fullImagePath);
close gcf;
savedcm = 1;
warning('off', 'all');
if savedcm == 1
    loctmp = find(dinfo{1}.Filename == '/', 2, 'last');
    if ~isempty(loctmp)
        loc = loctmp(1);
    else
        loctmp = find(dinfo{1}.Filename == '\', 2, 'last');
        loc = loctmp(1);
    end    
    logMag_mask = logMag_mask * 32768 / max(logMag_mask(:));               %MOD
    mag_kw_B1 = mag_kw_B * 32768 / max(mag_kw_B(:)) .* mask2 .* mask1gauss;%MOD
    uid1 = dicomuid;
    uid2 = dicomuid;
    uid3 = dicomuid;
    mag_cor_dir = [dinfo{1}.Filename(1:loc) 'mag_cor'];
    proc_dir = [dinfo{1}.Filename(1:loc) 'processed'];
    if ~exist(mag_cor_dir, 'dir')
    mkdir(mag_cor_dir);
    end
    if ~exist(proc_dir, 'dir')
    mkdir(proc_dir);
    end
    for m = 1:nslice
    %slice = adapthisteq(int16(squeeze(mag_kw_B1(m, :, :))),'NumTiles', [8 8], 'ClipLimit', 0.01); 
    slice = int16(squeeze(mag_kw_B1(m,:,:)));
    fn = sprintf('%s/mag.%.3d.dcm', mag_cor_dir, m);
    dinfo{m}.SeriesNumber = dinfo{m}.SeriesNumber * 10 + 4000;
    dinfo{m}.SeriesInstanceUID = uid3;
    dicomwrite(slice, fn, dinfo{m});    
    %slice = adapthisteq(int16(squeeze(logMag_mask(m, :, :))),'NumTiles', [8 8], 'ClipLimit', 0.01);
    slice = int16(squeeze(logMag_mask(m,:,:)));
    fn = sprintf('%s/mag.%.3d.dcm', proc_dir, m);
    dinfo{m}.SeriesNumber = dinfo{m}.SeriesNumber + 1;
    dinfo{m}.SeriesInstanceUID = uid1;
    dicomwrite(slice, fn, dinfo{m});
    end
end