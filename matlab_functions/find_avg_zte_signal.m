function [avg_org SI_thres binRange N binC] = find_avg_zte_signal(img, Nbins, Npkmin)
%Condition to check number of arguments
if nargin < 3                     
    Npkmin = 100;
end
if nargin < 2
    [N,edges] = histcounts(img);
else
    [N,edges] = histcounts(img, Nbins);
end
neg_thres = 2e-4;                                                          %DEL
N(1:3) = 1; 
Nidx = find(N>0);
binC = edges(1:end-1)+diff(edges)/2;                                       %To plot hist
N1 = smooth(N(Nidx));
N1= movmedian(N1, 3);
N1 = smooth(smooth(smooth(smooth(N1))));
if Npkmin < max(N(3:end))*0.01
    Npkmin = ceil(max(N)*0.01);
end
edges1 = edges(Nidx);
N1diff = diff(N1);
[pks , locs] = findpeaks(N1, 'MinPeakHeight', Npkmin, 'SortStr','descend');% ideally 2 pks

if isempty(pks)
    maxN1 = max(N1);
    maxN1idx = find(N1 == maxN1, 1,'first');
    avg_org = edges(maxN1idx);
    binRange = edges(end);
    SI_thres = 0;  
else
    locs1 = sort(locs);
    if length(pks) > 2 
        if locs(1) < locs(2)
            Negidx = locs(1);                                              %locs1(2);
        else
            Negidx = locs(2);
        end
    else
        Negidx = locs1(1);
    end
    idxvct = Negidx:length(N1diff);
    edges1tmp = edges1(idxvct);
    N1difftmp = N1diff(idxvct);
    zcd_N1 = sign(N1difftmp(1:end-1)) .* sign(N1difftmp(2:end));
    N1diffidx = 1;                                                         
    idx_st0 = find(zcd_N1>0, 1, 'first');
    idx_st1 = find(N1difftmp>0, 1,'first');
    idx_st = max([idx_st0, idx_st1]);
    N1diffidx = idx_st-1+find(zcd_N1(idx_st:end) < 0, 1, 'first');
    SI_thres = edges1tmp(N1diffidx);
    binRangeidx = find(N>Npkmin,1, 'last');
    binRange = [0, edges(binRangeidx)];
    avg_org = edges1(find(N1 == N1(locs1(end)),1, 'last'));
end
