function entropy = calc_entropy(points,boxsize,binsize)
%% normalize Point
Fmean = mean(points,2);
tmp = repmat(Fmean,1,size(points,2));
%RF = (points+1e-9)./(tmp+1e-9) -1;
RF = points./tmp -1;
RF_max = max(RF,[],2);
RF_min = min(RF,[],2);
%% set entropy box
tmpcut = 2*boxsize/binsize;
boxcuts = -boxsize + [0:binsize]*tmpcut;
%% calcuate pvalue
pvalues = zeros(size(points,1),binsize);
for i = 1:binsize
    pvalues(:,i) = sum(RF > boxcuts(i) & RF <= boxcuts(i+1),2);
end
validSize = sum(RF > -boxsize & RF < boxsize,2); 
%pvalues = pvalues/size(points,2);
pvalues = pvalues./repmat(validSize,1,size(pvalues,2));
pvalues = pvalues + 1e-14; % set minimal pvalue
%% calculate entropy
entropy = -pvalues .*log(pvalues);
entropy = sum(entropy,2);
%entropy(~vid) = 0;
    
    