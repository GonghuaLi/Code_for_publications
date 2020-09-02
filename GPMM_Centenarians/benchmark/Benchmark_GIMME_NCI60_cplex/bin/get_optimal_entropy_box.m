function ia = get_optimal_entropy_box(points)
Fmean = mean(points,2);
%%
tmp = repmat(Fmean,1,size(points,2));
RF = points./tmp -1;
RF_max = max(RF,[],2);
RF_min = min(RF,[],2);
%%
boxsize = [1:8000]/1000;
percRxns_inbox = zeros(length(boxsize),1);
for i = 1:length(boxsize)
    percRxns_inbox(i) = sum(RF_max < boxsize(i) & RF_min > -boxsize(i))/length(RF_max);
end

[xx,ia] = min(abs(percRxns_inbox - 0.95))
ia = ia/1000;