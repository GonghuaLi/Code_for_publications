function x = cal_CV(v)
x = std(v')'./abs(mean(v,2));