function [c,lags] = corrcoef_windows(wabp,wbfi)
%%
s_range = (1:ceil(length(wabp)/2))-1;
wabp_n = normalize(wabp,'range',[0 1]);
wbfi_n = normalize(wbfi,'range',[0 1]);
for s=4:length(s_range)
    wabp_s = wabp_n(s_range(s)+1:end);
    wbfi_s = wbfi_n(1:length(wbfi)-s_range(s));
    corr = corrcoef(wabp_s,wbfi_s);
    c(s) = corr(1,2);
end
lags = s_range;
end