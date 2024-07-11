function [t,f,fx,A_fx_out,p_fx_out,px,A_px_out,p_px_out] = freq_domain_analysis(gwin,fs)
% For implementation example see /tmu_code/under_development/2022_4_20/transfer_function_2022420.m

A_px_out = [];
A_fx_out = [];
p_px_out = [];
p_fx_out = [];

for w=1:height(gwin)
    bwin = gwin(w,:);
    wt = bwin.wt{:};
    wbfi = bwin.wbfi{:};
    wn_bfi = wbfi./mean(wbfi);
    wabp = bwin.wabp{:};
    wn_abp = wabp./mean(wabp);

    wnum = bwin.wnum(:);

    T = 1/fs;
    L = length(wbfi);
    t = (0:L-1)*T;

    % FFT FLOW
    fx = fft(wbfi./mean(wbfi));
    % One-sided amplitude
    afx2 = abs(fx/L);
    afx1 = afx2(1:floor(L/2)+1);
    afx1(2:end-1) = 2*afx1(2:end-1);
    f = fs*(0:floor(L/2))/floor(L);
    A_fx_out(w,:) = afx1;
    % One-sided phase
    pfx2 = angle(fx/L);
    pfx1 = pfx2(1:floor(L/2)+1);
    pfx1(2:end-1) = 2*pfx1(2:end-1);
    p_fx_out(w,:) = pfx1;    
    
    % FFT PRESSURE
    px = fft(wabp./mean(wabp));
    % One-sided amplitude
    apx2 = abs(px/L);
    apx1 = apx2(1:floor(L/2)+1);
    apx1(2:end-1) = 2*apx1(2:end-1);
    A_px_out(w,:) = apx1;
    % One-sided phase
    ppx2 = angle(px/L);
    ppx1 = ppx2(1:floor(L/2)+1);
    ppx1(2:end-1) = 2*ppx1(2:end-1);
    p_px_out(w,:) = ppx1;    
    
end

end