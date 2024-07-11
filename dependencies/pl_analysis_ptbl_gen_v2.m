function pl = pl_analysis_ptbl_gen_v2(ftbl,wtbl_n)
%%
wtbl_nk = wtbl_n;
 
pl = table();

for w=1:height(wtbl_nk)
    
    if wtbl_nk.wnum(w)==234
        fprintf('stop for debugging')
    end
    
    if isnan(wtbl_nk.maxcorr(w)) || wtbl_nk.maxcorr(w) < 0.3 % this is ABP calibration
        continue
    end
    
    t_abp = wtbl_nk.wt_125hz{w};
    t_bfi = wtbl_nk.wt{w};
    abp = wtbl_nk.wabp_125hz{w};
    abp_ds = wtbl_nk.wabp{w};
    bfi = wtbl_nk.wbfi{w};
    
    abp_p_idx = wtbl_nk.peaks(w).a_p;
    if isempty(abp_p_idx)
        continue
    elseif size(abp_p_idx,2) < 2
        continue
    end
    abp_pulse_bounds = abp_p_idx(1,:)';
    
    abp_p_s = abp_pulse_bounds;
    abp_p_e = circshift(abp_pulse_bounds,-1);
    
    abp_p_start = abp_p_s(1:end-1);
    abp_p_end = abp_p_e(1:end-1);
    abp_p_length = num2cell(abp_p_end - abp_p_start);
    abp_p_start = num2cell(abp_p_start);
    abp_p_end = num2cell(abp_p_end);
    n_pulses = size(abp_p_idx,2);
    
    
    % ABP pulse morphology
    a_onset = num2cell(wtbl_nk.peaks(w).ao(1:n_pulses-1))';
    a_peak = num2cell(wtbl_nk.peaks(w).ap(1:n_pulses-1))';
    a_dicron = num2cell(wtbl_nk.peaks(w).ad(1:n_pulses-1))';
    tpulse = struct2table(wtbl_nk.peaks(w).t);
    ypulse = struct2table(wtbl_nk.peaks(w).y,'AsArray',true);
    
    t_aop = num2cell(tpulse.aop(1:n_pulses-1))';
    t_aod = num2cell(tpulse.aod(1:n_pulses-1))';
    t_apd = num2cell(tpulse.apd(1:n_pulses-1))';
    
    y_aop = num2cell(ypulse.aop{1}(1:n_pulses-1));
    y_aod = num2cell(ypulse.aod{1}(1:n_pulses-1));
    y_apd = num2cell(ypulse.adp{1}(1:n_pulses-1));

    % Time series
    bfi_p_length = {};
    bfi_p_start = {};
    bfi_p_end = {};
    abp_pt = {};
    bfi_pt = {};
    abp_ln_pt = {};
    abp_p = {};
    bfi_p = {};
    abp_ln_p = {};
    abp_ds_p = {};
    bfi_delta_t = {};
    
    % BFI fit metrics
    bfi_avgcorr = {};
    bfi_avgcorr_sm = {};
    bfi_last_taus = {};
    bfi_curvefit = {};
    bfi_meanerror = {};
    bfi_bound = [];
    bfi_bound_loc = [];
    abp_pulses = ones(n_pulses,1);
 %%   
    for p = 1:n_pulses
        %%%% We will consider the overlap region to be the last BFI point before ABP 
        %onset to the last BFI point before the next ABP onset
        below_bound = find(t_bfi <= t_abp(abp_pulse_bounds(p))-wtbl_nk.tshift_maxcorr(w));
        if isempty(below_bound)
            abp_pulses(p) = 0;
            continue
        end
        bfi_bound(p) = t_bfi(below_bound(end));
        bfi_bound_loc(p) = find(t_bfi == t_bfi(below_bound(end)));

    end 
    
    bfi_s = bfi_bound_loc;
    bfi_e = circshift(bfi_s,-1);
    bfi_ends = bfi_e(1:end-1)';
    bfi_starts = bfi_s(1:end-1)';
    bfi_idx = [bfi_starts,bfi_ends];
    
    use_abp = find(abp_pulses);
    use_abp = use_abp(1:end-1);
 %%       
   for l = 1:length(use_abp);
        p = use_abp(l);
        bfi_pl = zeros(size(t_bfi));
        bfi_pl(bfi_idx(p,1):bfi_idx(p,end))=1;
        bfi_pt_range = logical(bfi_pl);
        abp_pt_range = (abp_p_start{p}:abp_p_start{p} + abp_p_length{p});

        abp_ptn = {t_abp(abp_pt_range)};
        abp_pt = [abp_pt; abp_ptn];
        bfi_ptn = {t_bfi(bfi_pt_range)};
        bfi_pt = [bfi_pt; bfi_ptn];
        abp_pn = {abp(abp_pt_range)};
        abp_p = [abp_p ; abp_pn];
        bfi_pn = {bfi(bfi_pt_range)};
        bfi_p = [bfi_p ; bfi_pn];
        if isempty(bfi_pn{1})
            bfi_delta_tn = {[0]};
            bfi_p_start(l,1) = {[0]};
            bfi_p_end(l,1) = {[0]}; 
            bfi_p_length(l,1) = {[0]};
        else
            bfi_delta_tn = num2cell(bfi_ptn{1}(1) - abp_ptn{1}(1));
            bfi_range = find(bfi_pt_range);
            bfi_p_start(l,1) = {bfi_range(1)};
            bfi_p_end(l,1) = {bfi_range(end)};
            bfi_p_length(l,1) = {bfi_range(end) - bfi_range(1)};
        end
        bfi_delta_t = [bfi_delta_t ; bfi_delta_tn];
        abp_ds_pn = {abp_ds(bfi_pt_range)};
        abp_ds_p = [abp_ds_p ; abp_ds_pn];
        
        % BFI fit metrics
        TimeAxisDCS_idx = find(ismember(ftbl.t_bfi{1},bfi_ptn{1}));
        bfi_avgcorr(l,1) = {ftbl.g2fit{1}.avg_corr_frame(TimeAxisDCS_idx)};
        bfi_avgcorr_sm(l,1) = {ftbl.g2fit{1}.corr_smooth(TimeAxisDCS_idx)};
        bfi_last_taus(l,1) = {ftbl.g2fit{1}.last_tau(TimeAxisDCS_idx)};
        bfi_curvefit(l,1) = {ftbl.g2fit{1}.curvefit2avg(TimeAxisDCS_idx)};
        bfi_meanerror(l,1) = {ftbl.g2fit{1}.meanerror(TimeAxisDCS_idx)};
      
    end
    
    pl_info_n = [abp_p_start(use_abp),...
        abp_p_end(use_abp),...
        abp_p_length(use_abp),...
        bfi_p_start,...
        bfi_p_end,...
        bfi_p_length,...
        bfi_delta_t,...
        a_onset(use_abp),...
        a_peak(use_abp),...
        a_dicron(use_abp),...
        t_aod(use_abp),...
        t_aop(use_abp),...
        t_apd(use_abp),...
        y_aod(use_abp),...
        y_aop(use_abp),...
        y_apd(use_abp),...
        abp_pt,...
        bfi_pt,...
        abp_p,...
        bfi_p,...
        abp_ds_p,...
        bfi_avgcorr,...
        bfi_avgcorr_sm,...
        bfi_curvefit,...
        bfi_last_taus,...
        bfi_meanerror];
 %%
    pl_info_n = cell2table(pl_info_n);
    pl_info_n.Properties.VariableNames = {'abp_p_start',...
        'abp_p_end',...
        'abp_p_length',...
        'bfi_p_start',...
        'bfi_p_end',...
        'bfi_p_length',...
        'bfi_delta_t',...
        'a_onset',...
        'a_peak',...
        'a_dicron',...
        't_aod',...
        't_aop',...
        't_apd',...
        'y_aod',...
        'y_aop',...
        'y_apd',...
        'abp_pt',...
        'bfi_pt',...
        'abp_p',...
        'bfi_p',...
        'abp_ds_p',...
        'bfi_avgcorr',...
        'bfi_avgcorr_sm',...
        'bfi_curvefit',...
        'bfi_last_taus',...
        'bfi_meanerror'};
 %%
    pl_winfo_n = repmat(wtbl_nk(w,1:8),length(use_abp),1);

    pl_n = [pl_winfo_n,pl_info_n];
    pl = [pl;pl_n];
end

end