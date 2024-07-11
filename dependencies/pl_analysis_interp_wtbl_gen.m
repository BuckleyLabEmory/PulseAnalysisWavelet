function wtbl_n = pl_analysis_interp_wtbl_gen(ftbl_n,wtbl_n,gamma,numnans)
% Please see /Users/taraurner/BuckleyLab/CrCP/CrCP_tmu/Notebooks/N2021_06_10_PulsatilityAnalysis.mlx
% and N2021_06_11_wtbl_analysis.mlx
% and associated folder for interactive testing of this workflow

    w_a = table();
    windows = height(wtbl_n);
    
    % interped data for xcorr
    t_interp = ftbl_n.t_bfi_interp{:};
    bfi_interp = ftbl_n.bfi_interp{:};
    abp_interp = ftbl_n.abp_binterp{:};
    
    % raw bfi and downsampled ABP for crcp
    t = ftbl_n.t_bfi{:};
    bfi = ftbl_n.bfi{:};
    abp = ftbl_n.abp_downsample{:};

    for w = 1:windows
        % Check for NaNs

        wstart = wtbl_n{w,'wstart'};
        wend = wtbl_n{w,'wend'};
        
        wt = t(find(t==wstart):find(t==wend));
        wbfi = bfi(find(t==wstart):find(t==wend));
        wabp = abp(find(t==wstart):find(t==wend));
        
        wt_interp = t_interp(wstart<=t_interp & t_interp<=wend);
        wbfi_interp = bfi_interp(wstart<=t_interp & t_interp<=wend);
        wabp_interp = abp_interp(wstart<=t_interp & t_interp<=wend);

        [~,~,wbfi,wabp,wkeep(w),discard_reason]...
        = window_test_generalized3(wbfi,wabp,wt,numnans);
        
%%
    
        if ~wkeep(w)
            fprintf('bad window')
            cpp = NaN; crcp = NaN; tau = NaN; cppmich = NaN; crcpmich = NaN;
            Fratio = NaN; Pratio = NaN; C = NaN; hr = NaN; AP = NaN; AF = NaN;
            f = NaN; I = NaN; Ihr = NaN; whr = NaN; Ppf = NaN; H = NaN; phi = NaN;
            peaks = getpeaks3(wabp_interp,ftbl_n.fs_abp);
            if contains(discard_reason,'NaNs Threshold - ABP')
                shift_maxcorr = NaN; tshift_maxcorr = NaN; maxcorr = NaN;
                %continue
            else
              fprintf('else') 
              non_nans = ~isnan(wabp) & ~isnan(wbfi);
              
              
              if strfind(ftbl_n.corr_method{:},'xcorr')
                % Xcorr lag calculation
                [c,lags] = xcorr(wabp(non_nans),wbfi(non_nans),'normalized');
                shift_maxcorr = lags(c==max(c));
                tshift_maxcorr = shift_maxcorr * mean(diff(t));
                maxcorr = max(c);                
              elseif strfind(ftbl_n.corr_method{:},'corrcoef')
                % Corrcoef (added in October 2021)
                [c,lags] = corrcoef_windows(wabp(non_nans),wbfi(non_nans));
                shift_maxcorr = lags(c==max(c));
                tshift_maxcorr = shift_maxcorr * mean(diff(t));
                maxcorr = max(c);  
              end
            
            end
            %continue
%             shift_maxcorr = NaN; tshift_maxcorr = NaN; maxcorr = NaN;

%%

        else

            if strfind(ftbl_n.corr_method{:},'xcorr')
                % Xcorr lag calculation
                [c,lags] = xcorr(wabp,wbfi,'normalized');
                shift_maxcorr = lags(c==max(c));
                tshift_maxcorr = shift_maxcorr * mean(diff(t));
                maxcorr = max(c);                
            elseif strfind(ftbl_n.corr_method{:},'corrcoef')
            % Corrcoef (added in October 2021)
                [c,lags] = corrcoef_windows(wabp,wbfi);
                shift_maxcorr = lags(c==max(c));
                tshift_maxcorr = shift_maxcorr * mean(diff(t));
                maxcorr = max(c);  
            end
          

            % CrCP calculation
            [cpp,crcp,tau,cppmich,crcpmich,Fratio,Pratio,C,hr,AP,AF,f,I,Ihr,whr,Ppf,H,phi]...
                = compute_cpp_windkessel_nofigs(wabp,wbfi,1/mean(diff(t)),gamma);
            
            % Waveforms
            
            peaks = getpeaks3(wabp_interp,ftbl_n.fs_abp);
            %peaks_20hz = getpeaks3(wabp,ftbl_n.fs_bfi);
            
        end

        % Add data to table
        w_a_dat = [{w},{wt},{wbfi},{wabp},{wkeep(w)},...
            {shift_maxcorr},{tshift_maxcorr}, {maxcorr},...
            {cpp},{crcp},{tau},{cppmich},{crcpmich},{Fratio},{Pratio},...
            {C},{hr},{AP},{AF},{f},{I},{Ihr},{whr},{Ppf},{H},{phi},...
            {wt_interp},{wabp_interp},{peaks},{ftbl_n.corr_method{:}}];
        w_a_names = {'wnum','wt','wbfi','wabp','wkeep',...
            'shift_maxcorr','tshift_maxcorr', 'maxcorr',...
            'cpp','crcp','tau','cppmich','crcpmich','Fratio','Pratio',...
            'C','hr','AP','AF','f','I','Ihr','whr','Ppf','H','phi',...
            'wt_125hz','wabp_125hz','peaks','corr_method'};
        w_a_dat = cell2table(w_a_dat);
        w_a_dat.Properties.VariableNames = w_a_names;
        
        if isfield(ftbl_n,'delineator_20hz')
            w_a_dat.delineator_20hz = 1;
        end
        
        w_a = [w_a ; w_a_dat];
    end
    wtemp = wtbl_n;
    wtbl_n = [wtbl_n,w_a];

end