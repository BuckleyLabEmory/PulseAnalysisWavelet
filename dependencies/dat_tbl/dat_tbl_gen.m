function dat_tbl = dat_tbl_gen(wtbl,pl)

%%
wtbl_n = wtbl;
%Length Normalized
wtbl = wtbl_n;
states = {'HC'};
wstates = categorical(wtbl.state);
wtbl = wtbl(wstates == states{1},:);
wstates = categorical(wtbl.state);
wsds = categorical(wtbl.state);

placehold = cell(height(wtbl),1);
clrs = lines(4);


bin_size = 0.05;

%Overlay pulses

wtbl.bfi_pstack = placehold;
wtbl.bfi_tstack = placehold;
wtbl.abp_pstack = placehold;
wtbl.abp_tstack = placehold;
wtbl.abp_dicron = placehold;
wtbl.clr = placehold;
wtbl.abp_peaks = placehold;
wtbl.abp_pk_next = placehold;
wtbl.abp_tpk_next = placehold;

w_empty = table();
empty_idx = zeros(height(wtbl),1);
for w=1:height(wtbl)
        
       wstart = wtbl.wstart(w);
       wshift = wtbl.tshift_maxcorr(w);
       wstate = wtbl.state{w};
       wsds = wtbl.sds(w);
       
       wt_125hz = wtbl.wt_125hz{w};
       
       if wsds==1 & strfind(wstate,states{1})
           wtbl.clr(w) = {clrs(1,:)};
       elseif wsds==2 & strfind(wstate,states{1})
            wtbl.clr(w) = {clrs(2,:)};
       elseif wsds==1 & strfind(wstate,states{2})
           wtbl.clr(w) = {clrs(3,:)};
       elseif wsds==2 & strfind(wstate,states{2})
           wtbl.clr(w) = {clrs(4,:)};
       end
       
       pl_w = pl(pl.wstart == wstart & pl.sds == wsds,:);
       len_abp = cellfun(@(x) length(x),pl_w.abp_pt);
       mlen_abp = max(len_abp);
       len_bfi = cellfun(@(x) length(x),pl_w.bfi_pt);
       mlen_bfi = max(len_bfi);
       bfi_tstack = nan(height(pl_w),mlen_bfi);
       bfi_pstack = nan(height(pl_w),mlen_bfi);
       abp_tstack = nan(height(pl_w),mlen_abp);
       abp_pstack = nan(height(pl_w),mlen_abp);
       abp_dicron = nan(height(pl_w),1);
       abp_peaks = cell(height(pl_w),1);
       abp_pk_next = nan(height(pl_w),1);
       abp_tpk_next = nan(height(pl_w),1);
       
       for p=1:height(pl_w)
            bfi_tstack(p,1:len_bfi(p)) = pl_w.bfi_pt{p} - pl_w.abp_pt{p}(1) + wshift; % bfi shifted by xcorr amount
            bfi_pstack(p,1:len_bfi(p)) = pl_w.bfi_p{p};
            abp_tstack(p,1:len_abp(p)) = pl_w.abp_pt{p} - pl_w.abp_pt{p}(1); % abp normalized to start at 0
            abp_pstack(p,1:len_abp(p)) = pl_w.abp_p{p};
            abp_dicron(p) = wt_125hz(pl_w.a_dicron(p)) - pl_w.abp_pt{p}(1);
            [~,pks] = findpeaks(abp_pstack(p,1:len_abp(p)),'MinPeakDistance',0.2);
            tpks = abp_tstack(p,pks);
            abp_peaks(p) = {tpks};
            next_pks = find(tpks > abp_dicron(p));
            if isempty(next_pks)
               next_pk = NaN;
               next_tpk = NaN;
            else
                next_pk = next_pks(1);
                next_tpk = tpks(next_pk);                
            end
            abp_pk_next(p) = next_pk;
            abp_tpk_next(p) = next_tpk;
       end
       if isempty(bfi_tstack)
           w_empty = [w_empty;wtbl(w,:)];
           empty_idx(w)=1;
           fprintf('empty')
       end
       wtbl.bfi_tstack(w) = {bfi_tstack};
       wtbl.bfi_pstack(w) = {bfi_pstack};
       wtbl.abp_tstack(w) = {abp_tstack};
       wtbl.abp_pstack(w) = {abp_pstack};
       wtbl.abp_dicron(w) = {abp_dicron};
       wtbl.abp_peaks(w) = {abp_peaks};
       wtbl.abp_pk_next(w) = {abp_pk_next};
       wtbl.abp_tpk_next(w) = {abp_tpk_next};
       
end

%%
wtbl = wtbl(~empty_idx,:);

%Normalize ABP on [0 1] (shape normalized)
placehold = cell(height(wtbl),1);
wtbl.abp_tnorm = placehold;
wtbl.bfi_tnorm = placehold;
wtbl.abp_dnorm = placehold;
wtbl.abp_pnorm = placehold;
wtbl.abp_avg_l = placehold;
for w=1:height(wtbl)
    abp_tstack = wtbl.abp_tstack{w};
    bfi_tstack = wtbl.bfi_tstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    bfi_pstack = wtbl.bfi_pstack{w};
    abp_dicron = wtbl.abp_dicron{w};
    abp_pnext = wtbl.abp_tpk_next{w};
    abp_ends = max(abp_tstack');
    norm_factor = (1./abp_ends)';
    abp_tnorm = abp_tstack.*norm_factor;
    bfi_tnorm = bfi_tstack.*norm_factor;
    abp_dnorm = abp_dicron.*norm_factor;
    abp_pnorm = abp_pnext.*norm_factor;
    wtbl.abp_tnorm(w) = {abp_tnorm};
    wtbl.bfi_tnorm(w) = {bfi_tnorm};
    wtbl.abp_dnorm(w) = {abp_dnorm};
    wtbl.abp_pnorm(w) = {abp_pnorm};
    wtbl.abp_avg_l(w) = {mean(abp_ends)};
%     if w==1
%         figure
%         for p=1:size(abp_tstack,1)
%             bfi_norm = normalize(bfi_pstack(p,:),'range',[0 1]);
%             abp_norm = normalize(abp_pstack(p,:),'range',[0 1]);
%             plot(abp_tnorm(p,:),abp_norm,'.-','Color',[1.00,0.50,0.50],'MarkerSize',10,'LineWidth',0.5)
%             hold on
%             plot(bfi_tnorm(p,:),bfi_norm,'.-','color',[0.70,0.90,1.00],'MarkerSize',10,'LineWidth',0.5)
%             plot(wtbl.abp_dnorm{w}(p),abp_norm(abp_tnorm(p,:) == wtbl.abp_dnorm{w}(p)),'o','MarkerSize',10,'color','black','Linewidth',2)
%             plot(wtbl.abp_pnorm{w}(p),abp_norm(abp_tnorm(p,:) == wtbl.abp_pnorm{w}(p)),'o','MarkerSize',10,'color',[0.08,0.00,0.5],'Linewidth',2)
% 
%             ylim([0 1.1])
%         end
%            %graph_beautify(gca,30,0.5,1)
% 
%     end
end
%Bin and average
%%
wtbl.bins_ln = placehold;
wtbl.tbin_ln = placehold;
wtbl.avg_bfi_ln = placehold;
wtbl.avg_abp_ln = placehold;
wtbl.std_bfi_ln = placehold;
wtbl.std_abp_ln = placehold;
wtbl.t_dicron_ln = placehold;
wtbl.loc_dicron_ln = placehold;
wtbl.t_next_ln = placehold;
wtbl.loc_next_ln = placehold;
wtbl.psf_ln = placehold;
wtbl.psp_ln = placehold;
wtbl.psf_loc_ln = placehold;
wtbl.psp_loc_ln = placehold;
wtbl.p_auc_ln = placehold;
wtbl.f_auc_ln = placehold;
wtbl.binsize_ln = placehold;


for w=1:height(wtbl)
    bins = -(bin_size/2):bin_size:(1+(bin_size/2));
    nbins = length(bins);

    bfi_tnorm = wtbl.bfi_tnorm{w};
    if isempty(bfi_tnorm)
        continue
    end
    abp_tnorm = wtbl.abp_tnorm{w};
    bfi_pstack = wtbl.bfi_pstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    bfi_clr = wtbl.clr{w}(:);
    state = wtbl.substate{w};
    sds_p = wtbl.sds(w);
    name = wtbl.name{w};
    wnum = wtbl.wnum(w);
    
    t_dicron = mean(wtbl.abp_dnorm{w});
    %t_next = mean(wtbl.abp_pnorm{w});
    t_next = median(wtbl.abp_pnorm{w});
    bin_dicron = zeros(length(bins)-1,1);
    bin_next = zeros(length(bins)-1,1);
    
    avg_bfi_ln = [];
    avg_abp_ln = [];
    std_bfi_ln = [];
    std_abp_ln = [];
    if w==1
        xline(bins(1))
        
    end
    for b = 2:nbins
        bfi_bin_ln = bfi_tnorm <= bins(b) & bfi_tnorm > bins(b-1);
        if t_dicron <= bins(b) & t_dicron > bins(b-1);
            bin_dicron(b-1) = 1;
        end
        if t_next <= bins(b) & t_next > bins(b-1);
            bin_next(b-1) = 1;
        end        
        avg_bfi_ln(b-1) = nanmean(bfi_pstack(bfi_bin_ln));
        std_bfi_ln(b-1) = nanstd(bfi_pstack(bfi_bin_ln));
        abp_bin_ln = abp_tnorm <= bins(b) & abp_tnorm > bins(b-1);
        avg_abp_ln(b-1) = nanmean(abp_pstack(abp_bin_ln));
        std_abp_ln(b-1) = nanstd(abp_pstack(abp_bin_ln));
       if w==1

        xline(bins(b))
       end
    end
    
    [fp,floc]=findpeaks(avg_bfi_ln);
    [pp,ploc]=findpeaks(avg_abp_ln);
    if isempty(floc)
        wtbl.psf_ln(w) = {NaN};
        wtbl.psf_loc_ln(w) = {NaN};
    else
        wtbl.psf_ln(w) = {fp(1)};
        wtbl.psf_loc_ln(w) = {floc(1)};
    end
        
    
    wtbl.psp_ln(w) = {pp(1)};
    
    wtbl.psp_loc_ln(w) = {ploc(1)};
    
    wtbl.binsize_ln(w) = {bin_size};
    wtbl.bins_ln(w) = {bins};
    wtbl.tbin_ln(w) = {bins(1:end-1) + bin_size/2};
    wtbl.t_dicron_ln(w) = {wtbl.tbin_ln{w}(logical(bin_dicron))};
    wtbl.loc_dicron_ln(w) = {find(bin_dicron)};
    wtbl.t_next_ln(w) = {wtbl.tbin_ln{w}(logical(bin_next))};
    wtbl.loc_next_ln(w) = {find(bin_next)};    
    wtbl.avg_bfi_ln(w) = {avg_bfi_ln};
    wtbl.std_bfi_ln(w) = {std_bfi_ln};
    wtbl.avg_abp_ln(w) = {avg_abp_ln};
    wtbl.std_abp_ln(w) = {std_abp_ln};
    wtbl.p_auc_ln(w) = {trapz(wtbl.tbin_ln{w},avg_abp_ln)};
    wtbl.f_auc_ln(w) = {trapz(wtbl.tbin_ln{w},avg_bfi_ln)};

end

gp = ~(cellfun(@isnan,wtbl.abp_avg_l));
wtbl.osf_ln = placehold;
wtbl.edf_ln = placehold;
wtbl.mf_ln = placehold;
wtbl.osp_ln = placehold;
wtbl.edp_ln = placehold;
wtbl.mp_ln = placehold;
wtbl.fph_ln = placehold;
wtbl.pph_ln = placehold;
wtbl.RI_ln = placehold;
wtbl.PI_ln = placehold;
wtbl.starts_ln = placehold;
wtbl.peaks_ln = placehold;

wtbl.osf_ln(gp) = num2cell(cellfun(@(x) x(1), wtbl.avg_bfi_ln(gp)));
wtbl.edf_ln(gp) = num2cell(cellfun(@(x) x(end), wtbl.avg_bfi_ln(gp)));
wtbl.mf_ln(gp) = num2cell(cellfun(@(x,y) x./(bin_size*length(y)), wtbl.f_auc_ln(gp),wtbl.tbin_ln(gp)));
%wtbl.mf_ln(gp) = num2cell(cellfun(@(y) (y(end)-y(1)),wtbl.avg_bfi_ln(gp)));

wtbl.fph_ln(gp) = num2cell(cell2mat(wtbl.psf_ln(gp)) - cell2mat(wtbl.edf_ln(gp)));
wtbl.PI_ln(gp) = num2cell(cell2mat(wtbl.fph_ln(gp))./cell2mat(wtbl.mf_ln(gp)));
wtbl.RI_ln(gp) = num2cell(cell2mat(wtbl.fph_ln(gp))./cell2mat(wtbl.psf_ln(gp)));

wtbl.osp_ln(gp) = num2cell(cellfun(@(x) x(1), wtbl.avg_abp_ln(gp)));
wtbl.edp_ln(gp) = num2cell(cellfun(@(x) x(end), wtbl.avg_abp_ln(gp)));
wtbl.mp_ln(gp) = num2cell(cellfun(@(x,y) x./(bin_size*length(y)), wtbl.p_auc_ln(gp),wtbl.tbin_ln(gp)));
wtbl.pph_ln(gp) = num2cell((cell2mat(wtbl.psp_ln(gp)) - cell2mat(wtbl.edp_ln(gp))));


cnames = categorical(wtbl.name);
names = categories(cnames);
binsdiff = diff(wtbl.bins_ln{1});
bin = binsdiff(1);

wtbl.starts_ln(gp) = num2cell(cellfun(@(x) mean(x(:,1)),wtbl.bfi_tnorm(gp)));
wtbl.peaks_ln(gp) = num2cell(cell2mat(wtbl.psf_loc_ln(gp)).*bin);
%% NOT LENGTH NORMALIZED

placehold = cell(height(wtbl),1);
wtbl.abp_t = placehold;
wtbl.bfi_t = placehold;
for w=1:height(wtbl)
    abp_tstack = wtbl.abp_tstack{w};
    bfi_tstack = wtbl.bfi_tstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    bfi_pstack = wtbl.bfi_pstack{w};

    abp_t = abp_tstack;
    bfi_t = bfi_tstack;
    wtbl.abp_t(w) = {abp_t};
    wtbl.bfi_t(w) = {bfi_t};
%     if w==1
%         figure
%         for p=1:size(abp_tstack,1)
%             plot(abp_t(p,:),normalize(abp_pstack(p,:),'range',[0 1]),'.-','Color','red','MarkerSize',10,'LineWidth',0.5)
%             hold on
%             plot(bfi_t(p,:),normalize(bfi_pstack(p,:),'range',[0 1]),'.-','color','blue','MarkerSize',10,'LineWidth',0.5)
%             graph_beautify(gca,30,0.5,1)
%             ylim([0 1.1])
%         end
%     end
end
%%
%Bin and average


wtbl.bins = placehold;
wtbl.tbin = placehold;
wtbl.avg_bfi = placehold;
wtbl.avg_abp = placehold;
wtbl.std_bfi = placehold;
wtbl.std_abp = placehold;
wtbl.psf = placehold;
wtbl.psp = placehold;
wtbl.psf_loc = placehold;
wtbl.psp_loc = placehold;
wtbl.p_auc = placehold;
wtbl.f_auc = placehold;
wtbl.binsize = placehold;


for w=1:height(wtbl)
    bins = -(bin_size/2):bin_size:(1+(bin_size/2));
    nbins = length(bins);

    bfi_t = wtbl.bfi_t{w};
    if isempty(bfi_t)
        continue
    end
    abp_t = wtbl.abp_t{w};
    bfi_pstack = wtbl.bfi_pstack{w};
    abp_pstack = wtbl.abp_pstack{w};
    bfi_clr = wtbl.clr{w}(:);
    state = wtbl.substate{w};
    sds_p = wtbl.sds(w);
    name = wtbl.name{w};
    wnum = wtbl.wnum(w);
        
    avg_bfi = [];
    avg_abp = [];
    std_bfi = [];
    std_abp = [];
    if w==1
        xline(bins(1))
        
    end
    % I disregard the first bin that just has BFI timepoints. This happens
    % because of sampling rate mismatching
    for b = 2:nbins 
        bfi_bin = bfi_t <= bins(b) & bfi_t > bins(b-1);
        avg_bfi(b-1) = nanmean(bfi_pstack(bfi_bin));
        std_bfi(b-1) = nanstd(bfi_pstack(bfi_bin));
        abp_bin = abp_t <= bins(b) & abp_t > bins(b-1);
        avg_abp(b-1) = nanmean(abp_pstack(abp_bin));
        std_abp(b-1) = nanstd(abp_pstack(abp_bin));
       if w==1

        xline(bins(b))
       end
    end
    
    [fp,floc]=findpeaks(avg_bfi);
    [pp,ploc]=findpeaks(avg_abp);
    if isempty(floc)
        wtbl.psf(w) = {NaN};
        wtbl.psf_loc(w) = {NaN};
    else
        wtbl.psf(w) = {fp(1)};
        wtbl.psf_loc(w) = {floc(1)};
    end
        
        
    wtbl.psp(w) = {pp(1)};
    
    wtbl.psp_loc(w) = {ploc(1)};
    
    wtbl.binsize(w) = {bin_size};
    wtbl.bins(w) = {bins};
    wtbl.tbin(w) = {bins(1:end-1) + bin_size/2};
    wtbl.avg_bfi(w) = {avg_bfi};
    wtbl.std_bfi(w) = {std_bfi};
    wtbl.avg_abp(w) = {avg_abp};
    wtbl.std_abp(w) = {std_abp};
    wtbl.p_auc(w) = {trapz(wtbl.tbin{w}(~isnan(avg_abp)),avg_abp(~isnan(avg_abp)))};
    wtbl.f_auc(w) = {trapz(wtbl.tbin{w}(~isnan(avg_bfi)),avg_bfi(~isnan(avg_bfi)))};

end

gp = ~(cellfun(@isnan,wtbl.abp_avg_l));

gpf = cellfun(@(x) ~isnan(x),wtbl.avg_bfi,'UniformOutput',false);
gpp = cellfun(@(x) ~isnan(x),wtbl.avg_abp,'UniformOutput',false)

wtbl.osf = placehold;
wtbl.edf = placehold;
wtbl.mf_nanmean = placehold;
wtbl.mf = placehold;
wtbl.osp = placehold;
wtbl.edp = placehold;
wtbl.mp_nanmean = placehold;
wtbl.mp = placehold;
wtbl.fph = placehold;
wtbl.pph = placehold;
wtbl.fpl = placehold;
wtbl.ppl = placehold;
wtbl.RI = placehold;
wtbl.PI = placehold;
wtbl.starts = placehold;
wtbl.peaks = placehold;

wtbl_gp = wtbl(gp,:);


for w=1:height(wtbl_gp)
tbin_bfi = wtbl_gp.tbin{w}(gpf{w,:});
tbin_abp = wtbl_gp.tbin{w}(gpp{w,:});
avg_bfi = wtbl_gp.avg_bfi{w}(gpf{w,:});
avg_abp = wtbl_gp.avg_abp{w}(gpp{w,:});
    
wtbl_gp.osf(w) = {avg_bfi(1)};
wtbl_gp.edf(w) = {avg_bfi(end)};
wtbl_gp.mf_nanmean(w) = {nanmean(avg_bfi)};
wtbl_gp.mf(w) = {wtbl_gp.f_auc{w}./(bin_size*length(wtbl_gp.tbin{w}))};
wtbl_gp.fpl(w) = {tbin_bfi(end) - tbin_bfi(1)};
%wtbl_gp.mf(w) = num2cell(cellfun(@(x,y) sum(x(:))./(y(end)-y(1)),wtbl_gp.avg_bfi(gp),wtbl_gp.tbin(gp)));
%wtbl.mf(gp) = num2cell(cellfun(@(y) (y(end)-y(1)),wtbl.avg_bfi(gp)));

wtbl_gp.fph(w) = {wtbl_gp.psf{w} - wtbl_gp.edf{w}};
wtbl_gp.ppl(w) = {tbin_abp(end) - tbin_abp(1)};
wtbl_gp.PI(w) = {wtbl_gp.fph{w}/wtbl_gp.mf{w}};
wtbl_gp.RI(w) = {wtbl_gp.fph{w}/wtbl_gp.psf{w}};

wtbl_gp.osp(w) = {avg_abp(1)};
wtbl_gp.edp(w) = {avg_abp(end)};
wtbl_gp.mp_nanmean(w) = {nanmean(avg_abp)};
wtbl_gp.mp(w) = {wtbl_gp.p_auc{w}./(bin_size*length(wtbl_gp.tbin{w}))};
wtbl_gp.pph(w) = {wtbl_gp.psp{w} - wtbl_gp.edp{w}};

end
wtbl(gp,:) = wtbl_gp;
wtbl.starts(gp) = num2cell(cellfun(@(x) mean(x(:,1)),wtbl_gp.bfi_t(gp)));
wtbl.peaks(gp) = num2cell(cell2mat(wtbl_gp.psf_loc(gp)).*bin);
%% DAT TBL

cnames = categorical(wtbl.name);
names = categories(cnames);
csds = unique(wtbl.sds);
cstate = categorical(wtbl.state);
states = flip(categories(cstate));

minutes = 3;
dat = {};

for s=1:length(names)
    dat{s,1} = char(names(s));
    dat{s,2} = wtbl(cnames==names(s),:);
    dat{s,3} = height(dat{s,2});
end

state = {'HC'};
statenames = {'hc'};
%fieldnames = {'name';'all';'length';'rs';'hc';'hc_bl';'hc_end';'rs_bl'};
fieldnames = {'name';'all';'length';'hc_bl';'hc_end'};

n = 2;

for s=1:size(dat,1)
    tbl_n = dat{s,2};
    %dat{s,3} = {height(dat{s,2})};

    for st = 1:length(states)

        tbl_sn = tbl_n(categorical(tbl_n.state)==state{st},:);
        % rs and hc
        %dat{s,3+st} = tbl_sn;
        if contains(states{st},'HC')
            % hc_end
            tbl_gs = tbl_sn(categorical(tbl_sn.substate)=='HCGasOn',:);
            wend = max(tbl_gs.wend);
            wstart = wend - (60*minutes);
            hc_tbln = tbl_gs(tbl_gs.wstart >= wstart,:);
            %dat{s,7} = hc_tbln;
            dat{s,5} = hc_tbln;
            fprintf('HCGasOn')
            
            % hc_bl
            gas_start = tbl_gs.wstart(1)
            
            hc_bl_n = tbl_sn(tbl_sn.wstart < gas_start,:);
            
            %dat{s,6} = hc_bl_n;
            dat{s,4} = hc_bl_n;
%         elseif contains(states{st},'Resting')
%             % rs_bl
%             rs_baseline = tbl_sn(categorical(tbl_sn.substate)=='RestingBaseline',:);
%             fprintf('rs')
%             dat{s,8} = rs_baseline;
        end
        
    end
end
%%
dat_tbl = cell2table(dat)
dat_tbl.Properties.VariableNames = fieldnames;  
end