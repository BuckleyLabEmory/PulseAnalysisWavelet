function peak_height_plots(stats_tbl,pstats_tbl,save_path)
stats_tbl

% IDs no identifiers
id1 = cellfun(@(x) strsplit(x,' '),stats_tbl.name,'UniformOutput',false);
ids = cellfun(@(x) x{3}, id1,'UniformOutput',false);
stats_tbl.ids = ids;
subs = categorical(ids);
sub_i = categories(subs);

states = categorical(stats_tbl.state);
state_i = categories(states);

clrs = lines(4);

ymax = [1.5,2,1.5,1.5,2,1.5,2.5,2.5,2];
%%
% SDS 1
sd = 1;
bl1 = pstats_tbl.hcbl1{categorical(pstats_tbl.param) == 'fph_ln'};
hc1 = pstats_tbl.hcend1{categorical(pstats_tbl.param) == 'fph_ln'};
dp1 = hc1./bl1;
inc = dp1;
bl_avg = bl1;
bl_name = 'hcbl1';
sclrs = [clrs(1,:);clrs(3,:)]


figure
n_subs = length(sub_i);
nrows = 3;
ncols = ceil(n_subs/nrows);

pplot

%%
% SDS 2
sd = 2;
bl2 = pstats_tbl.hcbl2{categorical(pstats_tbl.param) == 'fph_ln'};
hc2 = pstats_tbl.hcend2{categorical(pstats_tbl.param) == 'fph_ln'};
dp2 = hc2./bl2;

bl_avg = bl2;
inc = dp2;
bl_name = 'hcbl2';
sclrs = [clrs(2,:);clrs(4,:)];

pplot

%%
signrank(dp1,dp2)
% ABP
sd = 1;
bl2 = pstats_tbl.hcbl2{categorical(pstats_tbl.param) == 'pph_ln'};
hc2 = pstats_tbl.hcend2{categorical(pstats_tbl.param) == 'pph_ln'};
dp2 = hc2./bl2;

bl_avg = bl2;
inc = dp2;
bl_name = 'hcbl2';
sclrs = [clrs(2,:);clrs(4,:)];


pplot_abp

function pplot
    figure
    unsorted_ids = pstats_tbl.ids{1};
    [~,sort_order] = sort(unsorted_ids);

    osf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'osf_ln'};
    edf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'edf_ln'};
    psf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'psf_ln'};
    fph_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'fph_ln'};

    osf = osf_us(sort_order);
    edf = edf_us(sort_order);
    psf = psf_us(sort_order);
    fph = fph_us(sort_order);
    inc_s = inc(sort_order);

    for n=1:length(sub_i)
        for st = 1:length(state_i)
            subplot(nrows,ncols,n)
            tbl = stats_tbl(subs == sub_i(n) & stats_tbl.sds == sd & states == state_i(st),:);
            ln_pulses = tbl.ln_pulses{1};

            for p = 1:height(ln_pulses)
                pulse = ln_pulses.avg_bfi_ln{p};
                %plot(ln_pulses.tbin_ln{p},(pulse - pulse(end))./fph(n),'color',sclrs(st,:))
                %plot(ln_pulses.tbin_ln{p},pulse,'color',sclrs(st,:))
                plot(ln_pulses.tbin_ln{p},pulse,'color',sclrs(st,:))
                
                hold on
                %plot(ln_pulses.tbin_ln{p}(1),osf{p},'o')
                %plot(ln_pulses.tbin_ln{p}(end),edf{p},'o')
                %yline(psf{p})
            end
        end
        %yline(1,'linewidth',2,'color',sclrs(1,:))
        %yline(inc_s(n),'linewidth',2,'color',sclrs(2,:))

        %ylim([-0.5 ymax(n)])
       id = [sub_i{n} ' ' bl_name ' bfi ']
       title(id)
            
    end 
    

    saveas(gca,[save_path filesep bl_name ' pulse heights.png'],'png')
    savefig(gcf,[save_path filesep bl_name ' pulse heights.fig'])    
     
%     %%
%     figure
%     for n=1:length(sub_i)
%         for st = 1:length(state_i)
%             subplot(nrows,ncols,n)
%             tbl = stats_tbl(subs == sub_i(n) & stats_tbl.sds == sd & states == state_i(st),:);
%             ln_pulses = tbl.ln_pulses{1};
%             pt = ln_pulses.tbin_ln{1}
%             pulses = vertcat(ln_pulses.avg_bfi_ln{:});
%             %pm = median((pulses-pulses(:,end))./fph(n),1);
%             %pm_norm = normalize(pm,'range',[0 1]);
%             pm = median(pulses,1)
%             pm_norm = normalize(pm,'range',[0 1]);
%             %pstd = std((pulses-pulses(:,1))./fph(n));
%             pstd = std(pulses)./pm;
%             clr = sclrs(st,:);
% 
%             pulse_fillplot(pt,pm_norm,pstd,clr);
%             %ylim([-0.5 ymax(n)])
%             xticklabels('')
%             yticklabels('')
%             hold on
%             graph_beautify(gca,10,1,1)           
%         end
%     title([sub_i(n) ' ' bl_name])
% 
%     end
    
        %%
    figure
    for n=1:length(sub_i)
        figure
        for st = 1:length(state_i)
            tbl = stats_tbl(subs == sub_i(n) & stats_tbl.sds == sd & states == state_i(st),:);
            ln_pulses = tbl.ln_pulses{1};
            pt = ln_pulses.tbin_ln{1}
            pulses = vertcat(ln_pulses.avg_bfi_ln{:});
            %pm = median((pulses-pulses(:,end))./fph(n),1);
            %pm_norm = normalize(pm,'range',[0 1]);
            %pm = median(pulses,1)
            %pm_norm = normalize(pm,'range',[0 1]);
            %pstd = std((pulses-pulses(:,1))./fph(n));
            %pstd = std(pulses)./pm;
            pm = median(pulses,1);
            pstd = std(pulses);
            clr = sclrs(st,:);
                title(sub_i(n))

            pulse_fillplot(pt,pm,pstd,clr);
            %ylim([-0.5 ymax(n)])
            hold on
            graph_beautify(gca,10,1,1)           
        end

    id = [sub_i{n} ' ' bl_name ' bfi ']
    title(id)

    saveas(gca,[save_path filesep id ' pulse avgs.png'],'png')
    savefig(gcf,[save_path filesep id ' pulse avgs.fig'])  
    end
    

end


function pplot_abp
    figure
    unsorted_ids = pstats_tbl.ids{1};
    [~,sort_order] = sort(unsorted_ids);

    osf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'osp_ln'};
    edf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'edp_ln'};
    psf_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'psp_ln'};
    fph_us = pstats_tbl.(bl_name){categorical(pstats_tbl.param) == 'pph_ln'};

    osf = osf_us(sort_order);
    edf = edf_us(sort_order);
    psf = psf_us(sort_order);
    fph = fph_us(sort_order);
    inc_s = inc(sort_order);

    for n=1:length(sub_i)
        for st = 1:length(state_i)
            subplot(nrows,ncols,n)
            tbl = stats_tbl(subs == sub_i(n) & stats_tbl.sds == sd & states == state_i(st),:);
            ln_pulses = tbl.ln_pulses{1};

            for p = 1:height(ln_pulses)

                pulse = ln_pulses.avg_abp_ln{p};
                %plot(ln_pulses.tbin_ln{p},(pulse - pulse(end))./fph(n),'color',sclrs(st,:))
                plot(ln_pulses.tbin_ln{p},pulse,'color',sclrs(st,:))
                hold on

                %plot(ln_pulses.tbin_ln{p}(1),osf{p},'o')
                %plot(ln_pulses.tbin_ln{p}(end),edf{p},'o')
                %yline(psf{p})
            end
        end
        %yline(1,'linewidth',2,'color',sclrs(1,:))
        %yline(inc_s(n),'linewidth',2,'color',sclrs(2,:))

        %ylim([-0.5 ymax(n)])

    id = [sub_i{n} ' ' bl_name ' abp ']
    title(id)
    end  
    


    saveas(gca,[save_path filesep 'abp pulse heights.png'],'png')
    savefig(gcf,[save_path filesep id 'abp pulse heights.fig'])  
end
end