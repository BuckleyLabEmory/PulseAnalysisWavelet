function [dat_tbl_qc,stats_qc,pstats_qc] = dat_tbl_analysis(dat_tbl,savepath)
%% NEED TO ADD ONSET/PEAK/END MARKS TO PULSES TO FIGURE OUT WHAT IS GOING ON

%%

%%

% if there is a spike after 0.6 larger than the peak
% if pulse peak height is negative
try 
    name = dat_tbl.name;
catch
    name = dat_tbl.names;
end
names = categorical(name);
% IDs no identifiers
id1 = cellfun(@(x) strsplit(x,' '),name,'UniformOutput',false);
ids = cellfun(@(x) x{1},id1,'UniformOutput',false);
% ids = cellfun(@(x) x{3}, id1,'UniformOutput',false);
sds_i = [1 2];
states = {'hc_bl','hc_end'};
state_names = {'Resting','Hypercapnia'};
clrs = lines(4);
p_clrs = [{[clrs(1:2,:);[0.6 0 0]]},{[clrs(3:4,:);[0.6 0 0]]}];

% initialize the quality control table
dat_tbl_qc = table();
dat_tbl_qc.names = name;
for st=1:length(states)
   dat_tbl_qc.(states{st}) = cell(height(dat_tbl_qc),1); 
end


for n=1:length(names)
   for st=1:length(states)
       tbl = dat_tbl.(states{st}){n};
       id = [ids{n} ' ' state_names{st}];
       tbl_qc = pulse_tests(tbl);   
       %pulse_plots_ln(tbl_qc,[id ' length norm' ],p_clrs{st},savepath)
       pulse_plots(tbl_qc,[id ' length norm' ],p_clrs{st},1,savepath)
       pstatus = cell2mat(tbl_qc.wp_accept) == 1;
       fprintf(['Removing ' num2str(length(find(~pstatus))) ' pulses from ' id '\n'])
       dat_tbl_qc.(states{st}){n} = tbl_qc(pstatus,:);
       fprintf(['Lengths dat_tbl for ' id ' is ' num2str(height(dat_tbl.(states{st}){n})) '\n'...
           'Lengths dat_tbl_qc for ' id ' is ' num2str(height(dat_tbl_qc.(states{st}){n})) '\n'])
   end    
end

%close all
%%

params = {
            'onset systolic flow',      'osf_ln';...
            'peak systolic flow',       'psf_ln';...
            'end diastolic flow',       'edf_ln';...
            'flow peak height',         'fph_ln';...
            
            'onset systolic pressure',  'osp_ln';...
            'peak systolic pressure',   'psp_ln';...
            'end diastolic pressure',   'edp_ln';...            
            'pressure peak height',     'pph_ln';...
            
            'flow AUC',                 'f_auc_ln';...
            'pressure AUC',             'p_auc_ln';...
            

            'mean flow',                'mf_ln';...
            'mean pressure',            'mp_ln';...
            
            'flow resistive index',     'RI_ln';...
            'flow pulsatility index',   'PI_ln';...

            'Fhr/<F>',                      'Fratio';...
            'Phr/<P>',                      'Pratio';...
            'critical closing pressure',    'crcpmich';...
            'cerebral perfusion pressure',  'cppmich';...

      }
  
stats_qc = generate_sub_stats(params,dat_tbl_qc,states);

if height(dat_tbl) > 1
    pstats_qc = generate_pstats(stats_qc);


    psnames = {'hcbl1','hcend1','hcbl2','hcend2'};
    for pn=1:length(psnames)
       pstats_qc.([psnames{pn} ' mean']) = cellfun(@mean,pstats_qc.(psnames{pn})) ;
       pstats_qc.([psnames{pn} ' median']) = cellfun(@median,pstats_qc.(psnames{pn})) ;
       pstats_qc.([psnames{pn} ' std']) = cellfun(@std,pstats_qc.(psnames{pn})) ;

    end

%%

fun = @(m)sRGB_to_OSAUCS(m,true,true);
rgb = maxdistcolor(height(dat_tbl_qc),@sRGB_to_CIELab,'Cmin', 0.5, 'Cmax',0.6, 'Lmin',0.45, 'Lmax',0.8);
colors = rgb;

for p=1:length(params)
    figure_id = {};
    ax = {};
    
    param_i = params{p,1};
    p_select = categorical(pstats_qc.param_name) == param_i;
   %%%%%%%%%%%%%%%
   % BL1 vs. HC1 %
   %%%%%%%%%%%%%%%

   figs{1} = figure
   dat1 = pstats_qc.hcbl1{p_select}./pstats_qc.hcbl1{p_select};
   dat2 = pstats_qc.hcend1{p_select}./pstats_qc.hcbl1{p_select};
   name1 = "Baseline";
   name2 = "Hypercapnia";
   tt = [params{p,1}, {['% Change SDS1 p = ' num2str(pstats_qc.hcbl1_hcend1{p_select})]}];

   dat_plot(ids,dat1,dat2,name1,name2,tt,rgb)
   
   figure_id{1} = [params{p,2} ' % Change SDS1'];
   g1 = gca;
   l1 = get(g1,'ylim')
  %ylim([0 2.5])

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%
   % BL2 vs. HC2 %
   %%%%%%%%%%%%%%%

   figs{2} = figure
   dat1 = pstats_qc.hcbl2{p_select}./pstats_qc.hcbl2{p_select};
   dat2 = pstats_qc.hcend2{p_select}./pstats_qc.hcbl2{p_select};
   name1 = "Baseline";
   name2 = "Hypercapnia";
   tt = [params{p,1}, {['% Change SDS2 p = ' num2str(pstats_qc.hcbl2_hcend2{p_select})]}];

   dat_plot(ids,dat1,dat2,name1,name2,tt,rgb)
   
   figure_id{2} = [params{p,2} ' % Change SDS2'];
   g2 = gca;
   l2 = get(g2,'ylim')
   %ylim([0 2.5])
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
   lset_max = max([l1;l2]);
   lset_min = min([l1;l2]);
   set(g1,'ylim',[lset_min(1) lset_max(2)]);
   set(g2,'ylim',[lset_min(1) lset_max(2)]);
   %set(g1,'ylim',[0.5 2.5]) 
   ax = {g1,g2};
   for x=1:length(ax)
       fid = figure_id{x};
       saveas(ax{x},[savepath filesep fid '.png'],'png')
       figure(figs{x})
       savefig(gcf,[savepath filesep fid '.fig'])
   end

end

%%

%%

 

for p=1:length(params)
    figure_id = {};
    ax = {};
    
    param_i = params{p,1};
    p_select = categorical(pstats_qc.param_name) == param_i;
   %%%%%%%%%%%%%%%
   % BL1 vs. HC1 %
   %%%%%%%%%%%%%%%

   figs{1} = figure
   dat1 = pstats_qc.hcbl1{p_select};
   dat2 = pstats_qc.hcend1{p_select};
   name1 = "Baseline";
   name2 = "Hypercapnia";
   tt = [params{p,1}, {['Change SDS1 p = ' num2str(pstats_qc.hcbl1_hcend1{p_select})]}];

   dat_plot(ids,dat1,dat2,name1,name2,tt,rgb)
   
   figure_id{1} = [params{p,2} ' Change SDS1'];
   g1 = gca;
   l1 = get(g1,'ylim')
  %ylim([0 2.5])

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%
   % BL2 vs. HC2 %
   %%%%%%%%%%%%%%%

   figs{2} = figure
   dat1 = pstats_qc.hcbl2{p_select};
   dat2 = pstats_qc.hcend2{p_select};
   name1 = "Baseline";
   name2 = "Hypercapnia";
   tt = [params{p,1}, {[' Change SDS2 p = ' num2str(pstats_qc.hcbl2_hcend2{p_select})]}];

   dat_plot(ids,dat1,dat2,name1,name2,tt,rgb)
   
   figure_id{2} = [params{p,2} ' Change SDS2'];
   g2 = gca;
   l2 = get(g2,'ylim')
   %ylim([0 2.5])
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
   lset_max = max([l1;l2]);
   lset_min = min([l1;l2]);
   set(g1,'ylim',[lset_min(1) lset_max(2)]);
   set(g2,'ylim',[lset_min(1) lset_max(2)]);
   %set(g1,'ylim',[0.5 2.5]) 
   ax = {g1,g2};
   for x=1:length(ax)
       fid = figure_id{x};
       saveas(ax{x},[savepath filesep fid '.png'],'png')
       figure(figs{x})
       savefig(gcf,[savepath filesep fid '.fig'])
   end

end



%%
sub_params = {'osf_ln','psf_ln','edf_ln','fph_ln','osp_ln','psp_ln','edp_ln','pph_ln'};
fph_stats = stats_qc(any(categorical(stats_qc.param) == sub_params,2),:);
fph_pstats = pstats_qc(any(categorical(pstats_qc.param) == sub_params,2),:);
peak_height_plots(fph_stats,fph_pstats,savepath)

else
    pstats_qc = [];
end
%%
plot_state = 'hc_bl'
clrs = lines(4);
clr_s = clrs(1:2,:);
state_pulse_plot(dat_tbl_qc,plot_state,clr_s,savepath)
%%
clrs = lines(4);
clr_custom = [clrs(1,:);clrs(3,:);clrs(2,:);clrs(4,:)];
sstates = {'hc_bl','hc_end'};
ptimeseries(dat_tbl_qc,sstates,clr_custom,savepath)

% %% window figures
% swindow = dat_tbl_qc.hc_bl{2}(1,:);
% wprocess_plots(swindow,savepath)
end