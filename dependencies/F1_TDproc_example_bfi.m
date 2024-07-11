function F1_TDproc_example_bfi(dat_tbl, state, save_path,subject_id)
% 12/30/2022 - v6 updates
% First commit version taken directly from: /Users/taraurner/OneDrive - Emory University/Pulsatility_Paper/FIGURES AND DATA/1_F1 waveform analysis/Archive/Fig_TD_processing_f31.m

%clear
%close all

%load('/Volumes/labs/buckley-lab/Projects/Waveform_analysis/0_Papers/2023_PulsatilityPaper/40_PULSEANALYSIS_final/dat_tbl_qc.mat')
%save_path = '/Users/taraurner/OneDrive - Emory University/Pulsatility_Paper/FIGURES AND DATA/1_F1 waveform analysis/'
%%
close all

% fig_p1 = 100;
% fig_p2 = 100;
fig_h = 800;
fig_w = 800;
fig_p1 = 1.5
fig_p2 = 2,
fig_p3 = 8
fig_p4 = 8
fig_size = [fig_p1,fig_p2,fig_w,fig_h];
ax_size = [0 0 1 1];
fontsize = 25;
box_width = 5;

bfi_color = [0.10,0.25,1.00];
%abp_color = [0.6 0 0];
abp_color = [1 0.16 0]
%%%
% swindow = dat_tbl_qc.hc_bl{2}(37,:); subject 2 = PL, window 37

%swindow = dat_tbl_qc.hc_bl{5}(11,:); % subject 5 = SS, window 11

%state_name = dat_tbl.(state){1}.state{1};
all_dat = vertcat(dat_tbl.(state){:});
% IDs no identifiers
ids = string(cellfun(@(x) strsplit(x,'~'),all_dat.name,'UniformOutput',false));
subs = categorical(ids);
sub_i = categories(subs);

swindow = all_dat(subs == sub_i(1),:);

bfi = swindow.bfi_pstack{1};
abp = swindow.abp_pstack{1};
tbfi = swindow.bfi_tnorm{1};
tabp = swindow.abp_tnorm{1};
bins = swindow.bins_ln{1};
bfi_avg = swindow.avg_bfi_ln{1};
abp_avg = swindow.avg_abp_ln{1};
bfi_std = swindow.std_bfi_ln{1};
abp_std = swindow.std_abp_ln{1};
tbins = swindow.tbin_ln{1};

wt = swindow.wt{:};
wt_125hz = swindow.wt_125hz{:};
wbfi = swindow.wbfi{:};
wabp = swindow.wabp_125hz{:};

%%

%%%%%%%%%%%%%%
% BIN PLOTS
%%%%%%%%%%%%%%
close all

h_factor = 1;
w_factor = 1;

edge_alpha_bfi = 0.7;
patchlinewidth_bfi = 1.5;
edge_alpha_abp = 0.5;
patchlinewidth_abp = 1;
marker_size_bfi = 25;
marker_size_abp = 6;

line_size_avg = 4;
marker_size_avg = 50;

xt = [0, 0.25 0.5 .75 1];
xl = [-0.05 1.0];
%yl = [-0.1 1.5]
%yt = [0, 0.25 0.5 .75 1 1.25 1.5]

bin_color = [0.3 0.3 0.3]
bin_linewidth = 1.1;

% Overlay figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('Units','Inches','Position',[fig_p1,fig_p2,fig_p3*w_factor,fig_p4*h_factor])
ax = axes;
fname = strcat(subject_id,'TD_binned_overlay');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box on
box_w = 8

for b=1:length(bins)
    xline(bins(b),'LineWidth',bin_linewidth,'Color',bin_color)
    hold on
end

% abp_pmax = max(max(abp));
% abp_pmin = min(min(abp));
% bfi_pmax = max(max(bfi));
% bfi_pmin = min(min(bfi));

abp_pmax = max(abp_avg);
abp_pmin = min(abp_avg);
bfi_pmax = max(bfi_avg);
bfi_pmin = min(bfi_avg);

% BFI
yyaxis left
for p=1:size(bfi,1)
  tp_bfi = tbfi(p,:);
  tp_use = tp_bfi > bins(1);
  %p_bfi = (bfi(p,:) - bfi_pmin)./(bfi_pmax - bfi_pmin);
  p_bfi = bfi(p,:);
  patchline(tp_bfi(tp_use),p_bfi(tp_use),'LineStyle','-','EdgeColor',bfi_color,'EdgeAlpha',edge_alpha_bfi,'linewidth',patchlinewidth_bfi)%,'b.-')
  scatter(tp_bfi(tp_use),p_bfi(tp_use),marker_size_bfi,bfi_color,"filled")
end
%ylim([0 bfi_pmax+bfi_pmax*0.1])
set(gca,'YColor',bfi_color)
yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([0E-9 25E-9])

% ABP
yyaxis right
hold on
for p=1:size(bfi,1)
  tp_abp = tabp(p,:);
  %p_abp = (abp(p,:) - abp_pmin)./(abp_pmax - abp_pmin); % normalizing the pulse stack on [0 1]
  p_abp = abp(p,:);
  patchline(tp_abp,p_abp,'LineStyle','-','EdgeColor',abp_color,'EdgeAlpha',edge_alpha_abp,'linewidth',patchlinewidth_abp)%,'b.-')
  scatter(tp_abp,p_abp,marker_size_abp,abp_color,"filled")
end
set(gca,'YColor',abp_color)
yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([-5E-9 15E-9])

%plot(tbins,(abp_avg - abp_pmin)./abp_pmax,'.-','color',abp_color,'MarkerSize',marker_size_avg,'LineWidth',line_size_avg)

%plot(tbins,(bfi_avg - bfi_pmin)./bfi_pmax,'.-','color',bfi_color,'MarkerSize',marker_size_avg,'LineWidth',line_size_avg)


ax = gca;
%box(ax,'on')
set(ax,'linewidth',box_w)
set(ax,'FontSize',fontsize)
xticks(xt)
%yticks(yt)
xlim(xl)
%ylim(yl)
set(gca,'TickDir','out');

%%%%%
savefig(f,[save_path fname '.fig'])
yyaxis left
xticklabels('')
yticklabels('')
yyaxis right
xticklabels('')
yticklabels('')
set(ax,'OuterPosition',ax_size);
exportgraphics(f,[save_path fname '.png'],'Resolution',300)
%%%%%

%%

%%%%%%%%%%%%%%
% WINDOW PLOT
%%%%%%%%%%%%%%

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I can use fig_p1, p2, etc.. if I want to put the figures side by size in
% a row
%f = figure('Units','Inches','Position',[fig_p1,fig_p2,fig_w*4,fig_h]) 

%The below sizing was for the layout option of example window in top row of
%figure by itself
ws_size = [fig_p1, fig_p2, fig_p3*4,fig_p4*(2/3)]
ws_ylim = [0.27,2.32]
ws_yt = [0:0.5:2.5]
ws_xt = [0:5:15]
f = figure('Units','Inches','Position',ws_size)

ax = axes;

fname = strcat(subject_id,'_TD_sample_window');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box on
box_w = 4;


linewidth = 2;
markersize_abp = 6;
markersize_bfi = 20;
markers = 1;
wt_norm = wt - wt(1);
wt_125hz_norm = wt_125hz - wt_125hz(1);
xlim_seconds = [0 15]

yyaxis right
plot(ax,wt_125hz_norm,wabp,'-','Color',abp_color,'LineWidth',linewidth)
hold on
if markers
    plot(ax,wt_125hz_norm,wabp,'.','Color',abp_color,'MarkerSize',markersize_abp)
end
set(gca,'YColor',abp_color)
yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([-5E-9 15E-9])

yyaxis left
plot(ax,wt_norm,wbfi,'-','LineWidth',linewidth,'Color',bfi_color)
hold on
if markers
     plot(ax,wt_norm,wbfi,'.','Color',bfi_color,'MarkerSize',markersize_bfi)
end
set(gca,'YColor',bfi_color)

yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([0E-9 25E-9])

%box(ax,'off')
xlim(xlim_seconds)
%ylim(ws_ylim)
set(gca,'TickDir','out');
%yticks(ws_yt)
xticks(ws_xt)

%ax.XTick = wt(1):5:wt(end)
%box(ax,'on')
set(ax,'linewidth',box_w)
set(ax,'FontSize',fontsize)
pad = 0.07
%set(ax,'Position',)
set(ax,'OuterPosition',ax_size);
%set(ax,'Position',[0+pad, 0+pad, 1-pad,1-pad]);
set(ax,'Position',[0.05 0.05 0.9 0.9]);


%%%%%
savefig(f,[save_path fname '.fig'])
yyaxis left
xticklabels('')
yticklabels('')
yyaxis right
xticklabels('')
yticklabels('')
exportgraphics(f,[save_path fname '.png'],'Resolution',300)
%%%%%

%close all
%%

%%%%%%%%%%%%%%
% AVERAGE PLOT
%%%%%%%%%%%%%%

close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = figure('Units','Inches','Position',[fig_p1,fig_p2,fig_p3,fig_p4])
ax = axes;
fname = strcat(subject_id,'_normalized_pulses_std');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box on
box_w = 8;

linewidth = 4;
std_alpha = 0.2;
fontsize = 25;
xt = [0, 0.25 0.5 .75 1];
% yl is set in the binning plot so they will be the same
hold on
%%%

%%%%% ABP
pt = tbins;
%pm = normalize(abp_avg,'range',[0 1]);
%pm = (abp_avg - abp_pmin)./(abp_pmax - abp_pmin);
%pstd = abp_std./mean(abp_avg);
pm = abp_avg;
pstd = abp_std;
clr = abp_color;

curve_1 = (pm + pstd);
curve_2 = (pm - pstd);
x = pt;    

yyaxis right
plot(pt,pm,'color',clr,'LineWidth',linewidth)
hold on
x2 = [x, fliplr(x)];
inBetween = [curve_1, fliplr(curve_2)];
fill(x2,inBetween,'-','FaceColor',clr,'FaceAlpha',std_alpha,'EdgeColor',clr,'EdgeAlpha',std_alpha);
set(gca,'YColor',clr)
yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([-5E-9 15E-9])
%%%%%


%%%%% BFI
pt = tbins;
%pm = normalize(bfi_avg,'range',[0 1]);
% pm = (bfi_avg - bfi_pmin)./(bfi_pmax - bfi_pmin);
% pstd = bfi_std./mean(bfi_avg);
pm = bfi_avg;
pstd = bfi_std;
clr = bfi_color;

curve_1 = (pm + pstd);
curve_2 = (pm - pstd);
x = pt;    

yyaxis left
plot(pt,pm,'color',clr,'LineWidth',linewidth)
hold on
x2 = [x, fliplr(x)];
inBetween = [curve_1, fliplr(curve_2)];
fill(x2,inBetween,'-','FaceColor',clr,'FaceAlpha',std_alpha,'EdgeColor',clr,'EdgeAlpha',std_alpha);
set(gca,'YColor',clr)
yticks([5E-9 10E-9 15E-9 20E-9 25E-9])
%ylim([0E-9 25E-9])
%%%%%


ax = gca;
%ylim(yl)
%yticks(yt)
%box(ax,'off')
set(ax,'linewidth',box_w)
set(ax,'FontSize',fontsize)
set(gca,'TickDir','out');

%%%%%
savefig(f,[save_path fname '.fig'])
xticks(xt)
yyaxis left
xticklabels('')
yticklabels('')
yyaxis right
xticklabels('')
yticklabels('')
set(ax,'OuterPosition',ax_size);
exportgraphics(f,[save_path fname '.png'],'Resolution',300)
%%%%%

close all