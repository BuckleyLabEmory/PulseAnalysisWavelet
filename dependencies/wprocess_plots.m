function wprocess_plots(swindow,save_path)
swindow
%%
figure('Units','Inches','Position',[4,4,9,7])
bfi = swindow.bfi_pstack{1};
abp = swindow.abp_pstack{1};
tbfi = swindow.bfi_tnorm{1};
tabp = swindow.abp_tnorm{1};
bins = swindow.bins_ln{1};


for b=1:length(bins)
    xline(bins(b),'LineWidth',1.5,'Color','black')
end

for p=1:size(bfi,1)/2
    patchline(tabp(p,:),normalize(abp(p,:),'range',[0 1]),'EdgeColor',[0.6 0 0],'EdgeAlpha',0.7,'linewidth',1)%,'b.-')
      hold on

  plot(tabp(p,:),normalize(abp(p,:),'range',[0 1]),'.','color',[0.6 0 0],'MarkerSize',7)
  patchline(tbfi(p,:),normalize(bfi(p,:),'range',[0 1]),'EdgeColor',[0    0.4470    0.7410],'EdgeAlpha',1,'linewidth',1)%,'b.-')
  plot(tbfi(p,:),normalize(bfi(p,:),'range',[0 1]),'.','color',[0    0.4470    0.7410], 'MarkerSize',11)
  

end

xlim([-0.05 1])
graph_beautify(gca,1,2,1)

yticks('')
xticks('')
ylim([-0.2 1.2])
set(gca,'visible','off')

%%
figure('Units','Inches','Position',[4,4,5,6])
plot(swindow.tbin_ln{:},normalize(swindow.avg_bfi_ln{:},'Range',[0 1]),'-','color',[0    0.4470    0.7410])
hold on
plot(swindow.tbin_ln{:},normalize(swindow.avg_abp_ln{:},'Range',[0 1]),'-','color',[0.6 0 0])

xlim([-0.05 1])
graph_beautify(gca,1,4,1)

yticks('')
xticks('')
ylim([-0 1.1])
set(gca,'box','off')
set(gca,'visible','off')
%%
figure('Units','Inches','Position',[4,4,5,2])
wt = swindow.wt{:};
wbfi = swindow.wbfi{:};
wabp = swindow.wabp{:};
plot(wt,normalize(wbfi,'range',[0 1]),'.-','MarkerSize',5)
hold on
%plot(wn.wt_125hz{:},normalize(wn.wabp_125hz{:},'range',[0 1]),'.','Color',clrs(3,:))
plot(wt,normalize(wabp,'range',[0 1]),'.-','Color',[0.6 0 0],'MarkerSize',5)
graph_beautify(gca,10,0.5,1)

%yticks('')
%xticks('')
ylim([0 0.8])
xlim([1090 1095])
%set(gca,'visible','off')

end