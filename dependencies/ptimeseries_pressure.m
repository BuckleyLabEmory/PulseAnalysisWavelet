function ptimeseries_pressure(dat_tbl,sstates,clr_custom,save_path)
%% Two figures for each entry in dat_tbl with normalized abp and one sds in both states
id1 = cellfun(@(x) strsplit(x,' '),dat_tbl.names,'UniformOutput',false);
ids = cellfun(@(x) x{1},id1,'UniformOutput',false);

for s=1%:height(dat_tbl)
    c = 1;
    x = 1;
    
    for sd = 1:2  
        figure('Units','Pixels','Position',[162,439,1338,419]) 
        hold on
        sd;
        for st = 1:length(sstates)
            sstates{st};
            tbl_x = dat_tbl{s,sstates{st}};
            tbl_x = tbl_x{1};
            id = ids{s};
            tbl_n = tbl_x(tbl_x.sds==sd,:);
            ptbl_n = tbl_x(tbl_x.sds==2,:);
            % ABP
                ax = [x;length(horzcat(ptbl_n.avg_abp_ln{:}))];
                axp = (ax(st):(ax(st)+ax(st+1)-1))*(1/20);
                ayp = horzcat(ptbl_n.avg_abp{:});
                ayp_norm = ayp ./ mean(ayp); 
                size(axp);
                size(ayp);
                plot(axp,ayp_norm,'Color','red')
                graph_beautify(gca,20,2,2) 
                hold on
            % BFI
            x = [x;length(horzcat(tbl_n.avg_bfi_ln{:}))];
            xp = (x(st):(x(st)+x(st+1)-1))*(1/20);
            yp = horzcat(tbl_n.avg_bfi{:});
            yp_norm = yp ./ nanmean(yp);
            size(xp);
            size(yp);
            plot(xp,yp_norm,'Color',clr_custom(c,:))
            graph_beautify(gca,20,2,2)
            title([id ' each state normalized to mean '])
            c = c+1;
        end
       legend({'RS1','HC1','RS2','HC2'},'Location','North','Orientation','Horizontal')

    end
    saveas(gca,[save_path filesep id ' ptimeseries.png'],'png')
    savefig([save_path filesep id ' ptimeseries.fig'])
end
end