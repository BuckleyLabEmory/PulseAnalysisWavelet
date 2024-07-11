function ptimeseries(dat_tbl,sstates,clr_custom,save_path)
%% A figure for each entry in dat_tbl with sds1,2 in both states
id1 = cellfun(@(x) strsplit(x,' '),dat_tbl.names,'UniformOutput',false);
ids = cellfun(@(x) x{1},id1,'UniformOutput',false);

for s=1:height(dat_tbl)
    c = 1;
    figure('Units','Pixels','Position',[162,439,1338,419]) 
    hold on
    x = 1;
    for sd = 1:2  
        sd;
        for st = 1:length(sstates)
            sstates{st};
            tbl_x = dat_tbl{s,sstates{st}};
            tbl_x = tbl_x{1};

        
            id = ids{s};
            tbl_n = tbl_x(tbl_x.sds==sd,:);
            x = [x;length(horzcat(tbl_n.avg_bfi_ln{:}))];
            xp = (x(st):(x(st)+x(st+1)-1))*(1/20);
            yp = horzcat(tbl_n.avg_bfi{:});
            size(xp);
            size(yp);
            plot(xp,yp,'Color',clr_custom(c,:))
            graph_beautify(gca,20,2,2)
            title(id)
            c = c+1;
        end
       legend({'RS1','HC1','RS2','HC2'},'Location','North','Orientation','Horizontal')

    end
    saveas(gca,[save_path filesep id ' ptimeseries.png'],'png')
    savefig([save_path filesep id ' ptimeseries.fig'])
end
end