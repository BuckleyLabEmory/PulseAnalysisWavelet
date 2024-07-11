function pulse_plots(tbl,id,p_clrs,ln,savepath)
%%
    sds = unique(tbl.sds);
    %find unique windows
    winds = unique(tbl.wnum);
    row_ids = {'Accepted','Rejected'};
    column_ids = {'SDS1','SDS2','ABP'};
    
    fig = figure('Units','Pixels','Position',[488,279,907,579]);
    
    for w=1:length(winds)
        wselect = tbl.wnum==winds(w);
        tbl_rw = tbl(wselect,:);
        
        if ~ln
        tpulse = tbl_rw.tbin{1}; % the same across pulses and windows
        
        pulses = {};
        pulses{1} = tbl_rw.avg_bfi{1}; % sds 1 pulse from window
        pulses{2} = tbl_rw.avg_bfi{2}; % sds2 pulse from window
        pulses{3} = tbl_rw.avg_abp{1}; % abp pulse from window

        pulse_height = {};
        pulse_height{1} = tbl_rw.fph{1};
        pulse_height{2} = tbl_rw.fph{2};
        pulse_height{3} = tbl_rw.pph{1};
        
        
        elseif ln
            tpulse = tbl_rw.tbin_ln{1}; % the same across pulses and windows

            pulses = {};
            pulses{1} = tbl_rw.avg_bfi_ln{1}; % sds 1 pulse from window
            pulses{2} = tbl_rw.avg_bfi_ln{2}; % sds2 pulse from window
            pulses{3} = tbl_rw.avg_abp_ln{1}; % abp pulse from window

            pulse_height = {};
            pulse_height{1} = tbl_rw.fph_ln{1};
            pulse_height{2} = tbl_rw.fph_ln{2};
            pulse_height{3} = tbl_rw.pph_ln{1};
        end
        
        status = tbl_rw.wp_accept{1};
        
        if status == 1
            plts = [1,2,3];
        elseif status == 0
            plts = [4,5,6];
        end
        
        for p=1:3
            subplot(2,3,plts(p))
            plot(tpulse,pulses{p},'color',p_clrs(p,:))
            title([column_ids{p}])
            graph_beautify(gca,12,1,1)
            hold on
        end

    end  
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    han.Title.Visible ='on';
    ylabel(han,'Rejected                                       Accepted');
    xlabel(han,'Length normalized time (s)');
    title(han,[{id},{''},{''}]);
    graph_beautify(gca,15,1,1)
    
    saveas(fig,[savepath filesep id '.png'],'png')
    savefig(fig,[savepath filesep id '.fig'])
    
    
end