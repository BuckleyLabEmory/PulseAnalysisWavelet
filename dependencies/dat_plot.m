function dat_plot(names,dat1,dat2,name1,name2,tt,rgb)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   dat = [dat1;dat2];
   group = [repmat(name1,length(dat1),1) ; repmat(name2,length(dat2),1)];
   grp = unique(group);
   grp_dat = [dat1,dat2];      
   hl = boxplot(dat,group);
    hold on
    for ih=1:size(hl,1)
        set(hl(ih,:),'Color','k');
    end
   for ii=1:length(dat1)
       for nbox = 1:length(grp)
           x(1,nbox) = nbox-0.05+rand*0.1;
       end
       y = grp_dat(ii,:);           

       plot(x,y,'.-','color',rgb(ii,:),'MarkerSize',15)
   end
   legend(names(:,1),'Location','westoutside');
   graph_beautify(gca,15,2,2)
   title(tt);   
end