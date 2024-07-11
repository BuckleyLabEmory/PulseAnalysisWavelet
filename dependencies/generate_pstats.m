function pstats = generate_pstats(sub_stats)
%%
% Taken from N2021_11_8.mlx
state = categorical(sub_stats.state);
param = categorical(sub_stats.param_name);
sdss = categorical(sub_stats.sds);
params = categories(param);
states = categories(state);
sds = categories(sdss);

%varnames = {'param','param_name','data','names'};
%varnames = {};
for st = 1:length(states)
    for sd = 1:length(sds) 
        %vars{st,sd} = [states{st} '_' sds{sd}];
        vst{st,sd} = states{st};
        vsd{st,sd} = sds{sd};
        v = [states{st} '_' sds{sd}];
        v1 = strsplit(v,'_');
        vars{st,sd} = [v1{:}];
    %cellfun(@(x) strfind(x,'_'),v)
    end
end



if any(size(vars)~=[2,2])
    error('Wrong size for these stats')
else
    vars_c{1} = [vars{1,1} '_' vars{1,2}];
    vars_c{2} = [vars{1,1} '_' vars{2,1}];
    vars_c{3} = [vars{1,2} '_' vars{2,2}];
    vars_c{4} = [vars{2,1} '_' vars{2,2}];
end


varnames = [{'param','param_name','data'} reshape(vars,1,4) vars_c {'ids'}];


pstats = cell2table(cell(length(params),length(varnames)),'VariableNames',varnames);
% 
for p=1:length(params)
    ptbl = sub_stats(param==params{p},:);
    pstats.data{p} = ptbl;
    pstats(p,{'param','param_name'}) = ptbl(1,{'param','param_name'}); % all rows of ptbl should match
    for st = 1:length(states)
        for sd = 1:length(sds) 
            tvals = sub_stats(state == vst{st,sd} & sdss == vsd{st,sd} & param==params{p},:);
            pstats(p,vars{st,sd}) = {tvals.mean};
            
            id1 = cellfun(@(x) strsplit(x,' '),tvals.name,'UniformOutput',false);
            %ids = cellfun(@(x) x{3}, id1,'UniformOutput',false);
            %%
            ids = {};
            for ii=1:size(id1,1)
                t1 = any(cell2mat(strfind(id1{ii},'FC')));
                if t1
                    ids{ii} = str2num(id1{ii}{3}(3:4));
                else
                    ids{ii} = str2num(id1{ii}{1});
                end
            end
            %%
            pstats.ids(p) = {ids};
        end
    end

    for c=1:4
        comp = strsplit(vars_c{c},'_');
        cc1 = pstats{p,comp{1}};
        cc2 = pstats{p,comp{2}};
        pstats{p,vars_c{c}} = {signrank(cc1{:},cc2{:})};
    end
end
end