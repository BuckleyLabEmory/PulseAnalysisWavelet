function [load_table,cvr_db] = parse_cvr_db(db_dir)
% Tara Urner, Buckley Lab 2020
% load_table: the table of files as loaded
% patient_db: a table parsed on patient id and session id
% structure:
% patient_db{:,1} = patient id
% patient_db{:,2} = locations of patient data in load_table
% patient_db{:,3} = sessions data
% patient_db{:,3}{:,1} = session id
% patient_db{:,3}{:,2} = location in load table
% patient_db{:,3}{:,3} = filepaths
% patient_db{:,3}{:,3}.data_type = type of data contained in .csv
% patient_db{:,3}{:,3}.folder = folder location of data
% patient_db{:,3}{:,3}.name = name of file
% patient_db{:,3}{:,3}.path = full path for loading.
    cd(db_dir)
    pdb_dir = dir('*mat');
    pdb_t = struct2table(pdb_dir);
    fs = pdb_t.name;
    a = string(fs);
    %%
    % Find all patient/session names
    split_names = arrayfun( @(x) strsplit( x, '--' ), a, 'UniformOutput', false );
    try % original cvr database
        subject_idx = cellfun(@(x) strsplit(x,'-'),cellfun(@(x) x{1},split_names,'UniformOutput',false),'UniformOutput',false);
        pdb_t.sname = cellfun(@(x) x{1},split_names,'UniformOutput',false);
        pdb_t.subject_idx = cellfun(@(x) [x{length(x)-1} '-' x{length(x)}],subject_idx,'UniformOutput',false);
        pdb_t.session_idx = cellfun(@(x) [x{1} '_' x{2} '_' x{3}],subject_idx,'UniformOutput',false) ;
    catch % updated
        try
        subject_idx = cellfun(@(x) strsplit(x,'_'),cellfun(@(x) x{1},split_names,'UniformOutput',false),'UniformOutput',false);
        pdb_t.sname = cellfun(@(x) x{1},split_names,'UniformOutput',false);
        pdb_t.subject_idx = cellfun(@(x) [x{1}],subject_idx,'UniformOutput',false);
        pdb_t.session_idx = cellfun(@(x) [x{2} '_' x{3}],subject_idx,'UniformOutput',false) ; 
        catch % updated
        subject_idx = cellfun(@(x) x{1},split_names,'UniformOutput',false);
        pdb_t.sname = cellfun(@(x) x{1},split_names,'UniformOutput',false);
        pdb_t.subject_idx = subject_idx;
        pdb_t.session_idx =  subject_idx;  
    
        end
    end
    data_type = cellfun(@(x) strsplit(x,'.'),cellfun(@(x) x{2},split_names,'UniformOutput',false),'UniformOutput',false);
    pdb_t.data_type = cellfun(@(x) x{1},data_type,'UniformOutput',false);
    unique_sidx = unique(pdb_t.subject_idx);
    for i=1:length(unique_sidx)
        unique_sidx{i,2} = cell2mat(cellfun(@(x) any(strfind(x,unique_sidx{i})),fs,'UniformOutput',false));
        sessions = table2cell(unique(pdb_t(unique_sidx{i,2},'session_idx')));
        %session_id = pdb_t.session_idx;
        for j=1:size(sessions,1)
            % Unique sessions
            sessions{j,2} = cell2mat(cellfun(@(x) any(strfind(x,sessions{j})),pdb_t.session_idx,'UniformOutput',false));
            session_info = string(table2cell(pdb_t(sessions{j,2},'data_type')));
            session_info(:,2) = string(table2cell(pdb_t(sessions{j,2},'folder')));
            session_info(:,3) = string(table2cell(pdb_t(sessions{j,2},'name')));
            session_info(:,4) = strcat(session_info(:,2),filesep,session_info(:,3));
            session_info = array2table(session_info,'VariableNames',{'data_type','folder','name','path'});
            sessions{j,3} = session_info;
        end
        unique_sidx{i,3} = sessions;
        load_table = pdb_t;
        cvr_db = unique_sidx;
    end

end