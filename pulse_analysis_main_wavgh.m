%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Copyright (c) 2024 Srinidhi Bharadwaj & Tara Urner              
%  http://buckleylab.gatech.edu/           
%                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use ABP cardiac pulses and wavelet transformed BFI data (BFI_wv) to
% extract overlapping BFI data 
% 
% Inputs:
% 
% data: Folder of input data in the expected format.
% savedir: directory to save output data 

PATH_TO_DATA = 'ExampleData/'; % HC
numnans=2;
subject = '27_01_FC27';
savedir= ['ExampleOutput/' subject filesep];
if ~exist(savedir,"dir")
    mkdir(savedir)
end

load(strcat(PATH_TO_DATA,subject,'--Dbfit.mat'))
load(strcat(PATH_TO_DATA,subject,'--TimeAxis_DCS.mat'))
load(strcat(PATH_TO_DATA,subject,'--raw_ABP_masterdata.mat'))
load(strcat(PATH_TO_DATA,subject,'--Marks_index.mat'))

% BFI and ABP load
dt_ABP = pl_analysis_t_compare(TimeAxis_DCS,raw_ABP_masterdata,Dbfit,numnans);  
BFI = dt_ABP.BFI_interp; 
ABP = dt_ABP.ABP_binterp;
t_BFI = dt_ABP.t_BFI_interp;

% BFI and BFI_wv load
 
disp("Starting wavelet")
MRA_dim = 4;
init_BFI = Dbfit'; 
mra_BFI = modwtmra(modwt(init_BFI,8)); 
mra_BFI_wvwv = mra_BFI(MRA_dim,:) + mra_BFI(MRA_dim-1,:) + mra_BFI(MRA_dim-2,:); % change number of wavelets here

BFIwv_120Hz = interp1(1:numel(mra_BFI_wvwv),mra_BFI_wvwv, 1:1/6:numel(mra_BFI_wvwv))';
t_120Hz = interp1(1:numel(TimeAxis_DCS),TimeAxis_DCS, 1:1/6:numel(TimeAxis_DCS))'; 

t_120Hz(isnan(BFIwv_120Hz)) = [];
BFIwv_120Hz(isnan(BFIwv_120Hz)) = []; 


Dbfit_wv = [BFIwv_120Hz, t_120Hz];
dt_wv = pl_analysis_t_compare_wv(TimeAxis_DCS,Dbfit_wv,Dbfit,numnans);

dt_all = {dt_ABP,dt_wv};
state = 'Not Defined'; % hard coding state but remove if using marks
filenm = strcat(subject,'-wvwvwv-'); %uncomment for rs
type_ref = {'abp','bfi_wv'};

%%%%%% PIPELINE HERE IS SIMPLIFIED PULSE_ANALYSIS_PIPELINE https://github.com/BuckleyLabEmory/PulseAnalysis/

for k = 1:2 % k=1 analysis according to ABP waveform, k=2 analysis according to BFI_wv waveform
    % init
    pdo = struct();
    ftbl = table();
    wtbl = table();
    ptbl = table();
    dat_tbl = table();
    wv = table();
    rc = 1;

    % Full time series table
    dt = dt_all{k}; % if k=2, then all instances of ABP below are actually bfi_wv
    ftbl_n =pl_analysis_ftbl_gen(subject,state,dt.t_ABP_corr,dt.ABP_corr,dt.f_ABP,dt.ABP_downsample,...
        dt.t_DCS_corr,dt.DCS_corr,dt.f_DCS,dt.t_BFI_interp,dt.BFI_interp,dt.ABP_binterp);
    
    ftbl_n.dt = dt; % sampling information for raw ABP and BFI data
   
    ftbl_n.corr_method = {'xcorr'};

    sds =2;
    sds = repmat({sds},height(ftbl_n),1);
    tmp = cell2table(sds);  
    ftbl_n = [tmp, ftbl_n];

    ftbl = [ftbl ; ftbl_n];
    
    % Create windows
    window_size=15;
    %%% WINDOWING %%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Window entire signal with no overlap
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        
   [window_start, window_end] = window_abp_dcs(ftbl_n.t_bfi{:},ftbl_n.fs_bfi,window_size);    
    
    
    % Sort windows into states - either defined by marks or 'none'.
    wtbl_n = window_state_sort(ftbl_n,window_start,window_end);

    %%%%%%%%%%%%%%% Most generalized pulse analysis %%%%%%%%%%%%%%%%%%%%%%%
    % For anyone wanting to insert this pulse analysis method 
    % into an existing pipeline that handles all of the above 
    % loading/syncing/windowing for a different experimental setup, 
    % these three subfunctions are what you want to use.

    % Extract the window data, identify ABP pulses, calculate shift between ABP and BFI
    wtbl_n = pl_analysis_interp_wtbl_gen_v2(ftbl_n,wtbl_n,numnans);

    % Extract BFI pulses based on ABP reference, add to pulse table
   
    pl_n = pl_analysis_ptbl_gen_v3(ftbl_n,wtbl_n);
    ptbl = [ptbl ; pl_n];
   
    % Generate window-averaged pulses and analyze pulse morphology
    % If you want to change bin size that is hard-coded in here
    if k == 1
        [dat_tbl_n,wtbl_n] = dat_tbl_gen_v3_yc(wtbl_n,pl_n,"ABP");
    else
        [dat_tbl_n,wtbl_n] = dat_tbl_gen_v3_yc(wtbl_n,pl_n,"Dbfit_wv");
    end 
    wtbl = [wtbl ; wtbl_n];
    dat_tbl = [dat_tbl ; dat_tbl_n];

    % Save pulses
    states = string(dat_tbl.Properties.VariableNames(2:end));
    dat_tbl_qc = dat_tbl(:,'name');
    for s=1:length(states)
        dat_tbl_qc.(states{s}) = cell(height(dat_tbl),1);
        for subs=1:height(dat_tbl)
            % %%%%% Use quality control filters
            qc = dat_tbl{subs,states{s}}{1}.passed_pulse_tests == 1;
            dat_tbl_qc(subs,states{s}) = {dat_tbl{subs,states{s}}{1}(qc,:)};
        end
        if k==1
            state_pulse_plot_abp(dat_tbl_qc,states{s},savedir,filenm)
            save(strcat(savedir,filenm,'abp_dat_tbl_qc.mat'),"dat_tbl_qc");
        else
            state_pulse_plot_bfi(dat_tbl_qc,states{s},savedir,filenm)
            save(strcat(savedir,filenm,'bfiwv_dat_tbl_qc.mat'),"dat_tbl_qc");
        end
    end
end 

clear Dbfit TimeAxis_DCS raw_ABP_masterdata Dbfit_wv

