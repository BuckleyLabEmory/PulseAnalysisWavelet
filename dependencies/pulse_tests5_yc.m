function [tbl_qc,test_tbl_out] = pulse_tests5(tbl,p_i,ABP_or_Dbfit_wv)
% pulse_names = {p,p_ref} usually bfi and abp
    all_tests = {};
    all_result = {};
    wlabel = [];
    plabel = {};
    test_tbl = {};
%%
    tbl1 = tbl;
    winds = unique(tbl1.wnum);   
    tbl_qc = tbl1;
%%    

    for w=1:length(winds)
    tests = [];
    result = [];
    test_count = 0;
    % pulses are tested by window - if only and ABP pulse or BFI pulse fails they are both rejected to keep the dataset paired    

        wselect = tbl1.wnum==winds(w);
        tbl_rw = tbl1(wselect,:);
        
        tpulse = tbl_rw.tbin_ln{1}; % the same across pulses and windows
        pulse_names = {p_i,ABP_or_Dbfit_wv};
        % pulse time series
        pulses = {};
        pulses{1} = tbl_rw.avg_bfi_ln{1}; % sds2 pulse from window
        pulses{2} = tbl_rw.avg_abp_ln{1}; % abp pulse from window

        % pulse peak
        pulse_peak = {};
        pulse_peak{1} = tbl_rw.psf_ln{1};
        pulse_peak{2} = tbl_rw.psp_ln{1};
 
        % pulse amp
        pulse_height = {};
        pulse_height{1} = tbl_rw.fph_ln{1};
        pulse_height{2} = tbl_rw.pph_ln{1};

        % pulse end diastolic
        pulse_end = {};
        pulse_end{1} = tbl_rw.edf_ln{1};
        pulse_end{2} = tbl_rw.edp_ln{1};

        % pulse onset
        pulse_onset = {};
        pulse_onset{1} = tbl_rw.osf_ln{1};
        pulse_onset{2} = tbl_rw.osp_ln{1};

        % mean 
        pulse_mean = {}
        pulse_mean{1} = tbl_rw.mf_ln{1};
        pulse_mean{2} = tbl_rw.mp_ln{1};
        
        accept = ones(1,length(pulses));
        
        % this means for bfi and abp (or other reference)
        for p=1:length(pulses)
    
                pulse = pulses{p};
                ph = pulse_height{p};
                pp = pulse_peak{p};
                pe = pulse_end{p};
                po = pulse_onset{p};
                p_avg = pulse_mean{p};
                
                t_pp = tpulse(pulse == pp);
                t_pe = tpulse(pulse == pe);
                if find(pulse == pe) ~= length(pulse)
                    error('edf not last timepoint in pulse')
                end
                
                
                 %% test 1: are there spikes larger than the systolic peak later in the pulse?
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: spikes at t  > 0.6'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    
                 end_pulse = pulse(tpulse>0.6);
                if any(end_pulse > max(pulse))
                    accept(p) = 0;
                    fprintf('Result: True\n')
                    result(test_count) = 1;
                else 
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end
    
                %% test 2: is the pulse height negative?
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: negative amplitude'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
                if ph < 0 
                    accept(p) = 0;
                    fprintf('Result: True\n')
                    result(test_count) = 1;
                else 
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end
                
                %% test 3: If the slope of systolic period of pulse is negative
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: negative systolic slope - method 1'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
                pulse_start = pulse(tpulse<t_pp);
                if length(pulse_start) < 3
                   upstrk = pulse_start;
                else 
                    upstrk = pulse_start(1:3);
                end
                if mean(diff(upstrk)) < 0
                   accept(p) = 0;
                   fprintf('Result: True\n')
                   result(test_count) = 1;
                else 
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end
                
                %% test 4: If end diastolic flow greater than 2x onset
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: diastolic more than 2x onset'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
                
                if pe/po > 2
                    accept(p) = 0;
                    fprintf('Result: True\n')
                    result(test_count) = 1;
                else
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end
    
                %% test 6: is the edf greater than the mean of the pulse
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: edf greater than mean'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n']);
    
                if pe > p_avg
                    accept(p) = 0;
                    fprintf('Result: True\n')
                    result(test_count) = 1;
                else
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end            
              
                %% test 7: is the onset greater than the mean of the pulse
                test_count = test_count + 1;
                tests{test_count} = [pulse_names{p} ' pulse: onset greater than mean'];
                fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n']);
    
                if po > p_avg
                    accept(p) = 0;
                    fprintf('Result: True\n')
                    result(test_count) = 1;
                else
                    fprintf('Result: False\n')
                    result(test_count) = 0;
                end                   

                
                %% test 8: If reference pulse is abp and it is under 10mmHg
                if strcmp(ABP_or_Dbfit_wv, 'ABP')    
                    test_count = test_count + 1;
                    tests{test_count} = [pulse_names{p} ' pulse: ref pulse is abp and under 10 mmHg'];
                    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n']);
                    
                    if p==2 & all(pulse<10)
                        accept(p) = 0;
                        fprintf('Result: True\n')
                        result(test_count) = 1
                    else
                        fprintf('Result: False\n')
                        result(test_count) = 0;    
                    end
                
                         
                end

            
        end
        test_tbl = [test_tbl; [num2cell(repmat(winds(w),length(result),1)),tests',num2cell(result)']];
        test_tbl_out = cell2table(test_tbl);
        test_tbl_out.Properties.VariableNames = {'wnum','test','result'};
        if any(accept == 0)
            tbl_qc.passed_pulse_tests(wselect) = 0;
        elseif all(accept == 1)
            tbl_qc.passed_pulse_tests(wselect) = 1;
        end 
    end
end