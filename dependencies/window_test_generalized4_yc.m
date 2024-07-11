function [tests,result,window_BFI_update,window_ABP_update,window_pulsatility,window_keep,discard_reason] = window_test_generalized4(window_BFI,window_ABP,window_time,numnans, ABP_or_Dbfit_wv)
tests = [];
result = [];
test_count = 0;
window_keep = 1;
window_pulsatility = NaN;
discard_reason = 'none';
window_BFI_update = window_BFI;
window_ABP_update = window_ABP;

%% Window is all NaNs 
    test_count = test_count + 1;
    tests{test_count} = 'Window: All NaNs';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    if all(isnan(window_BFI)) || all(isnan(window_ABP))
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'All NaNs';
        return
    end
    
%% NaN threshold
    %%%%% ABP %%%%%%
    test_count = test_count + 1;
    tests{test_count} = 'Window: NaN threshold - ABP';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    % Pad the window with 0 on either side to prevent the edge case of a
    % window with a big chunk of NaNs at the end from getting through here
    
    [aB,aN,aIB] = RunLength(window_ABP);
    t = isnan(aB(aN > numnans));
    if any(t)
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'NaNs Threshold - ABP';
        return
    else
    result(test_count) = 0;
    fprintf('Result: False\n')
    end
    
    test_count = test_count + 1;
    tests{test_count} = 'Window: NaN threshold - BFI';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    %%%%% BFI %%%%%%
    [bB,bN,bIB] = RunLength(window_BFI);
    t = isnan(bB(bN > numnans));
    if any(t)
        result(test_count) = 1;
        fprintf('Result: True\n')
        window_keep = 0;
        discard_reason = 'NaNs Threshold - BFI';
        return
    else
        result(test_count) = 0;
        fprintf('Result: False\n')
    end

    %% ABP artifact (CNAP calibration, blood draw from A-line, etc)
    test_count = test_count + 1;
    tests{test_count} = 'Window: ABP Artifact';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])
    
    if strcmp(ABP_or_Dbfit_wv, 'ABP')    
       if length(find(window_ABP < 5)) > 10
           result(test_count) = 1;
           fprintf('Result: True - ABP Artifact...\n')
           window_keep = 0;
       else
           result(test_count) = 0;
           fprintf('Result: False \n')
       end
    end
%% Interpolate NaNs 
    
    %test_count = test_count + 1;
    %tests{test_count} = 'Window: Interpolate NaNs - BFI';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])   
    
    %%%%% BFI %%%%%%
    if any(isnan(window_BFI))
        %result(test_count) = 1;
        fprintf('Result: True - BFI. Interpolating...\n')
        window_BFI_update = interp_nans(window_BFI,window_time,0);
    else
        %result(test_count) = 0;
        fprintf('Result: False \n')
    end

    %test_count = test_count + 1;
    %tests{test_count} = 'Interpolate NaNs - ABP';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])

    %%%%% ABP %%%%%%
    if any(isnan(window_ABP))
        %result(test_count) = 1;
        fprintf('Result: True - ABP. Interpolating...\n')
        window_ABP_update = interp_nans(window_ABP,window_time,0);
    else
        %result(test_count) = 0;
        fprintf('Result: False \n')
    end
    
    %% Cardiac pulsatility in BFI
    test_count = test_count + 1;
    tests{test_count} = 'Window: BFI pulsatility';
    fprintf(['Test #' num2str(test_count) ': ' tests{test_count} '\n'])

    F0 = mean(window_BFI_update);%(1:100));
    fs = mode(1./(diff(window_time)));
    try
    [AF,~] = pwelch(window_BFI_update,[],[],[],fs);
    [AP,f] = pwelch(window_ABP_update,[],[],[],fs);
    catch
        % Special handling for NaN at start or end of BFI window
        use_idx = 1:length(window_BFI_update);
        if isnan(window_BFI_update(1)) | isnan(window_ABP_update(1))
            use_idx = use_idx(2:end);
        elseif isnan(window_BFI_update(end)) | isnan(window_ABP_update(end))
            use_idx = use_idx(1:end-1);
        end
        window_BFI_update = window_BFI_update(use_idx);
        window_ABP_update = window_ABP_update(use_idx)
        [AF,~] = pwelch(window_BFI_update,[],[],[],fs);
        [AP,f] = pwelch(window_ABP_update,[],[],[],fs);
    end
    AP = sqrt(AP(:)); %To get amplitude
    AF = sqrt(AF(:));
    I = (f > 0.5) & (f < 2);
    tmp = f(I);
    [~,loc] = max(AP(I));
    hr = tmp(loc);
    Ihr = find(f >= hr,1);
    F1 = AF(Ihr);
    Fratio = F1/F0;
    window_pulsatility = Fratio;
    if Fratio < 0.2
        fprintf('Result: True...\n')
        result(test_count) = 1;
        window_keep = 0;
    else
        result(test_count) = 0;
        fprintf('Result: False \n')
    end
end