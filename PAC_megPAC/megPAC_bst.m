function megPAC_bst( F, sRate, run_subject, lowfreqency, GoodChannel,OutputPath )
% This function generates the megPAC signal: a downsampled version of the
% input signal - the troughs and peaks of the given lowfrequency are
% determined for each input time series. The gamma amplitude from 80-150Hz
% is then calculated at the peaks and troughs of the lowfrequencz signal.
% The whole signal is in the end resampled to 10 Hz and saved

% Input: F(sources x time points)- source time series or sensor time series. This function does not
% transform from sensor to source space! In case sesor data are used as
% input the saving might not work in the brainstorm data base.
% sRate - sampling rate of F
% run_subject - string for the naming of the finally saved data
% lowfrequency -  vector of low frequencies in Hz with the size of the sources of
% F, which will be used to determine the peaks and troughs of each
% source/sensor time series. Ideally this is the nesting frequency from a
% Phase amplitude coupling for this particular time series for the gamma frequencies from 80-150 Hz.
% This is hard coded!!!!
% GoodChannel- Channel indices for the sensor data, which can be used - it
% is a needed (?) field in the results file for brainstorm
% OutputPath- Path, where the results should be saved. Ideally in the
% brainstorm database at the appropriate location.
% Caution: Depending on the data lengh, one cannot process all sources at
% the same time!

%   Author: Esther Florin

dl2=size(F,2);
%% Determine theta trough and gamma amplitude at trough

sp=get_signal_parameters('sampling_rate',sRate,...
    'number_points_time_domain',dl2);

% filter for low frequency phase extraction:
cf_list=make_center_frequencies(1,250,70,0.75); % min_freq, max_freq, numfreqs, min_freq_step
temp=(resample(F(1,:),1,round(sRate/10)));
sRate2 = sRate/round(sRate/10); %new sampling rate

HRF_conv=zeros(size(F,1),size(temp,2));
clear temp
for source=1:size(F,1) %Loop over all sources in F
    clear lf_g
    
    lf_g.center_frequency=lowfreqency(source); %  Hz
    lf_g.fractional_bandwidth=0.15;
    lf_g.chirp_rate=0;
    lf_g=make_chirplet(...
        'chirplet_structure',lf_g,...
        'signal_parameters',sp);  
    Data_sources=F(source,:);
    % Compute the maxima of the lowfrequency
    % precompute fft of raw_signal for later filtering loop:
    s=make_signal_structure(...
        'raw_signal',-Data_sources,...
        'output_type','analytic',...
        'signal_parameters',sp);
    % filter and extract triggering phase indicies for later averaging:
    fs=filter_with_chirplet(...
        'signal_structure',s,...
        'signal_parameters',sp,...
        'chirplet',lf_g);
    %Detection of phase maxima/minima of theta filtered signal
    phases=angle(fs.time_domain);
    [~,trigger_inds]=find_maxima(phases);
    
    
    % Compute the minima of the lowfrequency
    % precompute fft of raw_signal for later filtering loop:
    s=make_signal_structure(...
        'raw_signal',Data_sources,...
        'output_type','analytic',...
        'signal_parameters',sp);
    % filter and extract triggering phase indicies for later averaging:
    fs=filter_with_chirplet(...
        'signal_structure',s,...
        'signal_parameters',sp,...
        'chirplet',lf_g);
    %Detection of phase maxima/minima of lf_g filtered signal
    phases=angle(fs.time_domain);
    [~,trigger_inds2]=find_maxima(phases);
    
    trigger_inds=sort([trigger_inds trigger_inds2]);
    clear trigger_inds2
    % Calculate the gamma amplitude from 80-150 Hz
    ampmat=zeros(12,sp.number_points_time_domain);
    for f=51:62 %taking only frequencies from 80 to 150 Hz
        clear hf_g
        hf_g.center_frequency=cf_list(f); % Hz
        hf_g.fractional_bandwidth=0.15;
        hf_g.chirp_rate=0;
        hf_g=make_chirplet(...
            'chirplet_structure',hf_g,...
            'signal_parameters',sp);
        
        fs=filter_with_chirplet(...  % this type of filtering does not preserve amplitude units, only relative size within time series
            'signal_structure',s,...
            'signal_parameters',sp,...
            'chirplet',hf_g);
        temp=(abs(fs.time_domain));
        ampmat(f-50,:)=(temp-mean(temp))/std(temp); % zscore for normalization
    end
    
    
    ampmat=(mean(ampmat,1));
    temp=interp1(trigger_inds, ampmat(trigger_inds), 1:dl2);
    HRF_conv(source,:)=process_resample('Compute', temp, (0:dl2-1)/sRate, 10);
    
end

SaveResults(HRF_conv, [run_subject 'megPAC_lf_2-48_trough_gammaamp_resamples_bst'], 'lf 2-48 trough gamma amp bst' ,sRate2,GoodChannel,OutputPath)





%% ===== SAVE RESULTS =====
    function SaveResults(value, fileTag,Comment,sRate,GoodChannel,OutputPath)
        % Fill output structure
        kernelMat.ImageGridAmp = value;
        kernelMat.Comment      = Comment;
        kernelMat.sRate      = sRate;
        kernelMat.ImageGridTime         = 0:1/sRate:size(value,2)/sRate-1;
        kernelMat.DataFile=[];
        % kernelMat.ImagingKernel=[];
        % kernelMat.Whitener=[];
        kernelMat.Time         = 1:size(value,2);
        kernelMat.GoodChannel=GoodChannel;
        c = clock;
        strTime = sprintf('_%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
        % Output filename
        OPTIONS.ResultFile = fullfile(OutputPath, ['results_' fileTag, strTime, '.mat']);
        % Save file
        save(OPTIONS.ResultFile, '-struct', 'kernelMat');
        % Output filename
    end

end
