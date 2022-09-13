function [maxPAC , nestingFreq, directPAC,hfreq] = PACestimate_cluster(F, sRate, lfreq, bandNested,run_subject, varargin)
%% Note 20.7.2015:  The generation of the frequencies was not correct. It is now corrected!


%function [maxPAC , nestingFreq, nestedFreq, directPAC] = PACestimate(F, sRate, bandNesting, bandNested)
%function [maxPAC , nestingFreq, nestedFreq, directPAC] = PACestimate(F, sRate, bandNesting, bandNested, ImagingKernel)
% Estimate of directPAC measures in signal time series
% INPUTS:
% F: an array or vector containing signal time series. If an array, each row is a time series
% sRate: signal sampling rate (in Hz)
% bandNesting: candidate frequency band of phase driving oscillations e.g.,
% [0.5 48] Hz - Note that cycle minimal frequency bandNesting(1) needs to be
% at least 10 times smaller than signal length (duration)
% bandNested: candidate frequency band of nested oscillatiosn e.g., [48,
% 300] Hz
% Optional:
% ImagingKernel: an imaging kernel from a Brainstorm EEG or MEG imaging model
%
% Outputs:
% for each signal (row of F)
% maxPAC: Value of maximum PAC
% nestingFreq = optimal nesting frequency (frequency for phase)
% nestedFreq = optimal nested frequency (frequency for amplitude)
% directPAC = full array of direct PAC measures for all frequyency pairs

% Control minimal fnesting value to avoid artifacts due to too short signal
% durations
signalDuration = size(F,2) / sRate;
if signalDuration < 10/lfreq(1)
    disp(sprintf('Warning: signal duration too short for requested minimal frequency for phase of %3.2f Hz: adjusting to %3.2f Hz',...
        lfreq(1), 10/signalDuration))
    lfreq(1) = 10/signalDuration;
end


allsources=size(F,1);

% Define bounds of nested high-frequency oscillations of interest
minGamma = bandNested(1);
maxGamma = bandNested(end); % Hz

% Initialize variables
dl = size(F,2);
maxPAC = zeros(1,allsources);
nestingFreq = zeros(1,allsources);
PAC_low = zeros(1,allsources);
PAC_frequency_low = zeros(1,allsources);
sp = get_signal_parameters('sampling_rate',sRate,...
    'number_points_time_domain',dl);
lfreq_vector = [0.5:0.5:12 14:2:48];
lfreq=lfreq_vector(lfreq_vector>=lfreq(1) & lfreq_vector<=lfreq(2));
chirpCenterFreqs = ...
make_center_frequencies(1,250,70,0.75); % min_freq, max_freq, numfreqs, min_freq_step chirpCenterFreqs =chirpCenterFreqs(chirpCenterFreqs >= minGamma & chirpCenterFreqs <= maxGamma);
hfreq = find( chirpCenterFreqs >= minGamma & chirpCenterFreqs <= maxGamma);% indices of center frequencies in the upper frequency range
hfreq=chirpCenterFreqs(hfreq);

nfHigh = length(hfreq); % number of cf bins to evaluate for PAC with lower-frequency oscillations
nfLow = length( lfreq );
chirpCenterFreqs=[lfreq hfreq];
s = make_signal_structure(...
    'raw_signal',F(1,:),...
    'output_type','analytic',...
    'signal_parameters',sp);

%matlabpool

chirpF = zeros(1, length(s.frequency_domain), length(chirpCenterFreqs));

% Make set of chirplets
for iif = 1:length(chirpCenterFreqs)
    f_g=[];
    f_g.fractional_bandwidth = 0.15;
f_g.chirp_rate = 0;

    f_g.center_frequency = chirpCenterFreqs(iif); % Hz
    chirp_f(iif) = make_chirplet(...
        'chirplet_structure',f_g,...
        'signal_parameters',sp);
    chirpF(1 , chirp_f(iif).signal_frequency_support_indices , iif) = ...
        chirp_f(iif).filter;
    
end
clear chirp_f

% Will contain all scores of PAC between all pairs of low-f and high-f bins
directPAC = zeros(allsources, nfLow, nfHigh);

%matlabpool

% Loop on all signal traces (rows of F)
for source = 1:allsources
 % Transform sensor time series into analytic signals
    s.frequency_domain = fft(F(source,:),...
        sp.number_points_frequency_domain,2);
    
    s.frequency_domain(:,sp.frequency_support<0)=0;
    s.frequency_domain(:,sp.frequency_support>0)=...
        2*s.frequency_domain(:,sp.frequency_support>0);
    
    
   % Define minimal frequency support with non-zeros chirplet coefficients
[~,scol] = find(s.frequency_domain ~= 0);
scol = max(scol)+1;
[chirprow] = find((chirpF(1,:,:)) ~= 0);
chirprow = max(chirprow)+1;
nfcomponents = min(chirprow,scol); % Minimal number of frequency coefficients
clear chirprow scol

BLKchirpF = chirpF( :,1:nfcomponents, :);

     
    src_fs = s.frequency_domain(1:nfcomponents);
    fs = ifft( src_fs(:, : , ones(1,length(chirpCenterFreqs) ) ) .* BLKchirpF,...
        sp.number_points_frequency_domain, 2);
    clear src_fs
    
    AMP = abs( fs(: , 1:sp.number_points_time_domain, nfLow+1:end ) );
    PHASE = exp (1i * angle( fs(: , 1:sp.number_points_time_domain , 1:nfLow ) ) );
    clear fs
   % Compute direct PAC index for each high-freq and low-freq pair
    for ihf = 1:nfHigh
        
        tmp2 =  AMP(: , 1:sp.number_points_time_domain , ihf);
        directPAC(source,:,ihf) = squeeze( sum( PHASE .* ...
            tmp2(:,:,ones(1,nfLow)), 2 ) );
        
    end
    clear tmp2
    
    % Finalize scaling of direct PAC metric
    tmp2 = sqrt( ( sum(AMP.*AMP, 2) ) );
    directPAC(source,:,:) = abs( directPAC(source,:,:) ) ...
        ./ tmp2(:,ones(1,size(directPAC,2)),:);
    directPAC(source,:,:) = directPAC(source,:,:)/sqrt(sp.number_points_time_domain); 
     clear AMP PHASE
    
end
    save([run_subject '.mat'], 'directPAC')


end


