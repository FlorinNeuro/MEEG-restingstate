function HilbertEnvelope_wideBand_iir(Data,sRate,frequency_bins, run_subject, OutputPath)
%Aim reproduce Brookes et al. 2011 PNAS
%Determine envelope across whole bands (delta, theta, alpha,...) without a
%correction for 1/f and with the hilbert transform across the whole band!
%This seems to be the way Brookes et al. performed the analysis
%zscore is taken out of function. Makes it closer to
%Brookes approach

% Esther Florin 03.05.2018
dl=length(Data);

for ff=1:length(frequency_bins)
        temp=bst_bandpass_filtfilt(Data, sRate,frequency_bins(ff,1),frequency_bins(ff,2),0,'iir'); %%iir filter ensures that low-pass has not a too large passband
       
        temp=(abs(hilbert((temp)')))';
    [Data_f(ff,:,:),timeout]=process_resample('Compute',temp,(0:dl-1)/sRate,1,'resample-cascade');
    clear temp
end
save([OutputPath 'HilbertEnvelope_wideBand_iir' run_subject ],'Data_f' ,'timeout','frequency_bins','-v7.3')

