function HilbertEnvelope_smallband_iir(Data,sRate,frequency_bins, run_subject, OutputPath)
%Aim reproduce Brookes et al. 2011 PNAS
%Determine envelope in frequency_bins from Data and resample to 1Hz

%Update 8.11.2016: zscore is put back in for the filtered data, because otherwise
% there would be no correction for the 1/f decrease
% IN addition the frequency band width was corrected. Before it read frequency_bins(ff,2)<30 && frequency_bins(ff,2)>4Brookes approach


% Esther Florin 22.10.2016
dl=length(Data);
for ff=1:length(frequency_bins)
    if frequency_bins(ff,2)<=30 && frequency_bins(ff,2)>4
        freqbands=frequency_bins(ff,1):2:frequency_bins(ff,2);
    elseif frequency_bins(ff,2)>=30
        freqbands=frequency_bins(ff,1):5:frequency_bins(ff,2);
    elseif frequency_bins(ff,2)<=4
        freqbands=frequency_bins(ff,1):1:frequency_bins(ff,2);
        
    end
    temp=zeros(length(freqbands)-1,size(Data,1),size(Data,2));
    for ff_2=1:length(freqbands)-1
        
        temp(ff_2,:,:)=bst_bandpass_filtfilt(Data, sRate,freqbands(ff_2),freqbands(ff_2+1),0,'iir'); %%iir filter ensures that low-pass has not a too large passband
                temp(ff_2,:,:)= bsxfun(@rdivide,...
            bsxfun(@minus, squeeze(temp(ff_2,:,:)),mean(squeeze(temp(ff_2,:,:)),2)),std(squeeze(temp(ff_2,:,:)),[],2));%/std(temp(ff_2,:,:));

        temp(ff_2,:,:)=(abs(hilbert(squeeze(temp(ff_2,:,:))')))';
    end
    temp=squeeze(mean(temp,1));
    [Data_f(ff,:,:),timeout]=process_resample('Compute',temp,(0:dl-1)/sRate,1,'resample-cascade');
    clear temp
end
save([OutputPath 'HilbertEnvelope_small_band_iir' run_subject ],'Data_f' ,'timeout','frequency_bins','-v7.3')

