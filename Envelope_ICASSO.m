% % This script calculates independent compoents with ICASSO. Can then be used to obtain the resting state networks
% This can be done along the lines of the Envelope_ICA script

clear all

pathStr = 'F:\Esther\Matlab_scripts_toolboxes\FastICA_25\'; %path to the canolty code
addpath(genpath(pathStr));

pathStr = 'F:\Esther\Matlab_scripts_toolboxes\icatb\';
addpath(genpath(pathStr));

brainstorm_db='F:\Esther\RS_MEEG_fMRI\brainstorm_db\resting_young\';
Condition='wide_band_iir';

subject={'S001_MEG'; 'S002_MEG';'S003_MEG'; 'S004_MEG'; 'S006_MEG'; 'S007_MEG'; 'S008_MEG'; 'S009_MEG'; 'S010_MEG'; 'S011_MEG'; 'S012_MEG'; 'S013_MEG'; 'S014_MEG'; ...
    'S015_MEG'; 'S017_MEG'; 'S019_MEG'; 'S020_MEG'; 'S021_MEG'; 'S022_MEG'; 'S023_MEG'; 'S024_MEG'; 'S025_MEG'; 'S026_MEG';};

Frequency={'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma'};
for rep=1:4
for Frequenzitter=1:(length(Frequency))
    diffTimePoints=[];
    Dat_concat=[];
    isubjects=[];
    for subjects=1:length(subject)
        TMP=load([brainstorm_db 'data\Group_analysis\Envelope_' Condition '\results_' subject{subjects} '_' Frequency{Frequenzitter} '_Envelope_' subject{subjects} '_ssmooth_zscore' ]);
        
        TMP= bsxfun(@rdivide,...
            bsxfun(@minus,TMP.ImageGridAmp,mean(TMP.ImageGridAmp,1)),std(TMP.ImageGridAmp,[],1));
        
        Dat_concat=[Dat_concat TMP];
        diffTimePoints=[diffTimePoints size(TMP, 2)];
        isubjects=[isubjects subjects];
    end
    
    % %% Extracting ICAs according to Brookes et al. 2011
    %
    % % 1. pca with 30 components
    % % 2. prewhitening the data
    % % 3. temporal ICA with 25 components
    % % 4. Pearson correlation between each tICA and each voxel time series
    
    [iq, A, W, S, sR]=icasso(Dat_concat, 500,   'approach', 'defl', 'numOfIC', 20, 'epsilon', 1E-4,'maxNumIterations', 100000, 'maxFinetune', 5, ...
        'verbose', 'on',  'finetune', 'pow3','g', 'pow3', 'stabilization', 'off', 'lastEig',30, 'vis','off');
    
    [B,I]=sort(iq, 'descend');
    save(['ICASSO_run_' num2str(rep) '_'  Frequency{Frequenzitter} '_' Condition '500_repetitions'], 'iq', 'A', 'W', 'S', 'sR', '-v7.3');

   
    
end
end

