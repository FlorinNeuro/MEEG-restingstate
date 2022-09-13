% This script calculates the resting state networks with ICA and performs a
% bootstrap across the subjects. For each bootstrap itteration the
% resulting networks are saved in the specified brainstorm database

%Esther Florin 2020
clear


pathStr = '/gpfs/project/florine/icatb/';
addpath(genpath(pathStr));

brainstorm_db='/gpfs/project/florine/resting_young/';
Condition='wide_band_fir'; %Select for which data the calculation is done.

%%

test=dir([ brainstorm_db '\data\']);
subject={};

subject={'S001_MEG'; 'S002_MEG';'S003_MEG'; 'S004_MEG'; 'S006_MEG'; 'S007_MEG'; 'S008_MEG'; 'S009_MEG'; 'S010_MEG'; 'S011_MEG'; 'S012_MEG'; 'S013_MEG'; 'S014_MEG'; ...
    'S015_MEG'; 'S017_MEG'; 'S019_MEG'; 'S020_MEG'; 'S021_MEG'; 'S022_MEG'; 'S023_MEG'; 'S024_MEG'; 'S025_MEG'; 'S026_MEG';};
nbsubject=length(subject);
Frequenz={'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma'};
for rep=1:100
    
for Frequenzitter=1:length(Frequenz);
    diffTimePoints=[];
    Dat_concat=[];
%     isubjects=[];
subjectindex=randi(nbsubject,nbsubject,1);
    for subjects=1:length(subject)
        TMP=load([brainstorm_db 'data/Group_analysis/Envelope_' Condition '/results_' subject{subjectindex(subjects)} '_' Frequenz{Frequenzitter} '_Envelope_' subject{subjectindex(subjects)} '_zscore' ]);
        
        TMP= bsxfun(@rdivide,...
            bsxfun(@minus,TMP.ImageGridAmp,mean(TMP.ImageGridAmp,1)),std(TMP.ImageGridAmp,[],1));
        
        Dat_concat=[Dat_concat TMP];
        diffTimePoints=[diffTimePoints size(TMP, 2)];
%         isubjects=[isubjects subjects];
    end
    
    % %% Extracting ICAs according to Brookes et al. 2011
    %
    % % 1. pca with 30 components
    % % 2. prewhitening the data
    % % 3. temporal ICA with 25 components
    % % 4. Pearson correlation between each tICA and each vocel time series
        pca_opts.stack_data='yes';
    pca_opts.storage='full';
    pca_opts.precision='double';
    pca_opts.eig_solver='selective';
    
    [V, Lambda] = icatb_calculate_pca(Dat_concat, 30, 'type', 'standard', 'whiten', 0, 'verbose', 1, 'preproc_type', 'remove mean per timepoint', ...
        'pca_options', pca_opts, 'varToLoad', [], 'dims', size(Dat_concat));
    
    % data voxels x time (time for all subjects together Afterwards the mean
    % per time point is removed
    
    % Output
    %  V - Eigen vectors of dimensions columns x numComp (i.e. time points by
    %  components
    %   Lambda - Diagonal matrix of eigen values sorted in ascending order.
    %Group ICA spatial
    
    [icasig_tmp, A, W] = icatb_fastICA(V', 'epsilon', 1E-4,'maxNumIterations', 10000, 'maxFinetune', 5, 'sampleSize', 1, ...
        'numOfIC', 20, 'verbose', 'on', 'approach', 'defl', 'finetune', 'pow3', 'stabilization', 'off', 'g', 'pow3', 'only', 'all');
 
    icasig = icatb_remove_mean(icasig_tmp)';
   
    
    SpatialMap_group=corr( icasig, Dat_concat');

    kernelMat.ImageGridAmp =SpatialMap_group';
    kernelMat.Comment      = ['Groupresult_' Frequenz{Frequenzitter} '_Envelope_TC_cor_data' Condition ' ' num2str(rep)];
    kernelMat.sRate      = 1;
    kernelMat.ImageGridTime         =1:size(SpatialMap_group,1);
    kernelMat.DataFile=[];
    kernelMat.Time         = 1:size(SpatialMap_group,1);
    kernelMat.SurfaceFile=['@default_subject\tess_cortex_pial_low.mat'];
    kernelMat.HeadModelType='surface';
    
    % Output filename
    OPTIONS.ResultFile = fullfile([brainstorm_db '/data/Group_analysis/Group_ICA_Brookes_based_' Condition '_23Subjects_MEG/'], ['results_Group_' Frequenz{Frequenzitter} 'Brookes_GroupICA_backprojected_TC_corr_ICASSO_' Condition '_' num2str(rep)] );
    % Save file
    save(OPTIONS.ResultFile, '-struct', 'kernelMat');
    
    
end
end

