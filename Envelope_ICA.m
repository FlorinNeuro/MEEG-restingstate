% This script calculates the resting state networks with ICA and saves it
% into a brainstorm data base
%Esther Florin 2020

clear
rng('shuffle');


pathStr = 'C:\Users\Admin\Documents\icatb';
addpath(genpath(pathStr));

brainstorm_db='F:\Esther\RS_MEEG_fMRI\brainstorm_db\resting_young/';
Condition='wide_band_iir';


subject={'S001_MEG'; 'S002_MEG';'S003_MEG'; 'S004_MEG'; 'S006_MEG'; 'S007_MEG'; 'S008_MEG'; 'S009_MEG'; 'S010_MEG'; 'S011_MEG'; 'S012_MEG'; 'S013_MEG'; 'S014_MEG'; ...
    'S015_MEG'; 'S017_MEG'; 'S019_MEG'; 'S020_MEG'; 'S021_MEG'; 'S022_MEG'; 'S023_MEG'; 'S024_MEG'; 'S025_MEG'; 'S026_MEG';};
nbsubject=length(subject);
Frequency={'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma'};

    
for Frequenzitter=1:length(Frequency);
    diffTimePoints=[];
    Dat_concat=[];
    icasig_tmp=[];
    
subjectindex=1:nbsubject;
    for subjects=1:length(subject)
        TMP=load([brainstorm_db 'data/Group_analysis/Envelope_' Condition '/results_' subject{subjectindex(subjects)} '_' Frequency{Frequenzitter} '_Envelope_' subject{subjectindex(subjects)} '_ssmooth_zscore' ]);
        
        TMP= bsxfun(@rdivide,...
            bsxfun(@minus,TMP.ImageGridAmp,mean(TMP.ImageGridAmp,1)),std(TMP.ImageGridAmp,[],1));
        
        Dat_concat=[Dat_concat TMP];
        diffTimePoints=[diffTimePoints size(TMP, 2)];
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
        'numOfIC', 20, 'verbose', 'on', 'approach', 'defl', 'finetune', 'off', 'stabilization', 'off', 'g', 'pow3', 'only', 'all');
     if ~isempty(icasig_tmp)

    icasig = icatb_remove_mean(icasig_tmp)';
    count=0;
    for nSet=1:length(subject)
        TMP=load([brainstorm_db 'data/Group_analysis/Envelope_' Condition '/results_' subject{subjectindex(nSet)} '_' Frequency{Frequenzitter} '_Envelope_' subject{subjectindex(nSet)} '_ssmooth_zscore' ]);

        data= bsxfun(@rdivide,...
            bsxfun(@minus,TMP.ImageGridAmp,mean(TMP.ImageGridAmp,1)),std(TMP.ImageGridAmp,[],1));
        
        count=count+1;
        data = icatb_remove_mean(data);
        
        data = data'; % should be timepoints by sources
        % Temporal ica
        if (count == 1)
            startT = 1;
        else
            startT = sum(diffTimePoints(1:count-1)) + 1;
        end
        endT = sum(diffTimePoints(1:count));
        
        tc{nSet} = icasig(startT:endT, :);
        SpatialMap=corr( tc{nSet}, data);
        
        
    
        
        kernelMat.ImageGridAmp =SpatialMap';
        kernelMat.Comment      = [subject{nSet} '_' Frequency{Frequenzitter} '_Envelope_TC_cor_data'  Condition ' ONe run'];
        kernelMat.sRate      = 1;
        kernelMat.ImageGridTime         =1:size(SpatialMap,1);
        kernelMat.DataFile=[];
        kernelMat.Time         = 1:size(SpatialMap,1);
        kernelMat.SurfaceFile=['@default_subject\tess_cortex_pial_low.mat'];
        kernelMat.HeadModelType='surface';
        
        % Output filename
        OPTIONS.ResultFile = fullfile([brainstorm_db '/data/Group_analysis/Group_ICA_Brookes_based_' Condition '_23Subjects_MEG/'], ['results_' subject{nSet} '_' Frequency{Frequenzitter} 'Brookes_GroupICA_backprojected_TC_corr_ICA_' Condition '_oneRun'] );
        % Save file
        save(OPTIONS.ResultFile, '-struct', 'kernelMat');
        %         end
        
    end
    
    SpatialMap_group=corr( icasig, Dat_concat');

    kernelMat.ImageGridAmp =SpatialMap_group';
    kernelMat.Comment      = ['Groupresult_' Frequency{Frequenzitter} '_Envelope_TC_cor_data' Condition ' One Run'];
    kernelMat.sRate      = 1;
    kernelMat.ImageGridTime         =1:size(SpatialMap_group,1);
    kernelMat.DataFile=[];
    kernelMat.Time         = 1:size(SpatialMap_group,1);
    kernelMat.SurfaceFile=['@default_subject\tess_cortex_pial_low.mat'];
    kernelMat.HeadModelType='surface';
    
    % Output filename
    OPTIONS.ResultFile = fullfile([brainstorm_db '/data/Group_analysis/Group_ICA_Brookes_based_' Condition '_23Subjects_MEG/'], ['results_Group_' Frequency{Frequenzitter} 'Brookes_GroupICA_backprojected_TC_corr_' Condition '_one_Run']) ;
    % Save file
    save(OPTIONS.ResultFile, '-struct', 'kernelMat');
    
     end
     clear icasig
    
end


