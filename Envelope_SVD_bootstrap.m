% Within this file resting state networks are extracted based on the
% envelope data on each vertex within the defined frequencies (Frequency). 
% Across subjects a bootstrap is performed and the resting state networks for each repetition are saved. 
% Input:    The subjects IDs as stored in brainstorm
%           The Condition, where the files are saved und in the brainstorm
%           database
%           The brainstorm path and the path to the megPAC srcipts for
%           network extraction
pathStr = '/gpfs/project/florine/brainstorm3/'; %path to brainstorm
addpath(genpath(pathStr));
pathStr = '/gpfs/project/florine/hilbert-ica-rsn/'; %path to megPAC srcipts 
% for network extraction
addpath(genpath(pathStr));
FWHM =0.005; % Spatial smoothing (7mm)
Condition='wide_band_fir';
brainstorm_db='/gpfs/project/florine/resting_young/';

dataPath = [brainstorm_db '/data/Group_analysis/Envelope_' Condition '/'];
OutputPath='/home/florine/Output/';

subject={'S001_MEG'; 'S002_MEG';'S003_MEG'; 'S004_MEG'; 'S006_MEG'; 'S007_MEG'; 'S008_MEG'; 'S009_MEG'; 'S010_MEG'; 'S011_MEG'; 'S012_MEG'; 'S013_MEG'; 'S014_MEG'; ...
    'S015_MEG'; 'S017_MEG'; 'S019_MEG'; 'S020_MEG'; 'S021_MEG'; 'S022_MEG'; 'S023_MEG'; 'S024_MEG'; 'S025_MEG'; 'S026_MEG';};

Frequency ={'Delta'; 'Theta'; 'Alpha'; 'Beta'; 'Gamma'};
resultFiles={};
for Frequenzitter=1:length(Frequency) %loop over all frequencies
    resultFiles={};
    diffTimePoints=[];
    Dat_concat=[];
    isubjects=[];
    for subjects=1:length(subject)
        resultFiles{end+1}=  ['results_' subject{subjects} '_' Frequency{Frequenzitter} '_Envelope_' subject{subjects} '_ssmooth_zscore.mat' ];
     end
    
    resultFiles={resultFiles};
   
    t00 = tic;
    for bootstrap_index=1:100 %bootstrap repetitions
                
        t0 = tic;
        
        t1 = tic;
        [avgCorrCoefFile,nTime] = computeCorrCoef_all_boostrap(dataPath, OutputPath, resultFiles{1},FWHM);
        
        disp(sprintf('Done in %3.2f',toc(t1)))
        
        if bootstrap_index==1
            [avgRandCorrCoefFile] = averageRandCorrCoef_all_MEG_elekta_bootstrap(dataPath,OutputPath, resultFiles{1} ,nTime, FWHM);
        end
        
              
        disp(sprintf('Extracting principal modes of connectivity...'))
        t1 = tic;
    [networkModesFile] = extractNetworkModes_bootstrap(avgCorrCoefFile,bootstrap_index, [Frequency{Frequenzitter} '_' Condition], avgRandCorrCoefFile );
        
       
    end
    
    
end %loop across frequencies

return

