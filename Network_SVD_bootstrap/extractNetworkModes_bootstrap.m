function [networkModesFile] = extractNetworkModes_bootstrap(avgCorrCoefFile, bootstrap_index,frequency,varargin)
%function [networkModesFile] = extractNetworkModes(avgCorrCoefFile, avgRandCorrCoefFile)
% Computes the main PCA modes of the connectivity patterns in the array
% stored in avgCorrFile. If avgRandCorrFile is defined, the signal modes
% are projected away from the first mode of the noise connectivity patterns
% found from the array in avgRandCorrCoefFile.


% Written by the Neurospeed & MEG Program at McGill
% Contact: sylvain.baillet@mcgill.ca
% 2012

% Esther Florin 2020
% FWHM should be in m. Therefore if FWHM is 7mm it should be defined as 0.007 
% Works if sourcefiles are within the dabase. Otherwise the surface file
% cannot be found.
% Added bootstrap to the script.


iSave = 0; %Save intermediate files iSave = 1

nDims = 1175; % Number of dimensions used to compute principal modes.

if nargin > 2
    avgRandCorrCoefFile = varargin{1};
    isnoiseProj = 1;
else
    isnoiseProj = 0;
end

connArray = load(avgCorrCoefFile);

% Define basis for connectivity similarity analysis
% Here = random selection of nDims cortical vertices
Scout = round(linspace(1,size(connArray.ImageGridAmp,2),nDims));

% Consider only absolute values of correlation coefs, which denote connectivity
%connArray.ImageGridAmp = abs(connArray.ImageGridAmp); % Comment if using
%z-score maps 

% Average connectivity pattern; for the record
MeanConn = mean(connArray.ImageGridAmp, 2);

% Generate results structure where to later save the results
clusterResults = connArray; 
clusterResults.ImageGridAmp = []; 

% Compute connectivity similarity metrics
ttmp = connArray.ImageGridAmp * connArray.ImageGridAmp(Scout, :)';
%ttmp = normr( bsxfun(@minus, ttmp , mean(ttmp, 2) ) );

if isnoiseProj
    % Projection away from noise modes
    
    sprintf(' - Compute projector(s) away from principal noise connectivity pattern(s)')
    
    % Same data prep as for data time series
    noise = load(avgRandCorrCoefFile, 'ImageGridAmp');
  
    noise_ttmp = noise.ImageGridAmp * noise.ImageGridAmp(Scout, :)';
    
    %Compute singular values and vectors to identify principal noise modes
    [noiseU S V] = svd(noise_ttmp' * noise_ttmp , 0); clear V
    S = diag(S);
    countS = cumsum(S)/sum(S);clear S; 
    
    disp(sprintf('-- First singular values, noise modes:'))
    countS(1:10)'
    
    nprojNoise = 1; %max(find(countS<=.9)); % was set to 1
    
    if isempty(nprojNoise)
        nprojNoise = 1;
    end
   
    % Just keep first nprojNoise principal noise modes(s)
    noiseU = noiseU(:,1:nprojNoise);
    
    % Generate orthogonal projector away from principal noise mode(s)
    noisePattern = noise_ttmp * noiseU;
    noiseProj = eye(length(Scout)) - noiseU * noiseU';
   
    % Project noise corrcoef maps away from nprojNoise noise modes
    noise_ttmp = noise_ttmp * noiseProj;
%     noise_ttmp = normr(noise_ttmp); 
        noise_ttmp = normmatrix(noise_ttmp); 

    if iSave
        % And save for the record
        clusterResults.ImageGridAmp = noise_ttmp;
        clusterResults.Comment = ...
            sprintf('All subjects | NOISE MODES | %d noise modes removed | %d modes', nprojNoise, nDims);
        clusterResults.Time =  1:size(clusterResults.ImageGridAmp, 2);
        clusterResults.ImageGridTime = 1:size(clusterResults.ImageGridAmp, 2);
        clusterFile = strrep(avgCorrCoefFile,'results_', 'results_AllSubjects_0000_noiseModes.mat');
        save(clusterFile, '-struct', 'clusterResults')
        clusterResults.ImageGridAmp = [];
    end
  
    if iSave
        % Project away noise correlation coefficients from their principal
        % mode(s) and save residuals
        Proj = normmatrix(noisePattern')';
        Proj = eye(size(Proj,1)) - Proj * Proj';
        %     clusterResults.ImageGridAmp = bsxfun(@minus, noise.ImageGridAmp,...
        %         mean(noise.ImageGridAmp, 2));
        
        clusterResults.ImageGridAmp = Proj * noise.ImageGridAmp;
        % And save
        clusterResults.Comment = ...
            sprintf('All subjects | NOISE RESIDUALS | %d noise modes removed', nprojNoise);
        clusterResults.Time =  1:size(clusterResults.ImageGridAmp, 2);
        clusterResults.ImageGridTime = 1:size(clusterResults.ImageGridAmp, 2);
        clusterFile = strrep(avgCorrCoefFile,'results_', 'results_AllSubjects_0100_projectedNoiseCorrCoef.mat');
        save(clusterFile, '-struct', 'clusterResults')
        clusterResults.ImageGridAmp = [];
    end
    
else % No projection away from noise
    
    noiseProj = 1; noisePattern =  []; nprojNoise = 0;

end

% Project away signal corrcoef from noise modes
ttmp = ttmp * noiseProj;

disp(sprintf('-- First singular values, signal modes:'))
[U S V] = svd( ttmp'*ttmp ,0); clear V
S = diag(S);
countS = cumsum(S)/sum(S);
countS(1:10)'

nproj = 5; %max(find(countS<.995)); % was set to 5

if isempty(nproj)
    nproj = 1;
end

% Keep only first nproj signal modes
Proj = U(:,1:nproj);
% Compute orthogonal projectors away from these modes, to study residuals
Proj = eye(size(Proj,1)) - Proj * Proj';

% Projections
U2 = ttmp * U(:,1:nproj);
ttmp2 = ttmp * Proj;

%Save mean connectivity pattern, residuals, and maps of first nproj modes
clusterResults.ImageGridAmp = [MeanConn, sqrt(sum((ttmp2.*ttmp2), 2)),  U2];

clusterResults.Comment = sprintf('All subjects | Mean Conn. pattern, Residuals, SIGNAL MODES | %d noise modes, %d signal modes', nprojNoise, nproj);
clusterResults.Comment= [clusterResults.Comment ' | ' frequency '|' num2str(bootstrap_index)];
clusterResults.Time =  1:size(clusterResults.ImageGridAmp, 2);
clusterResults.ImageGridTime = 1:size(clusterResults.ImageGridAmp, 2);
networkModesFile = strrep(avgCorrCoefFile,'results_', 'results_AllSubjects_100_signalModes');
networkModesFile=[networkModesFile '_' frequency '_' num2str(bootstrap_index) '.mat'];
save(networkModesFile, '-struct', 'clusterResults')


% Now save positive and negative valued mode maps separately
ntmp = U2;
ptmp = U2; clear U2

% Keep only negative mode values
ntmp(ntmp>0) = 0;
% And rectify for display
ntmp = abs(ntmp);

% Keep only positive mode values
ptmp(ptmp<0) = 0;

%And save
clusterResults.ImageGridAmp = [ntmp, ptmp];
clusterResults.ImageGridAmp =  bsxfun(@rdivide, clusterResults.ImageGridAmp , max(abs(clusterResults.ImageGridAmp), [] , 1));

clusterResults.Comment = sprintf('All subjects | positive & negative, SIGNAL MODES | %d signal modes', nproj);
clusterResults.Comment= [clusterResults.Comment ' | ' frequency '|' num2str(bootstrap_index)];
clusterResults.Time =  1:size(clusterResults.ImageGridAmp, 2);
clusterResults.ImageGridTime = 1:size(clusterResults.ImageGridAmp, 2);
networkModesFile = strrep(avgCorrCoefFile,'results_', 'results_AllSubjects_000_posandnegSignalModes.mat');
networkModesFile=[networkModesFile '_' frequency '_' num2str(bootstrap_index) '.mat'];

save(networkModesFile, '-struct', 'clusterResults')


[tmp, itmp] = max(clusterResults.ImageGridAmp,[],2);
indices = unique(itmp);
tmp = zeros(size(clusterResults.ImageGridAmp,1), length(indices));
for k = 1:length(indices)
    tmp(itmp==indices(k),k) = k;
end
clusterResults.ImageGridAmp = tmp;
clusterResults.ImageGridAmp = [itmp, tmp];

clusterResults.Comment = sprintf('All subjects | SIGNAL MODES LABELS | %d modes', nproj);
clusterResults.Comment= [clusterResults.Comment ' | ' frequency '|' num2str(bootstrap_index)];
clusterResults.Time =  1:size(clusterResults.ImageGridAmp, 2);
clusterResults.ImageGridTime = 1:size(clusterResults.ImageGridAmp, 2);
networkModesFile = strrep(avgCorrCoefFile,'results_', 'results_AllSubjects_001_labelsSignalModes.mat');
networkModesFile=[networkModesFile '_' frequency '_' num2str(bootstrap_index) '.mat'];

save(networkModesFile, '-struct', 'clusterResults')



