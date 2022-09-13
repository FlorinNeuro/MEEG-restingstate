function [avgRandCorrCoefFile] = averageRandCorrCoef_all_MEG_elekta_bootstrap(dataPath, OutputPath,sourceFiles,nTime, FWHM)
% Generate surrogate correlation coefficient arrays from random sensor data time series
% to assess spurious spatial correlations due to the spatial resolution of
% MEG imaging. The random time series are low-pass filtered and sampled
% equivalently to the megBOLD and megPAC signal models
% ZavgRandCorrCoefFile contains the z-scores of each row of the average
% corrcoef array
% ZavgCorrCoefFile contains the z-scores of the data average
% corrcoef array, each row standardized with respect to the mean and std of
% the rows of the random corrcoef array


% Written by the Neurospeed & MEG Program at McGill
% Contact: sylvain.baillet@mcgill.ca
% 2012

% Esther Florin 2020
% FWHM should be in m. Therefore if FWHM is 7mm it should be defined as 0.007 
% Works if sourcefiles are within the dabase. Otherwise the surface file
% cannot be found.
% Added bootstrap to the script.


iSave = 1;
% Load ImagingKernel for template brain

load('results_MN_MEG_GRAD_MEG_MAG_KERNEL_180716_1021.mat') % inv

inv.ImagingKernel=ImagingKernel;
inv.Whitener=Whitener;
inv.nComponents=nComponents;

% Generate randomized time series of source series
% with characteristic akin to megPAC and megBOLD signal models
% and compute corresponding correlation coefficients




tmp = load(fullfile(dataPath, sourceFiles{1}));

% Sampling rate of low-passed versions of source series

sRate = 1./abs(diff(tmp.Time(1:2))); % Hz 


% Generate surrogate random time series
tmp.ImageGridAmp = randn(size(inv.ImagingKernel,2) , nTime );
% Low-pass filter below sRate<3
tmp.ImageGridAmp = bst_bandpass_fft(tmp.ImageGridAmp, sRate, 0, sRate/3, 1, 0);

% Compute corresponding source maps
tmp.ImageGridAmp  = inv.ImagingKernel * tmp.ImageGridAmp;
% apply spatial smoothing
tmp2=strfind(dataPath,'data');
database=dataPath(1:tmp2-1);
SurfaceMat=load(fullfile(database,'anat',tmp.SurfaceFile));

cortS.tri = SurfaceMat.Faces;
cortS.coord = SurfaceMat.Vertices';
% Get the average edge length
[vi,vj] = find(SurfaceMat.VertConn);
Vertices = SurfaceMat.VertConn;
meanDist = mean(sqrt((Vertices(vi,1) - Vertices(vj,1)).^2 + (Vertices(vi,2) - Vertices(vj,2)).^2 + (Vertices(vi,3) - Vertices(vj,3)).^2));
% FWHM in surfstat is in mesh units: Convert from millimeters to "edges"
FWHMedge = FWHM ./ meanDist;
tmp.ImageGridAmp=SurfStatSmooth(tmp.ImageGridAmp', cortS, FWHMedge)';


tmp.ImageGridAmp  = corrcoef( tmp.ImageGridAmp' );


if iSave
    tmp.Comment = sprintf('All Subjects | average | random corr maps');
    tmp.Time =  1:size(tmp.ImageGridAmp, 2);
    tmp.ImageGridTime = 1:size(tmp.ImageGridAmp, 2);
    avgRandCorrCoefFile = strrep(sourceFiles{1},'results_', 'results_000_NoiseAverageCorrmap');
    avgRandCorrCoefFile = fullfile(OutputPath,avgRandCorrCoefFile);
    save(avgRandCorrCoefFile, '-struct', 'tmp')
else
    avgRandCorrCoefFile = '';
end


