%This script loads data, which were computed in chunks in source space and
%puts them back together and saves them in brainstorm format
% Written by Esther Florin

clear 
GroupPath='\brainstorm_db\resting_young\data\'; % Brainstorm data base, where the data is supposed to be saved.
Inputpath='Chunks_of_data\';  %Path to the chunks which are supposed to be concatenated
Filename_append='Conv_megPAC_model_demeaned';
a=dir(fullfile(Inputpath));
Filepart='megPAC_lf_2-48'; % Here you can specify the ending of your files which make them unique

%%
subjectname={'S001' 'S002' 'S003' 'S004' 'S006' 'S007' 'S008' 'S009' 'S010' 'S011' 'S012' 'S013'};

for subjectindex=1:length(subjectname)
Conv=[];
for i = 1:101
    iCond=[];
    for k=1:length(a)
        
        if ~isempty(strfind(a(k).name,['results_' subjectname{subjectindex} '_sd_' num2str(i) Filepart])) %Here the filename is put together. Needs to be checked, if it is correct for the data used. 
            i
            iCond(end+1) = k;
        end
    end
       load([Inputpath  '/' a(iCond).name]);
    Conv=[Conv; ImageGridAmp];
end

  meanBaseline = mean(Conv, 2);
    % Compute zscore
   Conv =         bsxfun(@minus, Conv, meanBaseline);

kernelMat.ImageGridAmp = Conv;
kernelMat.Comment      = [  'megPAC_demeaned_wholebrain']; % This can be anything
kernelMat.sRate      = sRate;
kernelMat.ImageGridTime         = 1/sRate:1/sRate:size(Conv,2)/sRate;
kernelMat.DataFile=[];
kernelMat.Time         = 1:size(Conv,2);
kernelMat.GoodChannel=GoodChannel;
 kernelMat.SurfaceFile=[subjectname{subjectindex} '/tess_cortex_pial_low.mat'];
kernelMat.HeadModelType='surface';
OPTIONS.ResultFile = fullfile(GroupPath, subjectname{subjectindex}, ['/megPAC_EEG/results_'  Filename_append ] ); %Here it is important to check that is goes into the correct folter 
 % Save file
 save(OPTIONS.ResultFile, '-struct', 'kernelMat');

end