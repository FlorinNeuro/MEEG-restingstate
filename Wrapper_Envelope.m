
% With this script the data in sensor space are loaded and transformed to
% source space based on a pre-defined Imaging Kernel in brainstorm. It
% produces matlab files (101) per subject, which contain the envelope for
% each vertex in the defined frequency bins. Those need to be put back
% together so that it represents the whole brain again and can be viewed in
% brainstorm

% Written by Esther Florin 2018
clear
subject= {'S001_MEG' {'run2' 'run2' 'run3'};...
    };
condition='';
frequency_bins=[1 4;4 8; 8 13; 13 30; 30 50];
blocks=100;

database='/gpfs/project/florine/resting_young/';

pathStr = '/gpfs/project/florine/brainstorm3/'; %path to brainstorm
addpath(genpath(pathStr));
pathStr = '/gpfs/project/florine/hilbert-ica-rsn/'; %path to envelope scripts
addpath(genpath(pathStr));

OutputPath=['/home/florine/Output/'];

for isubject=1:size(subject,1) %loop over subjects
for sourcediv=1:blocks %loop over source_blocks to make it tracktable for the size of the data
    Bad_tmp=[];
    Data=[];
   
    
    for irun=1:length(subject{isubject,2}) % loop over runs
     iCond=[];
          sRate=[];

    % Loading the data
    
    a=dir(fullfile([database 'data/' subject{isubject,1} '/'  subject{isubject,2}{irun}  '/']));
    % loading the imaging Kernel from brainstorm
    for k=1:length(a)
        if ~isempty(strfind(a(k).name, 'results_PNAI')) || ~isempty(strfind(a(k).name, 'results_wMNE_MEG_GRAD_MEG_MAG_KERNEL_'))
            % This depends on the brainstorm version and source reconstruction technique used. Check in data base.
            Kernel = load ([database 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/' a(k).name]);
        end
    end
    
    for i=1:length(a)
        if ~isempty(strfind(a(i).name, 'block')) % Check in your database the name of the imported data
            iCond(end+1) = i;
        end
    end
    if size(iCond)>1
        disp('There were more than one datafile')
        quit
    end
     if isempty(iCond)
        disp('There was no datafile')
        quit
    end
    for i = 1:length(iCond)
        data=  load ([database '/data/'  subject{isubject,1} '/' subject{isubject,2}{irun} '/' a(iCond(i)).name]);
        if ~isempty(data)
            sRate = round(abs(1 / (data.Time(2) - data.Time(1)))); %new sampling rate
            if sourcediv<blocks+1 %check if last block, because that one is smaller and do the transformation of the loaded data into source-space
            data.F=Kernel.ImagingKernel((sourcediv-1)* floor(size(Kernel.ImagingKernel,1)/blocks)+1:sourcediv* floor(size(Kernel.ImagingKernel,1)/blocks),:)*data.F(Kernel.GoodChannel,:);            
            else
            data.F=Kernel.ImagingKernel((sourcediv-1)* floor(size(Kernel.ImagingKernel,1)/blocks)+1:end,:)*data.F(Kernel.GoodChannel,:);
            end
            for i=1:length(data.Events)
           
                if strcmp(data.Events(1,i).label,'BAD')
                    Bad_tmp=(data.Events(1,i).samples)-data.Time(1)*sRate;
                    
                end
            end
            
            % Remove bad segments from time series and change the events accordingly
            BAD=ones(1,size(data.F,2));
            if ~isempty(Bad_tmp)
            for ll=1:size(Bad_tmp,2)
                if (Bad_tmp(1,ll)-5)>5 && (Bad_tmp(2,ll)+5)<size(data.F,2)
                    BAD(Bad_tmp(1,ll)-5:Bad_tmp(2,ll)+5)=0;
                elseif (Bad_tmp(1,ll)-5)<=5
                    BAD(1:Bad_tmp(2,ll)+5)=0;
                elseif (Bad_tmp(2,ll)+5)>=size(data.F,2)
                    BAD(Bad_tmp(1,ll)-5:end)=0;
                end
            end
            % removes mean block by block
            
            for ll=2:size(Bad_tmp,2)
                if (Bad_tmp(1,ll)-5)>5 && (Bad_tmp(2,ll)+5)<size(data.F,2)
                    meanBaseline = mean(data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5), 2);
                    % Remove mean
                    data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5) =  bsxfun(@minus, data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5), meanBaseline);

                end
            end
            
            meanBaseline = mean(data.F(:,1:Bad_tmp(1,1)-5), 2);
            % Remove mean
            data.F(:,1:Bad_tmp(1,1)-5) =  bsxfun(@minus, data.F(:,1:Bad_tmp(1,1)-5), meanBaseline);
            meanBaseline = mean(data.F(:,Bad_tmp(2,size(Bad_tmp,2))+5:end), 2);
            % Remove mean
            data.F(:,Bad_tmp(2,size(Bad_tmp,2))+5:end) =  bsxfun(@minus, data.F(:,Bad_tmp(2,size(Bad_tmp,2))+5:end), meanBaseline);
            
        end
            meanBaseline = mean(data.F(:,BAD==1), 2);
            
            % Remove mean over all good segments. Not necessary
            data.F =  bsxfun(@minus, data.F, meanBaseline);
            Data=[Data data.F(:,BAD==1)];
        else
            sprintf(['missing data in ' a(iCond(i)).name])
        end
        
    end
    clear Kernel data
    end
%     HilbertEnvelope_wideBand_iir(Data,sRate,frequency_bins, [subject{isubject,1} '_' condition '_' num2str(sourcediv)], OutputPath);
     HilbertEnvelope_smallband_fir(Data,sRate,frequency_bins, [subject{isubject,1} '_' condition '_' num2str(sourcediv)], OutputPath);
% HilbertEnvelope_smallband_iir(Data,sRate,frequency_bins, [subject{isubject,1} '_' condition '_' num2str(sourcediv)], OutputPath);
% HilbertEnvelope_wideBand_fir(Data,sRate,frequency_bins, [subject{isubject,1} '_' condition '_' num2str(sourcediv)], OutputPath);
end

clear Data data sRate meanBaseline BAD Bad_tmp Trigger_Bad
end
