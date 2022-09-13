
% With this script the data in sensor space are loaded and transformed to
% source space based on a pre-computed Imaging Kernel in brainstorm. It
% produces matlab files (101) per subject, which contain the PAC matrix for
% each vertex in the defined frequency bins. From those the low-frequency
% of the maximum can be determined and used for the constructrion of the
% megPAC signal --> see Wrapper_megPAC

% Written by Esther Florin 2018
clear
subject= {'S001_MEG' {'run2' 'run2' 'run3'};...
    };
flow1 = [2 30]; %low frequency range
gamma_lowend=80; % lower end of the gamma range
gamma_low_high_end=150; %higer end of the gamma range
blocks=100;

database='resting_young/'; % path to brainstorm data base

pathStr = '/brainstorm3/'; %path to the brainstorm
addpath(genpath(pathStr));
pathStr = '/PAC_megPAC/'; %path to PAC_megPAC scripts
addpath(genpath(pathStr));

OutputPath=['/home/florine/Output/'];

for isubject=1:size(subject,1) %loop over subjects
for sourcediv=1:blocks+1 %loop over source_blocks
    Trigger_Bad=[];
    Data=[];
    iCond=[];
    
    for irun=1:length(subject{isubject,2}) % loop over runs
    iCond=[];
    % Loading the data
    
    a=dir(fullfile([database 'data/' subject{isubject,1} '/'  subject{isubject,2}{irun}  '/']));
    % loading the imaging Kernel
    for k=1:length(a)
        if ~isempty(strfind(a(k).name, 'results_PNAI_')) % this searach field needs to be adjusted according to the ImagingKernel used!
            Kernel = load ([database 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/' a(k).name]);
        end
    end
    
    for i=1:length(a)
        if ~isempty(strfind(a(i).name, 'block')) % this searach field needs to be adjusted according to a unique identifier in the data used!
            iCond(end+1) = i;
        end
    end
    if length(iCond)>1
        disp('There were more than one datafile')
        quit
    end
    for i = 1:length(iCond)
        data=  load ([database '/data/'  subject{isubject,1} '/' subject{isubject,2}{irun} '/' a(iCond(i)).name]);
        if ~isempty(data)
            sRate = round(abs(1 / (data.Time(2) - data.Time(1)))); %new sampling rate
            if sourcediv<blocks+1 %check if last block, because that one is smaller
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
            
            for ll=1:length(Bad_tmp)
                if (Bad_tmp(1,ll)-5)>5 && (Bad_tmp(2,ll)+5)<size(data.F,2)
                    BAD(Bad_tmp(1,ll)-5:Bad_tmp(2,ll)+5)=0;
                elseif (Bad_tmp(1,ll)-5)<=5
                    BAD(1:Bad_tmp(2,ll)+5)=0;
                elseif (Bad_tmp(2,ll)+5)>=size(data.F,2)
                    BAD(Bad_tmp(1,ll)-5:end)=0;
                end
            end
            % removes mean block by block
            
            for ll=2:length(Bad_tmp)
                if (Bad_tmp(1,ll)-5)>5 && (Bad_tmp(2,ll)+5)<size(data.F,2)
                    meanBaseline = mean(data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5), 2);
                    % Remove mean
                    data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5) =  bsxfun(@minus, data.F(:,Bad_tmp(2,ll-1)+5:Bad_tmp(1,ll)-5), meanBaseline);

                end
            end
            
            meanBaseline = mean(data.F(:,1:Bad_tmp(1,1)-5), 2);
            % Remove mean
            data.F(:,1:Bad_tmp(1,1)-5) =  bsxfun(@minus, data.F(:,1:Bad_tmp(1,1)-5), meanBaseline);
            meanBaseline = mean(data.F(:,Bad_tmp(2,length(Bad_tmp))+5:end), 2);
            % Remove mean
            data.F(:,Bad_tmp(2,length(Bad_tmp))+5:end) =  bsxfun(@minus, data.F(:,Bad_tmp(2,length(Bad_tmp))+5:end), meanBaseline);
            
            
            meanBaseline = mean(data.F(:,BAD==1), 2);
            
            % Remove mean over all good segments. Not necessary
            data.F =  bsxfun(@minus, data.F, meanBaseline);
            Data=[Data data.F(:,BAD==1)];
        else
            sprintf(['missing data in ' a(iCond(i)).name])
        end
    clear Kernel
        
    end
    end
PACestimate_cluster(Data,sRate,flow1,[gamma_lowend gamma_low_high_end],[OutputPath subject{isubject,1} '_PAC_rest_' num2str(sourcediv)])
end
clear Data data sRate meanBaseline BAD Bad_tmp Trigger_Bad
end
