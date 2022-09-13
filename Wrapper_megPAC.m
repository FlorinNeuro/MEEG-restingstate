% With this script the data in sensor space are loaded and transformed to
% source space based on a pre-computed Imaging Kernel in brainstorm. It
% produces matlab files (101) per subject, which contain the megPACPAC matrix for
% each vertex in the defined frequency bins. Those need to be put back
% together so that it represents the whole brain again and can be viewed in
% brainstorm
clear
subject= {'S001_MEG' {'run2' 'run2' 'run3'};...
    };

flow1 = [2:0.5:12 14:2:30];
gamma_lowend=80;
gamma_low_high_end=150;
blocks=100;
PAC_path='/home/florine/PAC/'; %path were the PAC matlab files are saved, which have been computed with the Wrapper_PAC
database='/resting_young/'; %path to brainstorm data base
pathStr = '/brainstorm3/'; %path to brainstorm
addpath(genpath(pathStr));
pathStr = 'PAC_megPAC'; %path to PAC_megPAC scripts
addpath(genpath(pathStr));

OutputPath=['/home/florine/Output/'];


for isubject=1:size(subject,1) %loop over subjects
    for sourcediv=1:blocks+1 %loop over source_blocks
        Trigger_Bad=[];
        Data=[];
        Kernel =[];
        
        
        for irun=1:length(subject{isubject,2}) % loop over runs
            iCond=[];
            % Loading the data
            
            a=dir(fullfile([database 'data/' subject{isubject,1} '/'  subject{isubject,2}{irun}  '/']));
            % loading the imaging Kernel
            Kernel =[];
            for k=1:length(a)
                if ~isempty(strfind(a(k).name, 'results_PNAI_'))  % this searach field needs to be adjusted according to the ImagingKernel used!
                    Kernel = load ([database 'data/' subject{isubject,1} '/' subject{isubject,2}{irun} '/' a(k).name]);
                end
            end
            
            for i=1:length(a)
                if ~isempty(strfind(a(i).name, 'block'))   % this searach field needs to be adjusted according to a unique identifier in the data used!
                    iCond(end+1) = i;
                end
            end
            if size(iCond)>1
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
                    
                    
                    meanBaseline = mean(data.F(:,BAD==1), 2);
                    
                    % Remove mean over all good segments. Not necessary
                    data.F =  bsxfun(@minus, data.F, meanBaseline);
                    Data=[Data data.F(:,BAD==1)];
                else
                    sprintf(['missing data in ' a(iCond(i)).name])
                end
                
            end
        end
        % clear Kernel nach letztem run muß gelöscht werden.
        %determine maxPAC for all sources
        load([PAC_path '/' subject{isubject,1} '_PAC_rest_' num2str(sourcediv) '.mat']);
        ind=find(flow1>=2);
        flow=flow1(ind);
        for PAC_ind=1:size(directPAC,1)
            pactmp=squeeze(directPAC(PAC_ind,ind,:));
            [max_pac,tmp]=max(pactmp,[],2);
            [~,ind2]=max(max_pac);
            
            lowfrequency(PAC_ind)=flow(ind2);
        end
        megPAC_bst( Data, sRate, [subject{isubject,1} '_sd_' num2str(sourcediv)], lowfrequency, Kernel.GoodChannel,OutputPath )
    end %loop over sources
end %loop over subjects