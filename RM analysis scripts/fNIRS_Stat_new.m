%written by Jun Yin

function statOutput = fNIRS_Stat(inputFileName, timeResolution,timeRange, subjName,groupFile)

% inputFileName: group matlab file for all subjects, must be defined
% timeResolution: the time resolution for averaging, must be defined
% timeRange:  the time length for computing, option;default: the time after 0
% subjName:  your interested subjects, option; default: all subjects
% groupFile: the brain region divided according to channels, option; default: no grouping



%==============================================================
% detect the nargin, and load matlab file
%==============================================================
close all
if nargin < 2
    disp('Error! Please input necessary variables, including file name of matlab data and time resolution for averaging');
end



%============================================================
% load matlab file and reorganize the data structure
%==========================================================

load(inputFileName);
dataStat.MeasList = group.procInput.procInputSubj.procInputRun.SD.MeasList(1:end,1:2);

%--------difine time information
dataStat.time.all = group.procResult.tHRF;
if isempty(timeRange)
    dataStat.time.anaIndex = find(group.procResult.tHRF>=0);
else
    dataStat.time.anaIndex = [(find(group.procResult.tHRF==timeRange(1))):find(group.procResult.tHRF==timeRange(2))];
end
timeFramesEpoch = timeResolution*round(1/(group.procResult.tHRF(2)-group.procResult.tHRF(1))); % how many frames should be computed 19/05/202 edited this manually to 23 because due to the rounding to 24 frames it ended up 3 frames short for the 5th time bin.

%------define conditions
%please change according to your design
dataStat.condition.names = {'Upright' 'Inverted'}; % please change according to names, 
conditionROI = [1 8]; %only be able to analyse two conditions, hence if more than 2, please choose which two conditions
legendLabel = {'Upright','Inverted','P-Upright','P-Inverted','P-Diff'}; % change accordingly


dataStat.condition.subjIndex = group.conditions.CondRunIdx;

%-----define subjects
if isempty(subjName)
    dataStat.subjInfo.number = group.nFiles;
    for i=1:dataStat.subjInfo.number
        dataStat.subjInfo.fileName{1,i} = group.subjs(1,i).name;
    end
else
    subjName = sort(subjName);
    dataStat.subjInfo.number = length(subjName);
    for i=1:dataStat.subjInfo.number
        if subjName(i)<10
            dataStat.subjInfo.fileName{1,i} = ['BabyFace' '0' num2str(subjName(i)) ]; % suffix shoud be changed accordingly
        else
            dataStat.subjInfo.fileName{1,i} = ['BabyFace' num2str(subjName(i))] ; % suffix shoud be changed accordingly
        end
    end
end

%------define data of each subject
subjCount = 0;
for i=1:length(group.subjs)
    for n=1:dataStat.subjInfo.number
        if strcmp(dataStat.subjInfo.fileName{1,n},group.subjs(1,i).name) == 1
            subjCount = subjCount + 1;
            dataStat.subjData(1,subjCount).name = group.subjs(1,i).name;
            dataStat.subjData(1,subjCount).dcAvgOld = group.subjs(1,i).procResult.dcAvg; %1-D: time; 2-D: signal types; 3-D: channels; 4-D: trials
            j = 0;
            
            % time points
            dataTemp = [];
            while j < (length(dataStat.time.anaIndex)/timeFramesEpoch-1)
                j = j + 1;
                indexStart = dataStat.time.anaIndex(1) +(j-1)*timeFramesEpoch;
                dataTemp(j,:,:,:) = mean(group.subjs(1,i).procResult.dcAvg([indexStart:(indexStart+timeFramesEpoch-1)],:,:,:),1);
                %1-D: time; 2-D: signal types; 3-D: channels; 4-D: conditions
                
            end
            
            % channels
            if isempty(groupFile)
                dataStat.subjData(1,subjCount).group.file = [1:size(dataStat.MeasList,1)];
                dataStat.subjData(1,subjCount).group.codeUsed = [1:size(dataStat.MeasList,1)];
                dataStat.subjData(1,subjCount).group.code = [dataStat.MeasList dataStat.subjData(1,subjCount).group.file'];
                dataStat.subjData(1,subjCount).dcAvgNew = dataTemp;
            else
                groupCodeTemp = importdata(groupFile);
                groupCode = groupCodeTemp.data;
                
                groupNum = unique(groupCode(:,3));
                dataStat.subjData(1,subjCount).group.code = groupCode;
                dataStat.subjData(1,subjCount).group.file = groupNum;
                count = 0;
                for groupLoop = 1:length(groupNum)
                    if groupNum(groupLoop)==0
                        continue;
                    else
                        count = count + 1;
                        dataStat.subjData(1,subjCount).group.codeUsed(count) = groupNum(groupLoop);
                        
                        if size(mean(dataTemp(:,:,find(groupCode(:,3)==groupNum(groupLoop)),:),3),1)==1
                           dataStat.subjData(1,subjCount).dcAvgNew(1,:,count,:) = mean(dataTemp(:,:,find(groupCode(:,3)==groupNum(groupLoop)),:),3);
                        else
                           dataStat.subjData(1,subjCount).dcAvgNew(:,:,count,:) = mean(dataTemp(:,:,find(groupCode(:,3)==groupNum(groupLoop)),:),3);
                        end
                        
                    end
                end
                
                
            end % end if for channels
        end % end if for subjData
    end % end for
     
    
end % end for 

%========================================
% statistics part
%=======================================

[timeN signalN channelN conditionN] = size(dataStat.subjData(1,1).dcAvgNew);
conditionN = length(conditionROI);
%for signalCount=1:signalN
    for channelCount=1:channelN
        for timeCount=1:timeN
            for conditionCount=1:conditionN
                for subjCount=1:dataStat.subjInfo.number                
                    dataStat.allSubjData.HbO(subjCount,conditionCount,timeCount,channelCount) = dataStat.subjData(1,subjCount).dcAvgNew(timeCount,1,channelCount,conditionROI(conditionCount));
                    dataStat.allSubjData.HbR(subjCount,conditionCount,timeCount,channelCount) = dataStat.subjData(1,subjCount).dcAvgNew(timeCount,2,channelCount,conditionROI(conditionCount));
                    dataStat.allSubjData.HbT(subjCount,conditionCount,timeCount,channelCount) = dataStat.subjData(1,subjCount).dcAvgNew(timeCount,3,channelCount,conditionROI(conditionCount));
                    % 1-D: subject; 2-D: conditions; 3-D: time; 4-D: channels
                end % end subj
                % do one-sample t test
                [H, dataStat.allSubjData.HbOTestOneSampleP(timeCount,channelCount,conditionCount), CI, STAT] = ttest(dataStat.allSubjData.HbO(:,conditionCount,timeCount,channelCount),0,0.05,'both');
                dataStat.allSubjData.HbOTestOneSampleT(timeCount,channelCount,conditionCount) = STAT.tstat;
                [H, dataStat.allSubjData.HbRTestOneSampleP(timeCount,channelCount,conditionCount), CI, STAT] = ttest(dataStat.allSubjData.HbR(:,conditionCount,timeCount,channelCount),0,0.05,'both');
                dataStat.allSubjData.HbRTestOneSampleT(timeCount,channelCount,conditionCount) = STAT.tstat;
                [H, dataStat.allSubjData.HbTTestOneSampleP(timeCount,channelCount,conditionCount), CI, STAT] = ttest(dataStat.allSubjData.HbT(:,conditionCount,timeCount,channelCount),0,0.05,'both');
                dataStat.allSubjData.HbTTestOneSampleT(timeCount,channelCount,conditionCount) = STAT.tstat;
            end % end condition
            % do paired t test
            
            [H, dataStat.allSubjData.HbOTestPairedP(timeCount,channelCount), CI, STAT] = ttest(dataStat.allSubjData.HbO(:,1,timeCount,channelCount),dataStat.allSubjData.HbO(:,2,timeCount,channelCount),0.05,'both');
            dataStat.allSubjData.HbOTestPairedT(timeCount,channelCount) = STAT.tstat;
            [H, dataStat.allSubjData.HbRTestPairedP(timeCount,channelCount), CI, STAT] = ttest(dataStat.allSubjData.HbR(:,1,timeCount,channelCount),dataStat.allSubjData.HbR(:,2,timeCount,channelCount),0.05,'both');
            dataStat.allSubjData.HbRTestPairedT(timeCount,channelCount) = STAT.tstat;
            [H, dataStat.allSubjData.HbTTestPairedP(timeCount,channelCount), CI, STAT] = ttest(dataStat.allSubjData.HbT(:,1,timeCount,channelCount),dataStat.allSubjData.HbT(:,2,timeCount,channelCount),0.05,'both');
            dataStat.allSubjData.HbTTestPairedT(timeCount,channelCount) = STAT.tstat;
        end % end time
    end % end channel
%end % end signal

outputFileName = inputFileName(1:end-4);
save([outputFileName 'AVG.mat'], 'dataStat');
statOutput = dataStat;

% output txt file, include t and p value for each group 


%=======================================
% plot statistics results
%=======================================
timePointsN = size(dataStat.allSubjData.HbO(:,:,:,:),3 );
figureHbO = figure('name', 'Results of HbO2 (Oxy)');
set(figureHbO, 'position', get(0,'ScreenSize'));
figureHbR = figure('name', 'Results of HHb (De-Oxy)');
set(figureHbR, 'position', get(0,'ScreenSize'));
figureHbT = figure('name', 'Results of HbT (Total)');
set(figureHbT, 'position', get(0,'ScreenSize'));

%po = get(figureHbO, 'Position' );


groupFinal = dataStat.subjData(1,1).group.codeUsed;
% HbO showing
figure(figureHbO);
xtick = [1:timePointsN];

for groupCount=1:length(groupFinal)
    subplot(ceil(channelN/5),5,groupCount);
    [figureHbOHandle(groupCount,:) graphH1 graphH2] = plotyy(xtick',[reshape(nanmean(dataStat.allSubjData.HbO(:,1,:,groupCount),1),timePointsN,1), ...
        reshape(nanmean(dataStat.allSubjData.HbO(:,2,:,groupCount),1),timePointsN,1)],...
        xtick',[reshape(dataStat.allSubjData.HbOTestOneSampleP(:,groupCount,1),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbOTestOneSampleP(:,groupCount,2),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbOTestPairedP(:,groupCount),timePointsN,1)]);
    set(figureHbOHandle(groupCount,1),'XColor','k','YColor','b');
    set(figureHbOHandle(groupCount,2),'XColor','k','YColor','r');
    set(get(figureHbOHandle(groupCount,1),'YLabel'),'String','\Delta Conc');
    set(get(figureHbOHandle(groupCount,2),'YLabel'),'String','P-Value');
    set(figureHbOHandle(groupCount,1),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbOHandle(groupCount,1),'YLim',[-.5E-6 .5E-6],'YTick',[-.5E-6 0 .5E-6]);
    set(figureHbOHandle(groupCount,2),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbOHandle(groupCount,2),'YTick',[0:0.5:1]);
    title(['HbO2 Channel-' num2str(groupFinal(groupCount)) ' (' num2str(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbO(:,1,1,groupCount)))) ')']);
    set(graphH1,'LineWidth',3,'Marker','o','MarkerSize',5)
    set(graphH2(1:2),'Marker','square','MarkerSize',5,'LineWidth',3,'LineStyle','--')
    set(graphH2(3),'Marker','^','MarkerSize',5,'LineWidth',3,'LineStyle',':')
    
    hold(figureHbOHandle(groupCount,1),'all');
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbO(:,1,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbO(:,1,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbO(:,1,1,groupCount)))),'Parent',figureHbOHandle(groupCount,1), 'Color',[0 0 0]);
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbO(:,2,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbO(:,2,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbO(:,2,1,groupCount)))),'Parent',figureHbOHandle(groupCount,1), 'Color',[0 0 0]);
    hold(figureHbOHandle(groupCount,2),'all');
    %bar(figureHbO, 'FaceColor',[0.90 0.90 0.90],'BarWidth',10,'Parent',figureHbOHandle(1,2));
    plot([xtick(1)-0.5:xtick(end)+0.5],0.05*ones(1,length([xtick(1)-0.5:xtick(end)+0.5])),'Parent',figureHbOHandle(groupCount,2),'LineWidth',1.5, 'Color',[0.50 0.50 0.50]);%'LineColor',[0.90 0.90 0.90]);
end
legend([graphH1' graphH2'],legendLabel,'Location',[0.85 0.1 0.1 0.15]);



% HbR showing
figure(figureHbR);
xtick = [1:timePointsN];
for groupCount=1:length(groupFinal)
    subplot(ceil(channelN/5),5,groupCount);
    [figureHbRHandle(groupCount,:) graphH1 graphH2] = plotyy(xtick',[reshape(nanmean(dataStat.allSubjData.HbR(:,1,:,groupCount),1),timePointsN,1), ...
        reshape(nanmean(dataStat.allSubjData.HbR(:,2,:,groupCount),1),timePointsN,1)],...
        xtick',[reshape(dataStat.allSubjData.HbRTestOneSampleP(:,groupCount,1),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbRTestOneSampleP(:,groupCount,2),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbRTestPairedP(:,groupCount),timePointsN,1)]);
    set(figureHbRHandle(groupCount,1),'XColor','k','YColor','b');
    set(figureHbRHandle(groupCount,2),'XColor','k','YColor','r');
    set(get(figureHbRHandle(groupCount,1),'YLabel'),'String','Amplitude');
    set(get(figureHbRHandle(groupCount,2),'YLabel'),'String','P-Value');
    set(figureHbRHandle(groupCount,1),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbRHandle(groupCount,1),'YLim',[-.5E-6 .5E-6],'YTick',[-.5E-6 0 .5E-6]);
    set(figureHbRHandle(groupCount,2),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbRHandle(groupCount,2),'YTick',[0:0.5:1]);
    title(['HHb Channel-' num2str(groupFinal(groupCount)) ' (' num2str(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbR(:,1,1,groupCount)))) ')']);
    set(graphH1,'LineWidth',3,'Marker','o','MarkerSize',5)
    set(graphH2(1:2),'Marker','square','MarkerSize',5,'LineWidth',3,'LineStyle','--')
    set(graphH2(3),'Marker','^','MarkerSize',5,'LineWidth',3,'LineStyle',':')
    
    hold(figureHbRHandle(groupCount,1),'all');
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbR(:,1,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbR(:,1,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbR(:,1,1,groupCount)))),'Parent',figureHbRHandle(groupCount,1), 'Color',[0 0 0]);
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbR(:,2,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbR(:,2,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbR(:,2,1,groupCount)))),'Parent',figureHbRHandle(groupCount,1), 'Color',[0 0 0]);
    hold(figureHbRHandle(groupCount,2),'all');
    %bar(figureHbR, 'FaceColor',[0.90 0.90 0.90],'BarWidth',10,'Parent',figureHbRHandle(1,2));
    plot([xtick(1)-0.5:xtick(end)+0.5],0.05*ones(1,length([xtick(1)-0.5:xtick(end)+0.5])),'Parent',figureHbRHandle(groupCount,2),'LineWidth',1.5, 'Color',[0.50 0.50 0.50]);%'LineColor',[0.90 0.90 0.90]);
end
legend([graphH1' graphH2'],legendLabel,'Location',[0.85 0.1 0.1 0.15]);

% HbT showing
figure(figureHbT);
xtick = [1:timePointsN];
for groupCount=1:length(groupFinal)
    subplot(ceil(channelN/5),5,groupCount);
    [figureHbTHandle(groupCount,:) graphH1 graphH2] = plotyy(xtick',[reshape(nanmean(dataStat.allSubjData.HbT(:,1,:,groupCount),1),timePointsN,1), ...
        reshape(nanmean(dataStat.allSubjData.HbT(:,2,:,groupCount),1),timePointsN,1)],...
        xtick',[reshape(dataStat.allSubjData.HbTTestOneSampleP(:,groupCount,1),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbTTestOneSampleP(:,groupCount,2),timePointsN,1), ...
        reshape(dataStat.allSubjData.HbTTestPairedP(:,groupCount),timePointsN,1)]);
    set(figureHbTHandle(groupCount,1),'XColor','k','YColor','b');
    set(figureHbTHandle(groupCount,2),'XColor','k','YColor','r');
    set(get(figureHbTHandle(groupCount,1),'YLabel'),'String','Amplitude');
    set(get(figureHbTHandle(groupCount,2),'YLabel'),'String','P-Value');
    set(figureHbTHandle(groupCount,1),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbTHandle(groupCount,1),'YLim',[-.5E-6 .5E-6],'YTick',[-.5E-6 0 .5E-6]);
    set(figureHbTHandle(groupCount,2),'XLim',[xtick(1)-0.5 xtick(end)+0.5],'XTick',xtick);
    set(figureHbTHandle(groupCount,2),'YTick',[0:0.5:1]);
    title(['HbT Channel-' num2str(groupFinal(groupCount)) ' (' num2str(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbT(:,1,1,groupCount)))) ')']);
    set(graphH1,'LineWidth',3,'Marker','o','MarkerSize',5)
    set(graphH2(1:2),'Marker','square','MarkerSize',5,'LineWidth',3,'LineStyle','--')
    set(graphH2(3),'Marker','^','MarkerSize',5,'LineWidth',3,'LineStyle',':')
    
    hold(figureHbTHandle(groupCount,1),'all');
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbT(:,1,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbT(:,1,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbT(:,1,1,groupCount)))),'Parent',figureHbTHandle(groupCount,1), 'Color',[0 0 0]);
    errorbar(xtick',reshape(nanmean(dataStat.allSubjData.HbT(:,2,:,groupCount),1),timePointsN,1),reshape(nanstd(dataStat.allSubjData.HbT(:,2,:,groupCount),0,1),timePointsN,1)/sqrt(dataStat.subjInfo.number-sum(isnan(dataStat.allSubjData.HbT(:,2,1,groupCount)))),'Parent',figureHbTHandle(groupCount,1), 'Color',[0 0 0]);
    hold(figureHbTHandle(groupCount,2),'all');
    %bar(figureHbT, 'FaceColor',[0.90 0.90 0.90],'BarWidth',10,'Parent',figureHbTHandle(1,2));
    plot([xtick(1)-0.5:xtick(end)+0.5],0.05*ones(1,length([xtick(1)-0.5:xtick(end)+0.5])),'Parent',figureHbTHandle(groupCount,2),'LineWidth',1.5, 'Color',[0.50 0.50 0.50]);%'LineColor',[0.90 0.90 0.90]);
end
legend([graphH1' graphH2'],legendLabel,'Location',[0.85 0.1 0.1 0.15]);




    
    