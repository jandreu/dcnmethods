%written by Jun Yin
%export data to txt for excel plotting
clear all
%please define the below parameters according to your design
fileName = 'groupResults_NIRSFaceAVG.mat';
groupFile = 'Channels.txt';
outputFile = 'NIRSFace_HbO_export2txt'; % specify output file name
channelsOutput = [1:20]; 
timeRange = [0 15];
%-----------------------------------------

load(fileName);
groupCodeTemp = importdata(groupFile);
groupCode = groupCodeTemp.data;


startingPoint = find(dataStat.time.all==timeRange(1));
endingPoint = find(dataStat.time.all==timeRange(2));
timeResultionSet = 0.1; 
timeResultionDefault = 0.1;%default from CW system
if timeResultionSet<timeResultionDefault
    timeResultionSet = timeResultionDefault;
end
sampleRate = round(timeResultionSet/timeResultionDefault);

%find index of channel in meanlist
SD_PairFile = [];
for i=1:length(channelsOutput)
    codeIndex = find(groupCode(:,3)==channelsOutput(i));
    for j=1:length(codeIndex)
        SD_PairFile = [SD_PairFile;groupCode(codeIndex,1:2),channelsOutput(i)];
    end
end
fid = fopen([outputFile,'_HbO.txt'],'w');
fprintf(fid,'%s\t%s\t%s\t','subjName','timePoint','condition');
for tempCount=1:size(SD_PairFile,1)
    fprintf(fid,'%s\t',['ChannelNO_' num2str(SD_PairFile(tempCount,3)) '_HbO']);           
end % end tempCount
fprintf(fid,'\n');
conditionNum = 1; 
signalNO = 1; % 1=HbO, 2=HbR, 3=HbT
for subjCount = 1:dataStat.subjInfo.number
    
    for timeCount=startingPoint:endingPoint
        for conditionCount = 1:conditionNum
            fprintf(fid,'%s\t%6.2f\t',dataStat.subjInfo.fileName{1,subjCount},dataStat.time.all(timeCount));
            fprintf(fid,'%s\t',num2str(conditionCount));            
            for tempCount=1:size(SD_PairFile,1)
                for channelCount = 1:length(dataStat.MeasList)/2
                    if dataStat.MeasList(channelCount,1)==SD_PairFile(tempCount,1) && ... 
                            dataStat.MeasList(channelCount,2)==SD_PairFile(tempCount,2)
                        fprintf(fid,'%6.5f\t',dataStat.subjData(1,subjCount).dcAvgOld(timeCount,signalNO,channelCount,conditionCount)*10^6);
                    end
                end % end channelCount
            end %end tempCount
            fprintf(fid,'\n');
        end %end conditionCount
    end %end timeCount
end %end subjCount

fclose(fid);
close all  
clear all
            

