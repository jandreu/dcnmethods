%written by Jun Yin
% save script
%save signal of HbO and HbR
filename = 'statResults-Amplitudes_NIRSFace_19072022.txt'; % specify output name
[subjN conditionN timeN channelN] = size(statOutput.allSubjData.HbO);
dataHbOAvg = statOutput.allSubjData.HbO; % changed from mean to nanmean
dataHbRAvg = statOutput.allSubjData.HbR; % changed from mean to nanmean
fid = fopen(filename,'w');
for s=1:subjN
    for i=1:channelN
        %HbO
        for j=1:timeN % time window
            for n=1:conditionN
                fprintf(fid,'%6.5f\t',dataHbOAvg(s,n,j,i)*(10^7));
            end
            fprintf(fid,'%6.5f\t',(dataHbOAvg(s,1,j,i)*(10^7)-dataHbOAvg(s,2,j,i)*(10^7)));
        end
        fprintf(fid,'\n');
        %HbR
        for j=1:timeN
            for n=1:conditionN
                fprintf(fid,'%6.5f\t',dataHbRAvg(s,n,j,i)*(10^7));
            end
            fprintf(fid,'%6.5f\t',(dataHbRAvg(s,1,j,i)*(10^7)-dataHbRAvg(s,2,j,i)*(10^7)));
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
close all


