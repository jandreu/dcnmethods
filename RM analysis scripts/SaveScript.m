%written by Jun Yin
% save script
%save signal of HbO and HbR
filename = 'statResults_200522.txt';
[subjN conditionN timeN channelN] = size(statOutput.allSubjData.HbO);
dataHbOAvg = mean(statOutput.allSubjData.HbO,1);
dataHbRAvg = mean(statOutput.allSubjData.HbR,1);
fid = fopen(filename,'w');
for i=1:channelN
    %HbO
    for j=1:timeN
        for n=1:conditionN
            fprintf(fid,'%6.4f%s%5.4f%s\t',dataHbOAvg(1,n,j,i)*(10^7),'(',statOutput.allSubjData.HbOTestOneSampleP(j,i,n),')');
        end
        fprintf(fid,'%6.4f%s%5.4f%s\t',(dataHbOAvg(1,1,j,i)*(10^7)-dataHbOAvg(1,2,j,i)*(10^7)),'(',statOutput.allSubjData.HbOTestPairedP(j,i),')');
    end
    fprintf(fid,'\n');
    %HbR
    for j=1:timeN
        for n=1:conditionN
            fprintf(fid,'%6.4f%s%5.4f%s\t',dataHbRAvg(1,n,j,i)*(10^7),'(',statOutput.allSubjData.HbRTestOneSampleP(j,i,n),')');
        end
        fprintf(fid,'%6.4f%s%5.4f%s\t',(dataHbRAvg(1,1,j,i)*(10^7)-dataHbRAvg(1,2,j,i)*(10^7)),'(',statOutput.allSubjData.HbRTestPairedP(j,i),')');
    end
    fprintf(fid,'\n');
end
fclose(fid);
close all        


