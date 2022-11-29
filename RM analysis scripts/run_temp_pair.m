%written by Jun Yin, will save the output of the fNIRS_Stat script in statOutput

%%Input for fNIRS stat is as follows: 
% FileName: group matlab file for all subjects, must be defined
% timeResolution: the time resolution for averaging, must be defined
% timeRange:  the time length for computing, option;default: the time after 0
% subjName:  your interested subjects, option; default: all subjects
% groupFile: the brain region divided according to channels, option; default: no grouping
statOutput = fNIRS_Stat_new('groupResults_NIRSFace.mat', 3,[],[],'Channels.txt');