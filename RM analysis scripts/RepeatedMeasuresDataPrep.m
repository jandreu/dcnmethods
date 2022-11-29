% This script will convert the statResults-Amplitude files produced by Jun
% Jin's fNIRS analysis scripts to a format the can be analysed using RM
% analyses in SPSS; questions? email c.deklerk@essex.ac.uk

clear all

%input basic aspects of design
chan=20; %number of channels
bins=1; %number of bins
chromophores=2; % Oxy and Deoxy
conditions=3; % Jun's scripts would always output three conditions: the two you wanted to compare and their difference
rowsperpart=chromophores * chan;

%read in datafile
datafile = 'statResults-Amplitudes_NIRSFace_19072022.xlsx'; 
[num, txt]= xlsread(datafile);
Nparticipants=(length(txt)-2)/chan; %to calculate the number of participants
Nparticipants=Nparticipants/2; %take length of txt file subtract 2 and divide by 2 (oxy and deoxy in file)

k=3;
for i=1:Nparticipants
    participant(i)= txt(k,1);
    k=k+rowsperpart;
end
participant=participant'; %transpose participant variable so that it is vertically instead of horizontally organised

chancounter=1;
bincounter=1;
realchancounter=1;
realcondcounter=1;
realbincounter=1;
condition={'Upright', 'Inverted', 'Diff'};
chromophore={'Hbo', 'Hbb'};
condcounter=1;
chromocounter=1;
%variablename=cell(chan*bins*conditions*chromophores);

%create the variable names
for j=1:chan*bins*conditions*chromophores %2484
    variablename(j)= strcat('bin',num2str(realbincounter),'_chan',num2str(realchancounter),chromophore(chromocounter),condition(realcondcounter));
    chancounter=chancounter+1;
    if chancounter==3 % to add 1 to the channel counter every other loop (there is oxy and deoxy data for each channel)
        realchancounter=realchancounter+1;
        if realchancounter>chan 
            realchancounter=1;
        end
        chancounter=1;
    end
    
    chromocounter=chromocounter+1; %to create a HbO and Hbb variable for each channel
    if chromocounter==3
        chromocounter=1;
    end
    
    condcounter=condcounter+1;
    if condcounter>rowsperpart % move to the next condition column once you have a variable for each row 
        realcondcounter=realcondcounter+1; 
        condcounter=1;
        if realcondcounter>3
            realcondcounter=1;
        end
    end
    
    bincounter=bincounter+1;
    if bincounter>3*rowsperpart % move to the next bin once you have a variable for each row for each of the three conditions
        realbincounter=realbincounter+1;
        bincounter=1;
    end
    
end

startpoint=1;
endpoint=chan*chromophores; %92
for p=1:Nparticipants
    columncount1=1;
    columncount2=chan*chromophores;
    for n= 3:bins*conditions+2 %for row 3 till the last row of data
        row=num(startpoint:endpoint,n); % e.g. get data from row 1:92 for column 3 
        transposedrow=row'; % transpose these data so that they are horizontal
        data(p,columncount1:columncount2)=transposedrow(1,1:chan*chromophores); %save horizontal data next to each other in file; new line for each participant
        columncount1=columncount1+chan*chromophores;
        columncount2=columncount2+chan*chromophores;
    end
    startpoint=startpoint+chan*chromophores;
    endpoint=endpoint+chan*chromophores;
end 

save('dataforSPSS.mat','data');





