
% totalStruct Load
load('FixedDisplacementInStruct-allNorm.mat');




% create emtpy container base, 21 measurement runs, max 22 different driving
% voltages/freq
% measurement,  8192 samples after resampling for 1 signal loop
base=zeros(21,22,8192);
% initializing containers for current, voltage and derivative of current
currentData=base;
voltageData=base;
derivCurrentData=base;
displacement=zeros(21,1);

% for all measurement runs
for m=1:length(totalStruct)
%     retrieve current,fs, frequencies and errorOffsets
    c=totalStruct(m).c;
    fs=totalStruct(m).fs;
    freq=totalStruct(m).freq;
%     error offset is a variable on how much to skip ahead if voltage
%     change was late and bled a little bit into measurement
    errOffset=totalStruct(m).errOffset;
%     for all driving voltages/freq
    for n=1:size(c,2)
%         get clean average one loop sinusoids for each variable
        [newCurr,newVolt, newDcurr]=signalPrep(n,totalStruct(m));
%         add to storage
        currentData(m,n,:)=newCurr;
        voltageData(m,n,:)=newVolt;
        derivCurrentData(m,n,:)=newDcurr;
        
    end
%     get total displacement as seen from the smaller part, upper.
    displacement(m)=sum(totalStruct(m).up);
end

SingleLoopData=struct;
SingleLoopData.voltage=voltageData;
SingleLoopData.current=currentData;
SingleLoopData.Dcurrent=derivCurrentData;
SingleLoopData.displacement=displacement;
save('singleLoopDataAllNorm.mat','SingleLoopData')
% save('singleLoopData.mat','SingleLoopData')