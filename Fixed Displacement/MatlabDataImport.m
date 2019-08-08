% addpath
addpath(pwd)

% load written data as base
writtenData=load('WrittenData.mat');
totalStruct=writtenData.startingStruct;

cd('FixedDisplacement\Up\')
for n=1:length(totalStruct)
% for n=6
%     generate folder name based on Data1 vector, frequency folder
    Position=string(totalStruct(n).MeasNum);
%     enter said folder
    cd(Position)
%     get data
    [v,c,fs]=PulseInterpreter(totalStruct(n));
    totalStruct(n).v=v;
    totalStruct(n).c=c;
    totalStruct(n).fs=fs;
    
    freq=totalStruct(n).freq;
    dvolt=totalStruct(n).driveVolts;
    
    freq(:,isnan(dvolt))=[];
    dvolt(:,isnan(dvolt))=[];
    
    totalStruct(n).freq=freq;
    totalStruct(n).driveVolts=dvolt;
    
    cd ..
end

cd ..\..

OrderByUp=[21,20,1:19];
totalStruct=totalStruct(OrderByUp);