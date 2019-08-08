function [v,c,fs]=PulseInterpreter(StructInUse)
    [voltage,fs]=ReadPulse('Expanded Time Voltage.txt');
    [Current,~]=ReadPulse('Expanded Time Current.txt');

%     size of which measurements can be taken from
    frameSize=ceil(9.5*fs);
%     headroom for transient start
    skipsizeFront=floor(0.35*fs);
%     room at end and 5 sek delay for next signal
    skipSizeDelay=ceil(5.15*fs);
%     function to simplify layout. given nth frame in measurement, finds
%     the delay before it happens, including the front headroom start.
    delay=@(n) n*skipsizeFront+(n-1)*skipSizeDelay+(n-1)*frameSize;
    
%     initialize matrices to store data based on frames
    v=zeros(frameSize,length(StructInUse.driveVolts));
    c=zeros(frameSize,length(StructInUse.driveVolts));
%     retrieve every frame
    for n=1:21
        v(:,n)=voltage(1+delay(n):frameSize+delay(n),2);
        c(:,n)=Current(1+delay(n):frameSize+delay(n),2);
    end
    
%     perform the same for additional measurement if exist
    if StructInUse.add
        [addVoltage, ~]=ReadPulse('addVoltage.txt');
        [addCurrent, ~]=ReadPulse('addVoltage.txt');
        for n=22:length(StructInUse.driveVolts)
            v(:,n)=addVoltage(1+delay(n-21):frameSize+delay(n-21),2);
            c(:,n)=addCurrent(1+delay(n-21):frameSize+delay(n-21),2);
        end
    end
    
%   remove NaN parts, failed measurement frames  
    v(:,isnan(StructInUse.driveVolts))=[];
    c(:,isnan(StructInUse.driveVolts))=[];
    
    
%   normalize current, resting point was not 0.
    restCurr=Current(ceil(10.25*fs):ceil(14.75*fs),2);
    c=c-mean(restCurr);
    
%     normalize voltage, in case resting point was not 0.
    restVolt=voltage(ceil(10.25*fs):ceil(14.75*fs),2);
    v=v-mean(restVolt);
    
end


% Function designed to convert the txt format from polytechs output
% to matlab compatible data, in the form of a 2xn matrix. also outputs
% check to ensure that for every folder we are getting both current and
% displacement, not just two of either and 0 of the other.
function [Data, fs]=ReadPulse(txtName)
% open file
    fid = fopen(txtName);
% pass through first 39 lines
    for n=1:38
        fgets(fid);
    end
%     get fs line
    fs=fgets(fid);
%     capture the numbers found in the string.
    fs=str2double(strrep(fs(min(regexp(fs,'\d')):max(regexp(fs,'\d'))),',','.'));
    
%     skip 45 lines
    for n=1:45
        fgets(fid);
    end
    
%     retrieve all the numerical data as strings
    Numbers=textscan(fid,'%s %s %s');
%     convert all the decimal commas to periods and convert the answer to
%     double.
    num1=str2double(strrep(Numbers{2},',','.'));
    num2=str2double(strrep(Numbers{3},',','.'));
%     close document
    fclose(fid);
%     combine numerical data
    Data=[num1,num2];
    
end