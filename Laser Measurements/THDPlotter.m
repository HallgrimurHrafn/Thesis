%% load and set up
clear all
close all
load('LaserMeasurementSingleFreqResults.mat')
%% Reorder based on growing VoltageLevels
fnames=fieldnames(baseStruct);
orderedBaseStruct=struct;
for n=1:length(fnames)
    structInUse=baseStruct.(fnames{n});
    cells=struct2cell(structInUse);
    sortvals = squeeze(cells(3,1,:));
    amp=squeeze(cells(4,1,:));
    for m=1:length(structInUse)
        sortvals{m}(end-4:end)='';
        sortvals{m}=str2double(sortvals{m});
    end
    sortvals=cell2mat(sortvals);
    amp=cell2mat(amp);
    volts=sortvals.*amp;
    [sortedVolts,ix]=sort(volts);
    structInUse=structInUse(ix);
    for m=1:length(structInUse)
        structInUse(m).volts=sortedVolts(m);
    end
    orderedBaseStruct.(fnames{n})=structInUse;
end


%% Calculate THD for di/dt and dx/dt


nmax=10;
for n=1:length(fnames)
    structInUse=orderedBaseStruct.(fnames{n});
    hz=fnames{n};
    hz=str2double(hz(4:end));
    for m=1:length(structInUse)
        current=structInUse(m).current(:,2);
        displacement = structInUse(m).displacement(:,2);
        fs=structInUse(m).fs;
        [dcurrSpec,dcurrTHD,dcurrF]=getDerivativeTHD(current,fs,hz,nmax);
        structInUse(m).dcurrSpec=dcurrSpec;
        structInUse(m).dcurrTHD=dcurrTHD;
        structInUse(m).dcurrF=dcurrF;
        [velSpec,velTHD,velF]=getDerivativeTHD(displacement,fs,hz,nmax);
        structInUse(m).velSpec=velSpec;
        structInUse(m).velTHD=velTHD;
        structInUse(m).velF=velF;
    end
    orderedBaseStruct.(fnames{n})=structInUse;
end

%% plot the THD levels

figure(1)
hold on
for n=1:length(fnames)
    structInUse=orderedBaseStruct.(fnames{n});
    v=zeros(1,length(structInUse));
    velTHD=zeros(1,length(structInUse));
    for m=1:length(structInUse)
        v(m)=structInUse(m).volts;
        velTHD(m)=structInUse(m).velTHD;
    end
    plot(v/1000,velTHD)
end
hold off
legend('20Hz','50Hz','80Hz','120Hz','250Hz')
xlabel('Driving Voltage [Vrms]')
ylabel('Acceleration THD  [%]')

figure(2)
hold on
for n=1:length(fnames)
    structInUse=orderedBaseStruct.(fnames{n});
    v=zeros(1,length(structInUse));
    dcurrTHD=zeros(1,length(structInUse));
    for m=1:length(structInUse)
        v(m)=structInUse(m).volts;
        dcurrTHD(m)=structInUse(m).dcurrTHD;
    end
    plot(v/1000,dcurrTHD)
end
hold off
legend('20Hz','50Hz','80Hz','120Hz','250Hz')
xlabel('Driving Voltage [Vrms]')
ylabel('Second Derivative of I w.r.t. time THD [%]')

%% Save ordered Struct
save('LaserMeasurementResultsOrdered.mat','orderedBaseStruct')


%%

function [Y,THD,f]=getDerivativeTHD(signal,fs,f0,nmax)
%     get signal length
    N=length(signal);
%     get frequency domain of signal
    fsignal=fft(signal);
%     only care about the first half
    fsignal=fsignal(1:N/2);
%     get frequency vector
    f=fs*(0:(N/2-1))/N;
%     radian vector
    w=2*pi*f;
%     get second derivative of signal
    Y=(1i*w').^2.*fsignal;

%     initialize denominator
    nom=0;
%     for all harmonic orders we examine
    for m=2:nmax
%         order indexes by approximity to a given harmony of frequency f0
        [~,inds] = sort(abs(f-m*f0)); 
%         add three closest frequencybins squared to denominator, as per
%         thd formula
    try
        nom= nom+sum(abs(Y(inds(1:3))).^2);
    catch
        nom;
    end
    end
    [~,inds] = sort(abs(f-f0));
    denom=sum(abs(Y(inds(1:3))).^2);
    THD = sqrt(nom)/sqrt(denom)*100;
    
end