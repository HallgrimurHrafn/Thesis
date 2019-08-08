%% Load Data

% get voltage-R*I, and R
load('R_v-ir-allNorm.mat')
% get voltage, current, derivative of current and displacements
load('singleLoopDataAllNorm.mat')
% get frequency and voltage data
load('FixedDisplacementInStruct-allNorm.mat')

%% Split between frequencies.

% get data on current, derivative of current and displacement from
% singleloopdata struct
current=SingleLoopData.current;
dcurrent=SingleLoopData.Dcurrent;
displacement=SingleLoopData.displacement;

% initialize storage for current
current10=zeros(21,16,8192);
current100=zeros(21,3,8192);
current1000=zeros(21,3,8192);

% initialize storage for derivative of current
dcurrent10=zeros(21,16,8192);
dcurrent100=zeros(21,3,8192);
dcurrent1000=zeros(21,3,8192);

% initialize storage for V-IR
vri10=zeros(21,16,8192);
vri100=zeros(21,3,8192);
vri1000=zeros(21,3,8192);

% initialize storage for DC voltage
vDC10=zeros(21,16);
vDC100=zeros(21,3);
vDC1000=zeros(21,3);

% initialize storage for voltage rms
Vrms10=zeros(21,16);
Vrms100=zeros(21,3);
Vrms1000=zeros(21,3);

% initialize storage for driving voltages
DriveVolt10=zeros(21,16);
DriveVolt100=zeros(21,3);
DriveVolt1000=zeros(21,3);

R10=zeros(21,16);
R100=zeros(21,3);
R1000=zeros(21,3);

for m=1:21
    h100Count=0;
    h1000Count=0;
    freq=totalStruct(m).freq;
    volt=totalStruct(m).driveVolts;
    for f=1:16
        current10(m,f,:)=current(m,f,:);
        dcurrent10(m,f,:)=dcurrent(m,f,:);
        vri10(m,f,:)=v_c(m,f,:);
        DriveVolt10(m,f)=volt(f);
        R10(m,f)=R(m,f);
        vDC10(m,f)=vDC(m,f);
        Vrms10(m,f)=Vrms(m,f);
    end
    
    for f=17:length(freq)
        if freq(f)==100
            h100Count=h100Count+1;
            current100(m,h100Count,:)=current(m,f,:);
            dcurrent100(m,h100Count,:)=dcurrent(m,f,:);
            vri100(m,h100Count,:)=v_c(m,f,:);
            DriveVolt100(m,h100Count)=volt(f);
            R100(m,h100Count)=R(m,f);
            vDC100(m,h100Count)=vDC(m,f);
            Vrms100(m,h100Count)=Vrms(m,f);
        else
            h1000Count=h1000Count+1;
            current1000(m,h1000Count,:)=current(m,f,:);
            dcurrent1000(m,h1000Count,:)=dcurrent(m,f,:);
            vri1000(m,h1000Count,:)=v_c(m,f,:);
            DriveVolt1000(m,h1000Count)=volt(f);
            R1000(m,h1000Count)=R(m,f);
            vDC1000(m,h1000Count)=vDC(m,f);
            Vrms1000(m,h1000Count)=Vrms(m,f);
            
        end
    end 
end
% displacement as a vector.
displacement=displacement(:,1);


% fixing the order of 100hz and 1000hz to match [0.3,1,4.5]
% 100hz
correctOrderFor100=CorrectOrder(DriveVolt100);
DriveVolt100=reOrder(DriveVolt100,correctOrderFor100);
current100=reOrder(current100,correctOrderFor100);
dcurrent100=reOrder(dcurrent100,correctOrderFor100);
vri100=reOrder(vri100,correctOrderFor100);
R100=reOrder(R100,correctOrderFor100);
vDC100=reOrder(vDC100,correctOrderFor100);
Vrms100=reOrder(Vrms100,correctOrderFor100);

% 1khz
correctOrderFor1000=CorrectOrder(DriveVolt1000);
DriveVolt1000=reOrder(DriveVolt1000,correctOrderFor1000);
current1000=reOrder(current1000,correctOrderFor1000);
dcurrent1000=reOrder(dcurrent1000,correctOrderFor1000);
vri1000=reOrder(vri1000,correctOrderFor1000);
R1000=reOrder(R1000,correctOrderFor1000);
vDC1000=reOrder(vDC1000,correctOrderFor1000);
Vrms1000=reOrder(Vrms1000,correctOrderFor1000);

% % plotting test to ensure/validify correct reordering.
% figure(1)
% plot(squeeze(current100(1,:,:))')
% hold on
% plot(squeeze(test3(1,:,:))')
% legend('1or','2or','3or','1t','2t','3t')
% hold off

% Errors in 100hz and 1khz measurement at measurement number 2,15,21
% found through looking at measured R values, those showed approximately R=1,
% likely short as the added resistor was 1 ohm.
noErrorInMeasurement=(1:21);
noErrorInMeasurement([2,15,21])=[];


h10=struct;
h10.displacement=displacement;
h10.current=current10;
h10.dcurrent=dcurrent10;
h10.v_ri=vri10;
h10.driveVolt=DriveVolt10;
h10.vDC=vDC10;
h10.Vrms=Vrms10;
h10.R=R10;
h10.freq=10;

h100=struct;
h100.displacement=displacement(noErrorInMeasurement,:,:);
h100.current=current100(noErrorInMeasurement,:,:);
h100.dcurrent=dcurrent100(noErrorInMeasurement,:,:);
h100.v_ri=vri100(noErrorInMeasurement,:,:);
h100.driveVolt=DriveVolt100(noErrorInMeasurement,:,:);
h100.vDC=vDC100(noErrorInMeasurement,:);
h100.Vrms=Vrms100(noErrorInMeasurement,:);
h100.R=R100(noErrorInMeasurement,:,:);
h100.freq=100;

h100=removeRepeatedMeasurement(h100,3,2); % two 0.3 Volt measurements in a row


h1000=struct;
h1000.displacement=displacement(noErrorInMeasurement,:,:);
h1000.current=current1000(noErrorInMeasurement,:,:);
h1000.dcurrent=dcurrent1000(noErrorInMeasurement,:,:);
h1000.v_ri=vri1000(noErrorInMeasurement,:,:);
h1000.driveVolt=DriveVolt1000(noErrorInMeasurement,:,:);
h1000.vDC=vDC1000(noErrorInMeasurement,:);
h1000.Vrms=Vrms1000(noErrorInMeasurement,:);
h1000.R=R1000(noErrorInMeasurement,:,:);
h1000.freq=1000;


h1000=removeRepeatedMeasurement(h1000,17,2); % two 0.3 volt measurements in a row.
h1000=removeRepeatedMeasurement(h1000,6,2); %accidental extra 10hz instead of 1khz

save('frequencySplitStruct-AllNorm.mat','h10','h100','h1000')


function order=CorrectOrder(fix)
    fitMe=[0.3,1,4.5];
    order=zeros(size(fix));
    for n=1:size(fix,1)
        inds=perms([1,2,3]);
        [~,correctInd]=min(sum((perms(fix(n,:))~=fitMe),2));
        order(n,:)=inds(correctInd,:);
    end
end

function fixed=reOrder(fixme,correctOrder)
    for n=1:21
        for m=1:3
            fixed(n,m,:)=fixme(n,correctOrder(n,m),:);
        end
    end
end






function structInUse=removeRepeatedMeasurement(structInUse,m,n)
    structInUse.current(m,n,:)=0*structInUse.current(m,n,:);
    structInUse.v_ri(m,n,:)=0*structInUse.v_ri(m,n,:);
    structInUse.driveVolt(m,n)=0;
    structInUse.vDC(m,n)=0;
    structInUse.Vrms(m,n,:)=0;
    structInUse.R(m,n)=0;
end






