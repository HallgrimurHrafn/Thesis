
txtName={'-5x24s-100mVrms.txt','-5x24s-0.5Vrms.txt','-5x24s-1Vrms.txt','-5x24s-1.5Vrms.txt'};
for n=1:4
    [f{n},impedance{n},phase{n}]=getDataInfo(txtName{n});
end
%%
figure(1)
hold on
for n=1:4
    plot(f{n},phase{n})
end
legend('100m','0.5','1','1.5')
set(gca,'XScale','log')
hold off
ylabel('Phase [Degrees]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',13);
set(gca,'LineWidth',1.5);
fig=gcf;
fig.Position=[680   558   640   420];
ax=gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];



figure(2)
hold on
for n=1:4
    ftemp=f{n};
    imptemp=impedance{n};
    plot(ftemp(ftemp>10 & ftemp<8000),imptemp(ftemp>10 & ftemp<8000))
end
legend('100m','0.5','1','1.5')
set(gca,'XScale','log', 'YScale', 'log')
hold off
ylabel('Impedance [\Omega]')
xlabel('Frequency [Hz]')
set(gca,'FontSize',13);
set(gca,'LineWidth',1.5);
fig=gcf;
fig.Position=[680   558   640   420];
ax=gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax_width = outerpos(3) - ti(1) - ti(3);
ax.Position = [left bottom ax_width ax_height];

%%

txtName1='Current-5x24s-1Vrms.txt';
txtName2='Voltage-5x24s-1Vrms.txt';
[voltData,fs]=ReadPulse(txtName2);
[currData,~]=ReadPulse(txtName1);

signalLength=find(isnan(voltData(:,2)),1)-1;
voltage=voltData(1:signalLength,2)-mean(voltData(1:signalLength,2));
current=currData(1:signalLength,2)-mean(currData(1:signalLength,2));

figure(3)
plot(voltage)
ylabel('Measured Voltage [V]')
xlabel('Sample [n]')
set(gca,'FontSize',13);
set(gca,'LineWidth',1.5);
fig=gcf;
fig.Position=[680   558   640   420];
ax=gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

figure(4)
plot(current)
ylabel('Measured Current [A]')
xlabel('Sample [n]')
set(gca,'FontSize',13);
set(gca,'LineWidth',1.5);
fig=gcf;
fig.Position=[680   558   640   420];
ax=gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];




function [f,impedance,phase]=getDataInfo(txtName)
    txtName1=strcat('Voltage',txtName);
    txtName2=strcat('Current',txtName);
    [voltData,fs]=ReadPulse(txtName1);
    [currData,~]=ReadPulse(txtName2);
    
    signalLength=find(isnan(voltData(:,2)),1)-1;
    voltage=voltData(1:signalLength,2)-mean(voltData(1:signalLength,2));
    current=currData(1:signalLength,2)-mean(currData(1:signalLength,2));

    [tf, f]=tfestimate(current,voltage,[],[],[],fs);
    impedance=smooth(abs(tf),100);
    phase=smooth(unwrap(angle(tf))*180/pi,200);

    ind=find(f==8,1):find(f==2600,1);
    f=f(ind);
    impedance=impedance(ind);
    phase=phase(ind);
end


function [Data, fs]=ReadPulse(txtName)
% open file
    fid = fopen(txtName);
% pass through first 38 lines
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
    Numbers=textscan(fid,'%s %s %s %s');
%     convert all the decimal commas to periods and convert the answer to
%     double.
    num1=str2double(strrep(Numbers{2},',','.'));
    num2=str2double(strrep(Numbers{3},',','.'));
    num3=str2double(strrep(Numbers{4},',','.'));
%     close document
    fclose(fid);
%     combine numerical data
    Data=[num1,num2,num3];
    
end