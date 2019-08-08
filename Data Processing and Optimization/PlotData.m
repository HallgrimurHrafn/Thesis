%% load data
clear all
load('phaseAlignedInFreqSplitStructs-AllNorm.mat')
PlotNumber=1;
close all

%% plotting the phaseDifference found in degrees (0:360)

% 10hz
figure(PlotNumber)
plotPhaseDifference(h10)
saveas(gcf,'Figures/Phase Difference For All 10Hz measurements.png')

% 100hz
PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotPhaseDifference(h100)
saveas(gcf,'Figures/Phase Difference For All 100Hz measurements.png')

% 1000hz
PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotPhaseDifference(h1000)
saveas(gcf,'Figures/Phase Difference For All 1000Hz measurements.png')

%% plotting resistor values found

% 10hz  resistor values
PlotNumber=PlotNumber+1;
figure(PlotNumber)
PlotResistance(h10)
saveas(gcf,'Figures/Measured DC Resistance For All 10Hz measurements.png')


% 100hz resistor values
PlotNumber=PlotNumber+1;
figure(PlotNumber)
PlotResistance(h100)
saveas(gcf,'Figures/Measured DC Resistance For All 100Hz measurements.png')

% 1000hz resistor values
PlotNumber=PlotNumber+1;
figure(PlotNumber)
PlotResistance(h1000)
saveas(gcf,'Figures/Measured DC Resistance For All 1000Hz measurements.png')
%% plotting v-ri after alignment.

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotAlignedSignal(h10)
saveas(gcf,'Figures/V-RI For a Few Measurement Over Few Displacements.png')

%% plot all rms for v_riAligned

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSV_riAligned(h10)
saveas(gcf,'Figures/Measured RMS of V-IR for 10hz.png')

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSV_riAligned(h100)
saveas(gcf,'Figures/Measured RMS of V-IR for 100hz.png')

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSV_riAligned(h1000)
saveas(gcf,'Figures/Measured RMS of V-IR for 1000hz.png')

%% plot all rms for dcurr
PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurr(h10)
saveas(gcf,'Figures/Dcurr Rms 10hz.png')

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurr(h100)
saveas(gcf,'Figures/Dcurr Rms 100hz.png')

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurr(h1000)
saveas(gcf,'Figures/Dcurr Rms 1000hz.png')

%% plot distorted dcurr versions

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotFewDcurr(h10)
saveas(gcf,'Figures/Dcurr Picked.png')


%% Plot RMS(V_ri)/Rms(dcurr)
PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurrVir(h10);

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurrVir(h100);

PlotNumber=PlotNumber+1;
figure(PlotNumber)
plotRMSDcurrVir(h1000);

%% plot rms for voltages
% PlotNumber=PlotNumber+1;
% figure(PlotNumber)
% plotRMS_v(h10)
% 
% PlotNumber=PlotNumber+1;
% figure(PlotNumber)
% plotRMS_v(h100)
% 
% PlotNumber=PlotNumber+1;
% figure(PlotNumber)
% plotRMS_v(h1000)


%% wtf in v_ri rms for 1khz at 17 displacement, 0.3drive volt (1st)
% figure
% plot(squeeze(h1000.R(17,1)*h1000.current(17,1,:)))
% hold on
% plot(squeeze(h1000.R(16,1)*h1000.current(16,1,:)))
% plot(squeeze(h1000.R(18,1)*h1000.current(18,1,:)))
% legend('17','16','18')
% hold off
% title('Current * R')

% figure
% plot(squeeze(h1000.current(17,1,:)))
% hold on
% plot(squeeze(h1000.current(16,1,:)))
% plot(squeeze(h1000.current(18,1,:)))
% legend('17','16','18')
% hold off
% title('current')

% figure
% plot(squeeze(h1000.v_ri(17,1,:)))
% hold on
% plot(squeeze(h1000.v_ri(16,1,:)))
% plot(squeeze(h1000.v_ri(18,1,:)))
% legend('17','16','18')
% hold off
% title('V-RI')

%% plotting Functions

function plotComparison(structInUse,m,n)

    v_riAligned=structInUse.v_riAligned;
    v_ri=structInUse.v_ri;
    dcurr=structInUse.dcurrent;
    
    vAtemp=squeeze(v_riAligned(m,n,:));
    vtemp=squeeze(v_ri(m,n,:));
    dctemp=squeeze(dcurr(m,n,:));

    subplot(2,1,1)
    plot(vAtemp)
    hold on
    plot(vtemp)
    legend('aligned','original')
    hold off

    subplot(2,1,2)
    plot(vAtemp/max(vAtemp))
    hold on
    plot(dctemp/max(dctemp))
    legend('aligned','derivative current')
    hold off

end

function plotPhaseDifference(structInUse) 
    displacement = getConvertedDisplacement(structInUse);
    
    phaseInDegrees=structInUse.phaseDiff/8192*360;
    R=structInUse.R;
    phaseDiffDegree = phaseInDegrees-360*(phaseInDegrees>180);
    vDC=structInUse.vDC;
    
    grad=colorGradient([0,1,0],[0,0,1],size(R,1));
    
    hold on
    for m=1:size(phaseInDegrees,1)
        plot(vDC(m,R(m,:)~=0),phaseDiffDegree(m,R(m,:)~=0),'-or','MarkerFaceColor', ...
            grad(m,:) ,'color',grad(m,:),'LineWidth',1.5)
        l(m)=string(displacement(m));
    end
    hold off
    xlabel('DC Voltage Offset [V]')
    ylabel(['Phase Difference [degrees]'])
    colormap(grad)
    clbr=colorbar('Ticks',linspace(displacement(1),displacement(end),size(R,1)), 'TickLabels',l);
    ylabel(clbr,'Displacement [mm]','FontSize',13)
    caxis([displacement(1),displacement(end)])    
    fig=gcf;
    fig.Position=[400,540,820,400];
    set(gca,'FontSize',13);
    set(gca,'LineWidth',1.5);
    
    ax=gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - 2*ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
end

function PlotResistance(structInUse)
    R=structInUse.R;
    grad=colorGradient([0,1,0],[0,0,1],size(R,1));
    vDC=structInUse.vDC;
    displacement =getConvertedDisplacement(structInUse);
    hold on
    for m=1:1:size(R,1)
        plot(vDC(m,R(m,:)~=0),R(m,(R(m,:)~=0)),'-or','MarkerFaceColor',...
            grad(m,:) ,'color',grad(m,:), 'LineWidth',1.5)
        l(m)=string(displacement(m));
    end
    hold off
    xlabel('DC Voltage Offset [V]')
    ylabel('Resistance Measured In [\Omega]')
    colormap(grad)
    clbr=colorbar('Ticks',linspace(displacement(1),displacement(end),size(R,1)), 'TickLabels',l);
    ylabel(clbr,'Displacement [mm]','FontSize',12)
    caxis([displacement(1),displacement(end)]) 
    fig=gcf;
    fig.Position=[400,540,820,400];
    set(gca,'FontSize',13);
    set(gca,'LineWidth',1.5);
    
    ax=gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - 1.75*ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end


function plotAlignedSignal(structInUse)
    v_riAligned=structInUse.v_riAligned;
    grad(1,:,:)=colorGradient([0.4,0,0],[1,0,0],4);
    grad(2,:,:)=colorGradient([0,0.4,0],[0,1,0],4);
    grad(3,:,:)=colorGradient([0,0,0.4],[0,0,1],4);
    dsiplacement=getConvertedDisplacement(structInUse);
    measurements=[2,10,20];
    hold on
    for m=measurements
        for n=1:5:16
            plot(1000*squeeze(v_riAligned(m,n,:)),'color',grad(ceil(m/9),ceil(n/5),:))
        end
    end
    xlabel('Sample [n]')
    ylabel('Voltage [mV]')
    xlim([1,8192])
    l(1)=plot(nan,nan,'r');
    l(2)=plot(nan,nan,'g');
    l(3)=plot(nan,nan,'b');
    legend(l,arrayfun(@(x) ['Displacement of ',num2str(x)] ...
        ,dsiplacement(measurements),'UniformOutput',false), 'Location','south')
    hold off
%     set(gca,'FontSize',13);
%     set(gca,'LineWidth',1.5);
    
end

function plotFewDcurr(structInUse)
    dcurrent=structInUse.dcurrent;
    grad=colorGradient([0,1,0],[0,0,1],3);
    dsiplacement=getConvertedDisplacement(structInUse);
    measurements=[2,10,20];
    hold on
    for m=measurements
            plot(1000*squeeze(dcurrent(m,10,:)), 'LineWidth', 1.5, ...
                'color',grad(find(measurements==m),:))
    end
    xlabel('Sample [n]')
    ylabel('Time-Derivative of the Current [A/s]')
    xlim([1,8192])
    legend(arrayfun(@(x) ['Displacement of ',num2str(x),'mm'] ...
        ,dsiplacement(measurements),'UniformOutput',false), 'Location','south')
    hold off
%     set(gca,'FontSize',13);
%     set(gca,'LineWidth',1.5);
    
end


function plotRMSDcurrVir(hfreq)
    vrms=hfreq.Vrms;
    grad=colorGradient([0,1,0],[0,0,1],length(hfreq.displacement));
    displacement=getConvertedDisplacement(hfreq);
    dcurr=hfreq.dcurrent;
    rmsdc=zeros(size(dcurr,1),size(dcurr,2));
    vDC=hfreq.vDC;
    for m=1:size(rmsdc,1)
        for n=1:size(rmsdc,2)
            rmsdcVir(m,n)=vrms(m,n)/rms(squeeze(dcurr(m,n,:)));
        end
    end
    hold on
    for m=1:size(rmsdcVir,1)
        plot(vDC(m,vDC(m,:)~=0),1000*rmsdcVir(m,vDC(m,:)~=0),'-or', 'LineWidth',1.5...
            ,'MarkerFaceColor',grad(m,:),'color',grad(m,:))
        l(m)=string(displacement(m));
    end
    hold off
    ylabel('Root Mean Square v\_iR/(dI/dt) [mH]')
    xlabel('DC Voltage Offset [V]')
    colormap(grad);
    clbr=colorbar('Ticks',linspace(displacement(1),displacement(end),size(hfreq.R,1)), 'TickLabels',l);
    ylabel(clbr,'Displacement [mm]','FontSize',12)
    caxis([displacement(1),displacement(end)]) 
    fig=gcf;
    fig.Position=[400,540,820,400];
    set(gca,'FontSize',13);
    set(gca,'LineWidth',1.5);
    
    ax=gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - 1.75*ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end

function plotRMSDcurr(hfreq)
    grad=colorGradient([0,1,0],[0,0,1],length(hfreq.displacement));
    displacement=getConvertedDisplacement(hfreq);
    dcurr=hfreq.dcurrent;
    rmsdc=zeros(size(dcurr,1),size(dcurr,2));
    vDC=hfreq.vDC;
    for m=1:size(rmsdc,1)
        for n=1:size(rmsdc,2)
            rmsdc(m,n)=rms(squeeze(dcurr(m,n,:)));
        end
    end
    hold on
    for m=1:size(rmsdc,1)
        plot(vDC(m,vDC(m,:)~=0),rmsdc(m,vDC(m,:)~=0),'-or', 'LineWidth',1.5...
            ,'MarkerFaceColor',grad(m,:),'color',grad(m,:))
        l(m)=string(displacement(m));
    end
    hold off
    ylabel('Root Mean Square dI/dt [A/s]')
    xlabel('DC Voltage Offset [V]')
    colormap(grad);
    clbr=colorbar('Ticks',linspace(displacement(1),displacement(end),size(hfreq.R,1)), 'TickLabels',l);
    ylabel(clbr,'Displacement [mm]','FontSize',12)
    caxis([displacement(1),displacement(end)]) 
    fig=gcf;
    fig.Position=[400,540,820,400];
    set(gca,'FontSize',13);
    set(gca,'LineWidth',1.5);
    
    ax=gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - 1.75*ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end

function plotRMS_v(hfreq)
    vrms=hfreq.Vrms;
    hold on
    for m=1:size(vrms,1)
        plot(vrms(m,:))
    end
    hold off

end

function plotRMSV_riAligned(hfreq)
    grad=colorGradient([0,1,0],[0,0,1],length(hfreq.displacement));
    displacement=getConvertedDisplacement(hfreq);
    vDC=hfreq.vDC;
    v_riAligned=hfreq.v_riAligned;
    rmsVri=zeros(size(v_riAligned,1),size(v_riAligned,2));
    for m=1:size(rmsVri,1)
        for n=1:size(rmsVri,2)
            rmsVri(m,n)=rms(squeeze(v_riAligned(m,n,:)));
        end
    end
    hold on
    for m=1:size(rmsVri,1)
        plot(vDC(m,vDC(m,:)~=0),1000*rmsVri(m,vDC(m,:)~=0), ...
            '-or','MarkerFaceColor',grad(m,:),'color',grad(m,:),'LineWidth',1.5)
        l(m)=string(displacement(m));
    end
    hold off
    ylabel('Root Mean Square Voltage [mVrms]')
    xlabel('DC Voltage Offset [V]')
    colormap(grad);
    clbr=colorbar('Ticks',linspace(displacement(1),displacement(end),size(hfreq.R,1)), 'TickLabels',l);
    ylabel(clbr,'Displacement [mm]','FontSize',12)
    caxis([displacement(1),displacement(end)]) 
    fig=gcf;
    fig.Position=[400,540,820,400];
    set(gca,'FontSize',13);
    set(gca,'LineWidth',1.5);
    
    ax=gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - 1.75*ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
end

function displacement =getConvertedDisplacement(structInUse)
    DisplaceMentHelper =[1.96,1.94,1.8,1.71,1.6,1.49,1.39,1.29,1.2,1.1, ...
            0.98,0.89,0.79,0.7,0.62,0.52,0.42,0.33,0.23,0.09,0];
    [~,ind,~]=intersect(round(DisplaceMentHelper,3),round(structInUse.displacement,3));
    spacersToDisplacement = linspace(1.95,1.9,21)';
    displacement = round((spacersToDisplacement(flipud(ind))-structInUse.displacement),2);
    
end

