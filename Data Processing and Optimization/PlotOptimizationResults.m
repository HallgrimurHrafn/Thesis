loadNames={'GDResult_10hz_3th_FullKlippel','GDResult_10hz_3th_ReducedKlippel', ...
    'GDResult_10hz_6th_FullKlippel','GDResult_allfreq_3th_FullKlippel'};
Labels={'10FK3-', '10RK3-','10FK6-','AFK3-'};

%%
error=cell(5,1);
x=cell(5,1);
L00=zeros(5,1);
T00=zeros(5,1);
for m=1:4
    [error{m},x{m},L00(m),T00(m)]=minErrorGraph(loadNames{m});
end


%%
figure(1)
clf
hold on
for n=1:length(error)
    plot(error{n},'LineWidth',1.5);
end
hold off
set(gca,'XScale','log','YScale','log');
xlabel('Step Number [n]')
ylabel('Error As Given By Costfunction')
ylim([5e-2 1e4])
xlim([1 20000])
legend(Labels)
grid on
set(gca,'FontSize',13);
set(gca,'LineWidth',1.5);
fig=gcf;
fig.Position=[680   558   640   420];
ax=gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
saveas(gcf, 'Optimization Figures/SearchResults_error.png');



%% Handpicked Full Klippel
% Deff value
Ag = 72/10^6;
Am = 144/10^6;
Dm= Ag/Am*10/1000/1.12;
D=1.9/1000;
Deff=D+Dm;
plotBool=0;


% %% 10hz 3 Rk        3
% m=3;
% [xi,f_L,permRatio,sx]=plotReducedKlippel(x{m},L00(m),T00(m),Deff,1-mod(m,2),plotBool);
% save(strcat('SimPrep/FluxFuncs-10hz-3Rk.mat'),'xi','f_L','permRatio','sx')
% %% alhz 3 Fk        9
% m=9;
% [xi,f_L,permRatio,sx]=plotFullKlippel(x{m},L00(m),T00(m),Deff,1-mod(m,2),plotBool);
% save(strcat('SimPrep/FluxFuncs-allhz-3Fk.mat'),'xi','f_L','permRatio','sx')
% %% 10hz 3 Fk        1
% m=1;
% [xi,f_L,permRatio,sx]=plotFullKlippel(x{m},L00(m),T00(m),Deff,1-mod(m,2),plotBool);
% save(strcat('SimPrep/FluxFuncs-10hz-3Fk.mat'),'xi','f_L','permRatio','sx')
% %% 10hz 6 Fk        5
% m=5;
% [xi,f_L,permRatio,sx]=plotFullKlippel(x{m},L00(m),T00(m),Deff,1-mod(m,2),plotBool);
% save(strcat('SimPrep/FluxFuncs-10hz-6Fk.mat'),'xi','f_L','permRatio','sx')
% %% 10hz 3 Fk Xi     2
% m=2;
% [xi,f_L,permRatio,sx]=plotFullKlippel(x{m},L00(m),T00(m),Deff,1-mod(m,2),plotBool);
% save(strcat('SimPrep/FluxFuncs-10hz-absXi-3Fk.mat'),'xi','f_L','permRatio','sx')


%% save nonlinear values
% Deff value
Ag = 72/10^6;
Am = 144/10^6;
Dm= Ag/Am*10/1000/1.12;
D=19/1000;
Deff=D+Dm;
plotBool=1;
absXi=[0,0,0,1,0];
RK=[0,1,0,0,0];
for n=1:4
    [xi,f_L,permRatio,sx]=plotModel(x{n},L00(n),T00(n),Deff,absXi(n),plotBool,RK(n),Labels{n});
    save(strcat('SimPrep/FluxFuncs-',Labels{n},'.mat'),'xi','f_L','permRatio','sx')
end


%% functions
function [error,x,L00,T00]=minErrorGraph(loadName)
    load(strcat('Results/',loadName,'.mat'));
    error=gradDescResults.error;
    x=gradDescResults.x;
    minlist=zeros(1,length(error));
    for m=1:length(error)
        minlist(m)=min(error{m});
    end
    [~,n]=sort(minlist);
    sx=-1;
    m=1;
%     while sx<0
%         m=m+1;
%         xtemp=x{n(m)};
%         [~,minErr]=min(error{n(m)});
%         sx=xtemp(minErr,1);
%     end
    x=x{n(m)};
    [~,minErr]=min(error{n(m)});
    x=x(minErr,:);
    error=error{n(m)};
    L00=gradDescResults.L00;
    T00=gradDescResults.T00;
end




function [xi,f_L,permRatio,sx]=plotModel(x,L00,T00,Deff,absXi,plotBool,RK,labels)
    D=1.85/1000;
    sx=x(1);
    sk=x(2:end);
    [f_L,xi,L_edi,T_di,f,p_a,R_a,permRatio]= nonLinearFunctions(L00,T00,Deff,absXi);
    f=@(x,i) f(x,i,sx,sk);
    f_L=@(x,i) f_L(x,i,sx,sk);
    xi=@(x,i) xi(x,i,sx,sk);
    L_edi=@(x,i) L_edi(x,i,sx,sk);
    T_di=@(x,i) -T_di(x,i,sx,sk);
    p_a=@(x,i) p_a(x,i,sx,sk);
    R_a=@(x,i) R_a(x,i,sx,sk);
    permRatio=@(x,i) permRatio(x,i,sx,sk);
    if RK
        f=@(x,i) L_edi(x,i);
    end
    g=[-D,D, -0.12,0.12];
    if plotBool
        plotMeshSpecifics(f_L,g,'fL','Amplitude',labels)
        plotMeshSpecifics(f,g,'Ind','Inductance [H]',labels)
        plotMeshSpecifics(R_a,g,'Ra','Magnetic Reluctance [H^{-1}]',labels)
        plotMeshSpecifics(xi,g,'xi','Amplitude ',labels)
    end

end

function plotMeshSpecifics(f,grid,titleName,zlbl,labels)
    name=strcat('Optimization Figures/',labels,titleName);
    figure
    s=fmesh(f, grid);
    s.FaceColor='interp';
    s.FaceAlpha=0.8;
    xlabel('Displacement [m]')
    ylabel('Current [A]')
    zlabel(zlbl)
    view(28,19)
    saveas(gcf,strcat(name,'.png'))
    view(0,90)
    colorbar
    saveas(gcf,strcat(name,'-flat','.png'))
%     title(titleName)

end



function [f_L,xi,L_edi,T_di,f,p_a,R_a,permRatio]= nonLinearFunctions(L00,T00,D,absXi)
    Ag = 72/10^6;
    u0=4*pi*1e-7;
    k=@(sk) 1:length(sk);
    % the Nonlinear flux-function with 4th order polynomial
    f_L=@(x,i_c,sx,sk) 1/(1-x^2/D^2+sum(sk.*(i_c+sx*x/D).^(2*k(sk))));
    
%     2u0A/Dpa
    permRatio=@(x,i_c,sx,sk) sum(sk.*(i_c+sx*x/D).^(2*k(sk)));
    
    % p_a, magnetic permeability of armature, R_a reluctance
    p_a=@(x,i_c,sx,sk) 2*u0*Ag/(D*sum(sk.*(i_c+sx*x/D).^(2*k(sk))));
    R_a=@(x,i_c,sx,sk) (D*sum(sk.*(i_c+sx*x/D).^(2*k(sk))))/(2*u0*Ag);
    
    % a term from the derivative of the nonlinear flux function with 4th order
    % polynomial
    if absXi
        xi =@(x,i_c,sx,sk) sum(k(sk).*sk.*(abs(i_c+sx*x/D)).^(2*k(sk)-1));
    else
        xi =@(x,i_c,sx,sk) sum(k(sk).*sk.*(i_c+sx*x/D).^(2*k(sk)-1));
    end
    
    % derivative of xi, eta (greek n) , (excludes internal derivatives)
    eta = @(x,i_c,sx,sk) sum(k(sk).*(2*k(sk)-1).*sk.*(i_c+sx*x/D).^(2*k(sk)-2));

    % L_edi function from klippels method: klippel reduced
    L_edi = @(x,i_c,sx,sk) L00*f_L(x,i_c,sx,sk)*(1-2*i_c*f_L(x,i_c,sx,sk)*xi(x,i_c,sx,sk));

    % full klippel model part
    T_di = @(x,i_c,sx,sk) T00*2*x*f_L(x,i_c,sx,sk)^2*xi(x,i_c,sx,sk);

    f = @(x,i_c,sx,sk) L_edi(x,i_c,sx,sk)-T_di(x,i_c,sx,sk);
end