% wBar=waitbar(0,'starting up');

% air gap+ equivalent distance of internal magnetic resistance.
Ag = 72/10^6;
Am = 144/10^6;
Dm= Ag/Am*10/1000/1.12;
D=1.9/1000;
Deff=D+Dm;
u0 = 4*pi*1e-7;
ur=1.12;
N = 1755;
Br = 0.41; % 0.4-0.42 T
% Fm=Br*Dm/u0; -> Br=Fm*u0/Dm;    Tm/(TM/A) = A.. correct
Fm =Br*Dm/(u0*ur) ;
L00 = 2*u0*Ag*N^2/Deff; 
T00 = 2*u0*Ag*Fm*N/Deff^2;
    
absXiStr={'','AbsXi_'};
f={'10hz_','allfreq_'};
rk={'th_FullKlippel','th_ReducedKlippel'};
settings=dec2bin(0:15)-'0';
baseString='Results/GDResult_';
for n= [1,2,3,5,9] %[1,2,3,5,9]
    n
    allfreq=settings(n,1);
    order=3+3*settings(n,2);
    reducedKlippel=settings(n,3);
    absXi=settings(n,4);
    saveName=strcat(baseString,absXiStr{absXi+1},f{1+allfreq},string(order),rk{1+reducedKlippel}) ;                  
    sxskInit=[0.5,0.1*ones(1,order)];
    optimzationFunctions(order,sxskInit,saveName,allfreq,reducedKlippel,L00,T00,Deff,absXi);
%     waitbar(0.125*n,wBar,string(n));
end

