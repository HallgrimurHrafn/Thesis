% get estimated mechanical stiffness function
load('MechanicalStiffness.mat')
% uses mm, gives mm. this converts to meters.
k_a=@(x) 1000*ka(x*1000);
% get all parameters for a particular model. Scale model Joe Jensen.


fs=1e6; % hz, sampling freq
measurementLength=5; %sec, length of simulation
freq=[20,50,80,120,250]; %hz, driving frequency
V=linspace(20/1000,1,20); % driving voltage

% load xi and fL for particular model
name='10Hz 3RK';
loadname=strcat('FluxFuncs-',name,'.mat');
load(loadname)
ReducedKlippel=1;

%%

Ts=1/fs;
egs=cell(0,0);
timings=(0:fs*measurementLength)/fs;
for n=1:length(freq)
    for m=1:length(V)
        egs{m,n}=@(t) V(m)*sin(2*pi*t*freq(n));
    end
end

[F,G]=StateMatrix(k_a,f_L,xi,permRatio,sx,ReducedKlippel);

%% Run Simulation
parpool(24)
parfor n=1:length(egs)
    [i_current{n},x{n},u{n}]=StateSpaceSimulation(F, G, egs, Ts,timings);
end

i_current=reshape(i_current,length(V),length(freq));
x=reshape(x,length(V),length(freq));
u=reshape(u,length(V),length(freq));
egs=reshape(egs,length(V),length(freq));
Results=struct;
Results.curr=i_current;
Results.x=x;
Results.u=u;
Results.egs=egs;

save(strcat('Results-',name,'.mat'),'Results');
