%% Load data
load('frequencySplitStruct-AllNorm.mat')

%% Phase Alignment

% perform for 10hz struct
dcurr=h10.dcurrent;
v_ri=h10.v_ri;
driveVolt=h10.driveVolt;
[v_riAligned,phaseDiff10]=alignPhases(v_ri,dcurr,driveVolt);
h10.v_riAligned=v_riAligned;
h10.phaseDiff=phaseDiff10;

%% perform for 100hz struct
dcurr=h100.dcurrent;
v_ri=h100.v_ri;
driveVolt=h100.driveVolt;
[v_riAligned,phaseDiff100]=alignPhases(v_ri,dcurr,driveVolt);
h100.v_riAligned=v_riAligned;
h100.phaseDiff=phaseDiff100;

%% perform for 1000hz struct
dcurr=h1000.dcurrent;
v_ri=h1000.v_ri;
driveVolt=h1000.driveVolt;
[v_riAligned,phaseDiff1000]=alignPhases(v_ri,dcurr,driveVolt);
h1000.v_riAligned=v_riAligned;
h1000.phaseDiff=phaseDiff1000;

%%
save('phaseAlignedInFreqSplitStructs-AllNorm','h10','h100','h1000')





function [v_riAlignedMatrix, phaseDiff]=alignPhases(v_ri,dcurr, driveVolt)


%     Initialize storage for the aligned vector
    v_riAlignedMatrix=zeros(size(v_ri));
    phaseDiff=zeros(size(dcurr,1),size(dcurr,2));
%     for every measurement
    for m=1:size(dcurr,1)
%         for every point in v_ri and dcurr
        for n=1:size(dcurr,2)
            if driveVolt(m,n)~=0
%             get the shifted version of v_ri for aligned phases, get the
%             phase difference as well.
                [v_riAlignedMatrix(m,n,:),phaseDiff(m,n)]=fundamentalPhaseAlignment(squeeze(v_ri(m,n,:)),squeeze(dcurr(m,n,:)));
            else
                phaseDiff(m,n)=NaN;
                v_riAlignedMatrix(m,n,:)=squeeze(v_ri(m,n,:));
            end
        end
    end
end


function [v_riAligned, phaseDiff]=fundamentalPhaseAlignment(v_ri,dcurr)

    v_riLong=repmat(v_ri,[10,1]);
    v_riNorm=v_riLong/max(v_riLong);
    ff_v_ri=fft(v_riNorm);
    ff_v_ri([1:10,12:end-10,end-8:end])=0;
    v_riFundaMental=real(ifft(ff_v_ri));

    dccNorm=dcurr/max(dcurr);
    dccNorm=repmat(dccNorm,[10,1]);
    ff_dcurr=fft(dccNorm);
    ff_dcurr([1:10,12:end-10,end-8:end])=0;
    dcurr_fundamental=real(ifft(ff_dcurr));

    v_riFundaMental=v_riFundaMental(1:8192*2);
    func2minimize=@(phaseDiff) sum(abs(dcurr_fundamental(1:8192)-v_riFundaMental(1+phaseDiff:phaseDiff+8192)));
    [~,phaseDiff]=min(arrayfun(func2minimize,0:8191));
    phaseDiff=mod(phaseDiff,8192);
    
    v_riAligned=v_riLong(phaseDiff+1:phaseDiff+8192);

end