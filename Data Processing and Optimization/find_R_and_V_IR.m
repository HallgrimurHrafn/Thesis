% load singleLoopData
load('singleLoopDataAllNorm.mat')
c=SingleLoopData.current;
v=SingleLoopData.voltage;
dCurr=SingleLoopData.Dcurrent;

current=zeros(21,22);
voltage=zeros(21,22);
derivativeCurrent=zeros(21,22);
R=zeros(21,22);
v_c=zeros(21,22,8192);
vDC=zeros(21,22);
Vrms=zeros(21,22);

for m=1:21
    for n=1:22
%         measurement m, voltage/freq set n
%         cc is current, vv is voltage, ddc is derivative of current
        cc=squeeze(c(m,n,:));
        vv=squeeze(v(m,n,:));
        ddc=squeeze(dCurr(m,n,:));

%         making sure the particular voltage/freq combo exists
        if ~mean(cc)==0
%            find what the resistance should be such that V-RI oscilliates
%            around 0
            R(m,n)=mean(vv)/mean(cc);
%             V-IR
            v_c(m,n,:)=vv-R(m,n)*cc;
            vDC(m,n)=mean(vv);
            Vrms(m,n)=rms(vv-mean(vv));
        end
    end
end



save('R_v-ir-allNorm.mat','v_c','R','vDC','Vrms')