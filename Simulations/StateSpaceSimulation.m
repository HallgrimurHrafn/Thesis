function [i,x,u]=StateSpaceSimulation(F, G, egs, Ts, L)
% follows the order [i,x,u]. i.e. states(1,n) is the estimated current at a
% given step n.
% F(x,i)*[i,x,u] & G(x,i)*egs
    f=@(state,t) F(state(2),state(1))*state+G(state(2),state(1))*egs(t); 
    states=zeros(3,length(L));
    dL=L(2)-L(1);
    for n=1:length(L)
        k1=Ts*f(states(:,n),L(n));
        k2=Ts*f(states(:,n)+k1/2,L(n)+dL/2);
        k3=Ts*f(states(:,n)+k2/2,L(n)+dL/2);
        k4=Ts*f(states(:,n)+k3,L(n)+dL);

        states(:,n+1)=states(:,n)+(k1+2*k2+2*k3+k4)/6;
    end
%     
%     for n=1:length(L)
%         k=states(:,n)+Ts*f(states(:,n),L(n));
%         states(:,n+1)=states(:,n)+Ts*f(states(:,n)+k/2,L(n)+dL/2);
%     end
    
%     translate to curr,disp and velocity, removing every other sample as
%     rk4 uses higher resolution
%     for n=1:length(L)
%         states(:,n+1)=states(:,n)+Ts*f(states(:,n),L(n));
%     end
% 
%         k1=Ts*(F(states(2,n),states(1,n))*states(:,n) ...
%             + G(states(2,n),states(1,n))*egs(L(n)));
%         
%         k2=Ts*(F(states(2,n)+k1(2)/2,states(1,n)+k1(1)/2)*states(:,n) ...
%             + G(states(2,n)+k1(2)/2,states(1,n)+k1(1)/2)*egs(L(n)+dn/2));
%         
%         k3=Ts*(F(states(2,n)+k2(2)/2,states(1,n)+k2(1)/2)*states(:,n) ...
%             + G(states(2,n)+k2(2)/2,states(1,n)+k2(1)/2)*egs(L(n)+dn/2));
%         
%         k4=Ts*(F(states(2,n)+k3(2),states(1,n)+k3(1))*states(:,n) ...
%             + G(states(2,n)+k3(2),states(1,n)+k3(1))*egs(L(n)+dn));

    i=states(1,:);
    x=states(2,:);
    u=states(3,:);
end