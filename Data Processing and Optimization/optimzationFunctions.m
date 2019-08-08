
function optimzationFunctions(order,sxskInit,saveName, allfreq,reducedKlippel,L00,T00,D,absXi)
    %% Load data
%     clear all
%     close all
    load('phaseAlignedInFreqSplitStructs-AllNorm.mat')

    %% Cost function
    % get cost and gradient functions for 10hz measurements
    [costFunction{1},gradFSx{1},gradFSkn{1}] ...
        =getCostAndGradientFunctions(h10,order,reducedKlippel,L00,T00,D,absXi);
    
    if allfreq
        % get cost and gradient functions for 100hz measurements
        [costFunction{2},gradFSx{2},gradFSkn{2}] ...
            =getCostAndGradientFunctions(h100,order,reducedKlippel,L00,T00,D,absXi);
        % 
        % % get cost and gradient functions for 1000hz measurements
        [costFunction{3},gradFSx{3},gradFSkn{3}] ...
            =getCostAndGradientFunctions(h1000,order,reducedKlippel,L00,T00,D,absXi);
    end

    %%
    % combine the costfunctions into 1 cell column
    costFunction=cat(1,costFunction{:});

    % combine the gradient functions into 1 cell columns
    gradFSx=cat(1,gradFSx{:});

    gradFSkn=cat(1,gradFSkn{:});


    % combine the gradient functions into a vector
    gradientFunctions = gradFSx;
    for o=1:order
        gradientFunctions=[gradientFunctions(:,:),gradFSkn(:,o)];
    end


    %% Gradient Descent
    epsilon=1e-6;
    maxNumSteps=20000;
    error=cell(1,16);
    x=cell(1,16);
    runs=zeros(8,1);
    
    m=20;
    fprintf('Progress:\n');
    fprintf(['\n' repmat('*',1,m) '\n\n']);
%     n=1;
    parfor n=1:m
        fprintf('\b|\n');
        [error{n},x{n},runs(n)]=GradientDescent(costFunction,gradientFunctions,[rand*3,rand(1,order)*6-3],epsilon,maxNumSteps);
    end
    fprintf('\n');
    gradDescResults=struct;
    gradDescResults.x=x;
    gradDescResults.error=error;
    gradDescResults.L00=L00;
    gradDescResults.T00=T00;
    gradDescResults.epsilon=epsilon;
    gradDescResults.maxNumSteps=maxNumSteps;
    save(saveName,'gradDescResults')
end



%% Gradient Descent!
function [error,x,n]=GradientDescent(costFunction,gradientFunctions,sxskInit, epsilon,maxNumSteps)


    % initial point:
    x(1,:)=sxskInit;
    
    %%
    n=1;
    epsilonCount=0;
    zn=zeros(1,length(sxskInit));
%     beta=0.15;
%     batchSize=ceil(linspace(150,200,maxNumSteps));
%     batchHelper=1:length(costFunction);
    error=zeros(maxNumSteps,1);
    error(1)=evaluateCostFunction(costFunction,x(1,:));
    improvement=1;
    while improvement>epsilon && n<maxNumSteps
    %     use back tracking for stepsize
%         Batch=datasample(batchHelper,333);
        gx=evaluateGradientFunctions(gradientFunctions,x(n,:));
        stepSize=getStepsize(costFunction,gx,x(n,:));
%         zn=gx'+beta*zn;
        x(n+1,:)=x(n,:)-stepSize*gx';
        error(n+1)=evaluateCostFunction(costFunction,x(n+1,:));
        n=n+1;
        improvement=error(n-1)-error(n);
%         if n>100
%             if error(n-1)-error(n)<epsilon
%                 epsilonCount=1;
% %             else
% %                 epsilonCount=0;
%             end
%         end
    end
    error=error(1:n);

end

function [costFunction,gradFSx,gradFSkn]=getCostAndGradientFunctions(hFreq,order,reducedKlippel,L00,T00,D,absXi)  
    %% important base functions
    % helper function to simplify writing the functions
    k=@(sk) 1:length(sk);
    % the Nonlinear flux-function with 4th order polynomial
    f_L=@(x,i_c,sx,sk) 1/(1-x^2/D^2+sum(sk.*(i_c+sx*x/D).^(2*k(sk))));
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

%     full klippel model part
    T_di = @(x,i_c,sx,sk) -T00*2*x*f_L(x,i_c,sx,sk)^2*xi(x,i_c,sx,sk);

    % gradient functions
    % L_edi function, derivative with respect to s_x
    L_edi_dsx = @(x,i_c,sx,sk) 2*L00*f_L(x,i_c,sx,sk)^2*x/D* ...
        (4*i_c*f_L(x,i_c,sx,sk)*xi(x,i_c,sx,sk)^2-xi(x,i_c,sx,sk)-i_c*eta(x,i_c,sx,sk));

    % L_edi function, derivative with respect to sn where n=[1,2,3,4] and by
    % inserting n=1, the function becomes L_edi derivative with respect to s_1
    % from the s_k series
    L_edi_dsk = @(x,i_c,sx,sk,n) L00*f_L(x,i_c,sx,sk)^2*(i_c+sx*x/D)^(2*n-1) * ...
        (4*i_c*f_L(x,i_c,sx,sk)*xi(x,i_c,sx,sk)*(i_c+sx*x/D)-(i_c+sx*x/D)-2*i_c*n);

    % gradient functions, T part. part of full klippel
    T_di_dsx = @(x,i_c,sx,sk) -T00*2*x^2/D*f_L(x,i_c,sx,sk)^2 * ...
        (eta(x,i_c,sx,sk)-4*f_L(x,i_c,sx,sk)*xi(x,i_c,sx,sk)^2);
    
    T_di_dsk = @(x,i_c,sx,sk,n) -T00*2*x*f_L(x,i_c,sx,sk)^2*(i_c+sx*x/D)^(2*n-1) * ...
        (n-2*f_L(x,i_c,sx,sk)*xi(x,i_c,sx,sk)*(i_c+sx*x/D));
    
    
    if reducedKlippel
        f = @(x,i_c,sx,sk) L_edi(x,i_c,sx,sk);
        f_dsx= @(x,i_c,sx,sk) L_edi_dsx(x,i_c,sx,sk);
        f_dsk = @(x,i_c,sx,sk,n) L_edi_dsk(x,i_c,sx,sk,n);
    else
        f = @(x,i_c,sx,sk) L_edi(x,i_c,sx,sk)+T_di(x,i_c,sx,sk);
        f_dsx= @(x,i_c,sx,sk) L_edi_dsx(x,i_c,sx,sk) + T_di_dsx(x,i_c,sx,sk);
        f_dsk = @(x,i_c,sx,sk,n) L_edi_dsk(x,i_c,sx,sk,n)+T_di_dsk(x,i_c,sx,sk,n);
        
    end

%% data point extraction

    current=hFreq.current;
    dcurrent=hFreq.dcurrent;
    v_ri=hFreq.v_riAligned;
%     displacement was noted in cm. Also as a distance from magnet, not
%     actual displacement. Not completely accurate, but close enough:
% create vector for all displacement used
    DisplaceMentHelper =[1.96,1.94,1.8,1.71,1.6,1.49,1.39,1.29,1.2,1.1, ...
        0.98,0.89,0.79,0.7,0.62,0.52,0.42,0.33,0.23,0.09,0];
%     find intersection, i.e. disregard missing displacement measurement if
%     relevant.
    [~,ind,~]=intersect(round(DisplaceMentHelper,3),round(hFreq.displacement,3));
%     conversion due from thickness of separators to armature location (arc
%     behavior effects what you can place around the reed)
    spacersToDisplacement = linspace(1.95,1.9,21)';
%     get the displacement in [m]
    displacement = (spacersToDisplacement(flipud(ind))-hFreq.displacement)/1000;
    
%     points=zeros(size(current,1),size(current,2),numPoints);
    displacement3dMatrix=zeros(size(current,1),size(current,2));
    
    maxVri=zeros(size(current,1),size(current,2));
    
    for m=1:size(current,1)
%         convert displacement vector to same shape as others
        displacement3dMatrix(m,:,:)=displacement(m);
        for n=1:size(current,2)
            vriInUse=squeeze(v_ri(m,n,:));
            [maxVri(m,n),maxInd]=max(vriInUse);
            cc=squeeze(current(m,n,:));
            measCurr(m,n)=rms(cc-mean(cc))+mean(cc);
            measDCurr(m,n)=rms(squeeze(dcurrent(m,n,:)));
            measV_ri(m,n)=rms(squeeze(v_ri(m,n,:)));
%             points=floor(mod(maxInd+(1:numPoints)*length(vriInUse)/numPoints,length(vriInUse)));
%             for p=1:length(points)
%                 measCurr(m,n,p)=current(m,n,points(p));
%                 measDCurr(m,n,p)=dcurrent(m,n,points(p));
%                 measV_ri(m,n,p)=v_ri(m,n,points(p));
%             end
            if maxVri(m,n)==0
%               preventing 0 division where no measurement is made later
%               on.
                maxVri(m,n)=1;
            end
            
        end
    end
%     for m=1:size(points,1)
%         for n=1:size(points,2)
%             for p=1:size(points,3)
%                 measCurr(m,n,p)=current(m,n,points(m,n,p));
%                 measDCurr(m,n,p)=dcurrent(m,n,points(m,n,p));
%                 measV_ri=v_ri(m,n,points(m,n,p));
%             end
%         end
%     end
%     measCurr=current(points);
%     measDCurr=dcurrent(points);
%     measV_ri=v_ri(points);


%     reshaping as a vector of meas*driveVolt*displacement,1,1
    measCurr=measCurr(:);
    measDCurr=measDCurr(:);
    measV_ri=measV_ri(:);
    displacementFlattened=displacement3dMatrix(:);
    WeightNormalized=measV_ri;
%     prevent zero division
    WeightNormalized(WeightNormalized==0)=1;

%% cost function and gradient functions generation
    costFunction=cell(length(measCurr),1);
    gradFSx=cell(length(measCurr),1);
    gradFSkn=cell(length(measCurr),order);
    
    for mnd=1:length(measCurr)
%         Le_di_atPoint function at this point
        f_atPoint=@(sx,sk) (measV_ri(mnd)- ... 
            f(displacementFlattened(mnd),measCurr(mnd),sx,sk)*measDCurr(mnd))/WeightNormalized(mnd);

%         calculating the cost function ((v-Ri -Le didt)/weight)^2
        costFunction{mnd}=@(sx,sk) f_atPoint(sx,sk)^2;
        
%         gradient function for sx: 2/weight*didt*Le_di*L_edi_dsx
        gradFSx{mnd}=@(sx,sk) -2/WeightNormalized(mnd)*f_atPoint(sx,sk)* ...
            f_dsx(displacementFlattened(mnd),measCurr(mnd),sx,sk)* ...
            measDCurr(mnd);
        
%         gradient functions for sk: 2/weight*didt*Le_di*L_edi_dsk
        for o=1:order
            gradFSkn{mnd,o}=@(sx,sk) -2/WeightNormalized(mnd)*f_atPoint(sx,sk)* ...
            f_dsk(displacementFlattened(mnd),measCurr(mnd),sx,sk,o)* ...
            measDCurr(mnd);
        end
        
    end
    
    
end

function result = evaluateCostFunction(costFunction, x)
    result=1/size(costFunction,1)*sum(cell2mat(cellfun(@(c) ...
        c(x(1),x(2:end)),costFunction,'UniformOutput',false)));
end

function Gr = evaluateGradientFunctions(gradientFunctions,x)
    [N,NumGrad]=size(gradientFunctions);
%         gr takes the shape of [0,0...] for [sx,sk1,sk2,sk3...skOrder]
    Gr=zeros(NumGrad,1);
%         for each derivative function in gradient
    for n=1:NumGrad
%             get gradient evaluation
        Gr(n)=1/N*sum(cell2mat(cellfun(@(c) ...
            c(x(1),x(2:NumGrad)),gradientFunctions(:,n),'UniformOutput',false)));
    end
end

function stepsize = getStepsize(costFunction, gx, x)
% beta value
    beta=0.5;
%     initial step size, gamma=t
    t=1;
%     get f(x)
    fx=evaluateCostFunction(costFunction,x);
%     get f(x-g(x))
    f_xgx = @(t) evaluateCostFunction(costFunction,x-t*gx');
%     get the euclidean norm squared of g(x)
    l2NormSquared_Gx=sum(gx.^2);
% if f_xgx>fx then f_xgx>fx-tgx for any positive t; Thus, reduce step size
% in f_xgx until f_xgx~>fx. then the rest of the algorithm can be run
%     while f_xgx>fx
%         t=beta*t;
%         f_xgx = evaluateCostFunction(costFunction,x-t*gx');
%         
%     end
%     while goldstein condition
    while f_xgx(t)>fx-t/2*l2NormSquared_Gx
        t=beta*t;
    end
    if t==1
        while f_xgx(t)<=fx-t/2*l2NormSquared_Gx
            t=t/beta;
        end
        t=t*beta;
    end
    stepsize=t;
end
