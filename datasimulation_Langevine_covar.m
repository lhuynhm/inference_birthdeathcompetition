% Linh Huynh, April 2022
function [X,t,dt] = datasimulation_Langevine(bfunc,dfunc,nrun)
    tic   
    tmax = 3000; 
    dt   = 1/30; %logistic growth
    t    = 0:dt:tmax;

    X=nan(nrun,length(t)); % State vector preallocation
    dW=randn(size(X))*sqrt(dt); % Wiener process increments
    
    %X(:,1)=poissrnd(3); % logistic growth
    X(:,1) = 10;

    for j=2:length(t)
        covarX = -bfunc(X(:,j-1)).*dfunc(X(:,j-1)).*X(:,j-1)*dt;
        varX   = X(:,j-1).*(bfunc(X(:,j-1))+dfunc(X(:,j-1))) + covarX;
        X(:,j)=X(:,j-1)+...
            dt*X(:,j-1).*(bfunc(X(:,j-1))-dfunc(X(:,j-1)))+...
            sqrt(varX).*dW(:,j);
        X(:,j) = max(X(:,j),0);
    end

    % Ought to protect against X<0 but probably won't happen for these params.
    runtime=toc;
    disp('Run time')
    disp(runtime)

    %% Store
    %save datasimulation_Langevine.mat
end

