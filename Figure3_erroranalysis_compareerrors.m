% Linh Huynh, April 2022
close all; clear all;
%note: dNlengthvec = vector of sample sizes
power_vec = [1.5;2;2.5;3;3.5;4;4.5];
binvec    = 10.^power_vec; 
K         = 1e5; 

%% Case 1: gamma = 0
err_birth_expect_binvec = [];
err_birth_var_binvec    = [];
err_birth_std_binvec    = [];
err_birth_sample_binvec = [];
err_death_expect_binvec = [];
err_death_var_binvec    = [];
err_death_std_binvec    = [];
err_death_sample_binvec = [];
gamma = 0;
nrun  = 100; %number of trajectories
b0    = 2/120; %instrinsic per capita birth rate
d0    = 0.1/120; %instrinsic per capita death rate
r     = b0-d0; %instrinsic per capita net growth rate 

% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
for i = 1:length(binvec)
    binsz_orig = binvec(i); %bin size
    
    % Estimation of Birth and Death Rates
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

    % True Birth and Death Rates
    N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;

    % Different averages of Delta N (Nk,Nki,empirical)
    pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
    pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
    sample_meanNki =  dNmeanvec; %sample mean of DeltaNki

    % Different variances of Delta N (Nk,Nki,empirical)
    pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
    pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                    +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                    -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                    +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                    -(pop_meanNki.^2);
    sample_varNki = dNvarvec;

    % Errors in estimating birth rate
    err_birth_expect        = ((pop_varNk-pop_varNki)+(pop_meanNk-pop_meanNki))./(2*dt);
    err_birth_expect_binvec = [err_birth_expect_binvec; norm(err_birth_expect,2)];
    err_birth_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_birth_var_binvec    = [err_birth_var_binvec; norm(err_birth_var,2)];
    err_birth_std           = sqrt(err_birth_var);
    err_birth_std_binvec    = [err_birth_std_binvec; norm(err_birth_std,2)];
    err_birth_sample        = brate_actual-brate_computed;
    err_birth_sample_binvec = [err_birth_sample_binvec; norm(err_birth_sample,2)];
    
    % Errors in estimating birth rate
    err_death_expect        = ((pop_varNk-pop_varNki)-(pop_meanNk-pop_meanNki))./(2*dt);
    err_death_expect_binvec = [err_death_expect_binvec; norm(err_death_expect,2)];
    err_death_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_death_var_binvec    = [err_death_var_binvec; norm(err_death_var,2)];
    err_death_std           = sqrt(err_death_var);
    err_death_std_binvec    = [err_death_std_binvec; norm(err_death_std,2)];
    err_death_sample        = drate_actual-drate_computed;
    err_death_sample_binvec = [err_death_sample_binvec; norm(err_death_sample,2)];
end
fg1 = figure(1);
plot(log10(binvec),err_birth_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_birth_std_binvec,'^-','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_birth_sample_binvec,'ko-','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Birth Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg1, {'(A)'},'FontSize',25)

fg2 = figure(2);
plot(log10(binvec),err_death_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_death_std_binvec,'^-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_death_sample_binvec,'ko-','LineWidth',1.5,'MarkerSize',10)
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Death Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg2, {'(B)'},'FontSize',25)

%% Case 2: gamma = 0.5
err_birth_expect_binvec = [];
err_birth_var_binvec    = [];
err_birth_std_binvec    = [];
err_birth_sample_binvec = [];
err_death_expect_binvec = [];
err_death_var_binvec    = [];
err_death_std_binvec    = [];
err_death_sample_binvec = [];
gamma = 0.5;
nrun  = 100; %number of trajectories
b0    = 2/120; %instrinsic per capita birth rate
d0    = 0.1/120; %instrinsic per capita death rate
r     = b0-d0; %instrinsic per capita net growth rate 

% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
for i = 1:length(binvec)
    binsz_orig = binvec(i); %bin size
    
    % Estimation of Birth and Death Rates
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

    % True Birth and Death Rates
    N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;

    % Different averages of Delta N (Nk,Nki,empirical)
    pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
    pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
    sample_meanNki =  dNmeanvec; %sample mean of DeltaNki

    % Different variances of Delta N (Nk,Nki,empirical)
    pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
    pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                    +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                    -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                    +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                    -(pop_meanNki.^2);
    sample_varNki = dNvarvec;

    % Errors in estimating birth rate
    err_birth_expect        = ((pop_varNk-pop_varNki)+(pop_meanNk-pop_meanNki))./(2*dt);
    err_birth_expect_binvec = [err_birth_expect_binvec; norm(err_birth_expect,2)];
    err_birth_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_birth_var_binvec    = [err_birth_var_binvec; norm(err_birth_var,2)];
    err_birth_std           = sqrt(err_birth_var);
    err_birth_std_binvec    = [err_birth_std_binvec; norm(err_birth_std,2)];
    err_birth_sample        = brate_actual-brate_computed;
    err_birth_sample_binvec = [err_birth_sample_binvec; norm(err_birth_sample,2)];
    
    % Errors in estimating birth rate
    err_death_expect        = ((pop_varNk-pop_varNki)-(pop_meanNk-pop_meanNki))./(2*dt);
    err_death_expect_binvec = [err_death_expect_binvec; norm(err_death_expect,2)];
    err_death_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_death_var_binvec    = [err_death_var_binvec; norm(err_death_var,2)];
    err_death_std           = sqrt(err_birth_var);
    err_death_std_binvec    = [err_death_std_binvec; norm(err_death_std,2)];
    err_death_sample        = drate_actual-drate_computed;
    err_death_sample_binvec = [err_death_sample_binvec; norm(err_death_sample,2)];
end
fg3 = figure(3);
plot(log10(binvec),err_birth_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_birth_std_binvec,'^-','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_birth_sample_binvec,'go-','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','g')
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Birth Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg3, {'(C)'},'FontSize',25)

fg4 = figure(4);
plot(log10(binvec),err_death_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_death_std_binvec,'^-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_death_sample_binvec,'go-','LineWidth',1.5,'MarkerSize',10)
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Death Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg4, {'(D)'},'FontSize',25)

%% Case 3: gamma = 1
err_birth_expect_binvec = [];
err_birth_var_binvec    = [];
err_birth_std_binvec    = [];
err_birth_sample_binvec = [];
err_death_expect_binvec = [];
err_death_var_binvec    = [];
err_death_std_binvec    = [];
err_death_sample_binvec = [];
gamma = 1;
nrun  = 100; %number of trajectories
b0    = 2/120; %instrinsic per capita birth rate
d0    = 0.1/120; %instrinsic per capita death rate
r     = b0-d0; %instrinsic per capita net growth rate 

% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
for i = 1:length(binvec)
    binsz_orig = binvec(i); %bin size
    
    % Estimation of Birth and Death Rates
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

    % True Birth and Death Rates
    N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;

    % Different averages of Delta N (Nk,Nki,empirical)
    pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
    pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
    sample_meanNki =  dNmeanvec; %sample mean of DeltaNki

    % Different variances of Delta N (Nk,Nki,empirical)
    pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
    pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                    +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                    -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                    +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                    -(pop_meanNki.^2);
    sample_varNki = dNvarvec;

    % Errors in estimating birth rate
    err_birth_expect        = ((pop_varNk-pop_varNki)+(pop_meanNk-pop_meanNki))./(2*dt);
    err_birth_expect_binvec = [err_birth_expect_binvec; norm(err_birth_expect,2)];
    err_birth_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_birth_var_binvec    = [err_birth_var_binvec; norm(err_birth_var,2)];
    err_birth_std           = sqrt(err_birth_var);
    err_birth_std_binvec    = [err_birth_std_binvec; norm(err_birth_std,2)];
    err_birth_sample        = brate_actual-brate_computed;
    err_birth_sample_binvec = [err_birth_sample_binvec; norm(err_birth_sample,2)];
    
    % Errors in estimating birth rate
    err_death_expect        = ((pop_varNk-pop_varNki)-(pop_meanNk-pop_meanNki))./(2*dt);
    err_death_expect_binvec = [err_death_expect_binvec; norm(err_death_expect,2)];
    err_death_var           = ((pop_varNki./dNlengthvec)+((2*pop_varNki.^2)./(dNlengthvec-1)))./(4*dt^2);
    err_death_var_binvec    = [err_death_var_binvec; norm(err_death_var,2)];
    err_death_std           = sqrt(err_death_var);
    err_death_std_binvec    = [err_death_std_binvec; norm(err_death_std,2)];
    err_death_sample        = drate_actual-drate_computed;
    err_death_sample_binvec = [err_death_sample_binvec; norm(err_death_sample,2)];
end
fg5 = figure(5);
plot(log10(binvec),err_birth_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_birth_std_binvec,'^-','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_birth_sample_binvec,'mo-','LineWidth',1.5,'MarkerSize',7,'MarkerFaceColor','m','MarkerEdgeColor','m')
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Birth Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg5, {'(E)'},'FontSize',25)

fg6 = figure(6);
plot(log10(binvec),err_death_expect_binvec,'s-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','b')
set(gca,'FontSize',20)
hold on
plot(log10(binvec),err_death_std_binvec,'^-','LineWidth',1.5,'MarkerSize',10,'MarkerEdgeColor','r')
hold on
plot(log10(binvec),err_death_sample_binvec,'mo-','LineWidth',1.5,'MarkerSize',10)
xlabel('log_{10}(Bin Size)')
ylabel('2-Norm Death Rate Error','Interpreter','tex')
legend('expected','std','empirical')
AddLetters2Plots(fg6, {'(F)'},'FontSize',25)


