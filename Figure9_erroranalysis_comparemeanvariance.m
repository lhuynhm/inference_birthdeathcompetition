% Linh Huynh, April 2022
close all; clear all;
%% Case 1: gamma = 0
gamma      = 0;
binsz_orig = 1e3; %bin size
nrun       = 100; %number of trajectories
b0         = 2/120; %instrinsic per capita birth rate
d0         = 0.1/120; %instrinsic per capita death rate
r          = b0-d0; %instrinsic per capita net growth rate 
K          = 1e5; %carrying capacity
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
%[nr,nc] = size(Xmat); %numbers of rows and columns of cell number ensemble

% Estimation of Birth and Death Rates
[brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

% True Birth and Death Rates
N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
% brate_actual = bfunc(N).*N;
% drate_actual = dfunc(N).*N;

% Comparison between different averages of Delta N (Nk,Nki,empirical)
pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
sample_meanNki =  dNmeanvec; %sample mean of DeltaNki
%plot
fg1 = figure(1);
plot(N,pop_meanNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_meanNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(N,sample_meanNki,'ko','MarkerSize',10)
xlabel('Number of Cells')
ylabel('Mean of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg1, {'(A)'},'FontSize',25)

%Comparison between different variances of Delta N (Nk,Nki,empirical)
pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                -(pop_meanNki.^2);
sample_varNki = dNvarvec;
%plot
fg2 = figure(2);
plot(N,pop_varNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_varNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(N,sample_varNki,'ko','MarkerSize',10)
xlabel('Number of Cells')
ylabel('Variance of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg2, {'(B)'},'FontSize',25)

%% Case 2: gamma = 0.5
gamma      = 0.5;
binsz_orig = 1e3; %bin size
nrun       = 100; %number of trajectories
b0         = 2/120; %instrinsic per capita birth rate
d0         = 0.1/120; %instrinsic per capita death rate
r          = b0-d0; %instrinsic per capita net growth rate 
K          = 1e5; %carrying capacity
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
% [nr,nc] = size(Xmat); %numbers of rows and columns of cell number ensemble

% Estimation of Birth and Death Rates
[brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

% True Birth and Death Rates
N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
% brate_actual = bfunc(N).*N;
% drate_actual = dfunc(N).*N;

% Comparison between different averages of Delta N (Nk,Nki,empirical)
pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
sample_meanNki =  dNmeanvec; %sample mean of DeltaNki
%plot
fg3 = figure(3);
plot(N,pop_meanNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_meanNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(N,sample_meanNki,'go','MarkerSize',10)
xlabel('Number of Cells')
ylabel('Mean of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg3, {'(C)'},'FontSize',25)

%Comparison between different variances of Delta N (Nk,Nki,empirical)
pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                -(pop_meanNki.^2);
sample_varNki = dNvarvec;
%plot
fg4 = figure(4);
plot(N,pop_varNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_varNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
xlabel('Number of Cells')
hold on
plot(N,sample_varNki,'go','MarkerSize',10)
ylabel('Variance of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg4, {'(D)'},'FontSize',25)

%% Case 3: gamma = 1
gamma      = 1;
binsz_orig = 1e3; %bin size
nrun       = 100; %number of trajectories
b0         = 2/120; %instrinsic per capita birth rate
d0         = 0.1/120; %instrinsic per capita death rate
r          = b0-d0; %instrinsic per capita net growth rate 
K          = 1e5; %carrying capacity
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); %analytical per capita birth rate
dfunc = @(n) d0+(1-gamma)*(r/K).*n; %analytical per capita death rate
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat'; %the whole ensemble of cell numbers
tvec = tvec'; %matrix of time vectors corresponding to the cell number trajectories
% [nr,nc] = size(Xmat); %numbers of rows and columns of cell number ensemble

% Estimation of Birth and Death Rates
[brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);

% True Birth and Death Rates
N = (N(2:end)+N(1:end-1))./2; %midpoints of the bins
% brate_actual = bfunc(N).*N;
% drate_actual = dfunc(N).*N;

% Comparison between different averages of Delta N (Nk,Nki,empirical)
pop_meanNk     = (bfunc(N)-dfunc(N)).*N*dt; %population/analytical mean of DeltaNk
pop_meanNki    = pop_meanNk + (b0-d0-2*(r/K).*N)*dt*(binsz_orig/2) - (r/K)*dt*(binsz_orig^2)*(1/3); %population/analytical mean of DeltaNki
sample_meanNki =  dNmeanvec; %sample mean of DeltaNki
%plot
fg5 = figure(5);
plot(N,pop_meanNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_meanNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
plot(N,sample_meanNki,'mo','MarkerSize',10)
xlabel('Number of Cells')
ylabel('Mean of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg5, {'(E)'},'FontSize',25)

%Comparison between different variances of Delta N (Nk,Nki,empirical)
pop_varNk     = (bfunc(N)+dfunc(N)).*N*dt;
pop_varNki    = pop_varNk + (b0+d0+2*(1-2*gamma)*(r/K).*N)*dt*(binsz_orig/2) + (1-2*gamma)*(r/K)*dt*(binsz_orig^2)*(1/3)...
                +(dt^2)*(r^2)*(1/(3*binsz_orig))*((N+binsz_orig).^3 - N.^3)...
                -(2/K)*(dt^2)*(r^2)*(1/(4*binsz_orig))*((N+binsz_orig).^4 - N.^4)...
                +(1/(K^2))*(dt^2)*(r^2)*(1/(5*binsz_orig))*((N+binsz_orig).^5 - N.^5)...
                -(pop_meanNki.^2);
sample_varNki = dNvarvec;
%plot
fg6 = figure(6);
plot(N,pop_varNk,'r','LineWidth',5)
set(gca,'FontSize',20)
hold on
plot(N,pop_varNki,'s','MarkerSize',4,'MarkerFaceColor','b','MarkerEdgeColor','b')
xlabel('Number of Cells')
hold on
plot(N,sample_varNki,'mo','MarkerSize',10)
ylabel('Variance of \DeltaN','Interpreter','tex')
%title(['\gamma=',num2str(gamma)])
legend('theoryNk','theoryNki','empirNki')
AddLetters2Plots(fg6, {'(F)'},'FontSize',25)


