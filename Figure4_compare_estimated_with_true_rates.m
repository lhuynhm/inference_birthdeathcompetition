% Linh Huynh, April 2022
close all; clear all;
% Parameters
b0         = 1.1/120; % intrinsic birth rate
d0         = 0.1/120; % intrinsic death rate
r          = 1/120; % intrinsic net growth rate
K          = 1e5; % carrying capacity
binsz_orig = 1e3; % bin size
nrun       = 100; % number of trajectories

%% Case 1: gamma = 0
gamma = 0;
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n); 
dfunc = @(n) d0+(1-gamma)*(r/K).*n; 
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat';
tvec = tvec';
[~,nc] = size(Xmat);
for j = 1:nc
    fg1 = figure(1);
    plot(tvec,Xmat(:,j))
    set(gca,'FontSize',22)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg1, {'(A)'},'FontSize',25)
end
% Estimation of Birth and Death Rates
[brate_computed,drate_computed,~,~,~,N,~,~,~,~,~] = separatebirthdeathrates(Xmat,dt,binsz_orig);
N = (N(2:end)+N(1:end-1))./2;
brate_actual = bfunc(N).*N;
drate_actual = dfunc(N).*N;
% Plot
fg2 = figure(2);
plot(N,brate_computed,'k+','MarkerSize',15)
hold on
plot(N,brate_actual,'b-','LineWidth',1.5)
set(gca,'FontSize',22)
hold on
plot(N,drate_computed,'ko','MarkerSize',15)
hold on
plot(N,drate_actual,'r-','LineWidth',1.5)
xlabel('Number of Cells')
ylabel('Total Rates')
grid on
axis tight
AddLetters2Plots(fg2, {'(B)'},'FontSize',25)

%% Case 2: gamma = 0.5
gamma = 0.5;
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n);
dfunc = @(n) d0+(1-gamma)*(r/K).*n;
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat';
tvec = tvec';
[~,nc] = size(Xmat);
for j = 1:nc
    fg3 = figure(3);
    plot(tvec,Xmat(:,j))
    set(gca,'FontSize',22)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg3, {'(C)'},'FontSize',25)
end
% Estimation of Birth and Death Rates
[brate_computed,drate_computed,~,~,~,N,~,~,~,~,~] = separatebirthdeathrates(Xmat,dt,binsz_orig);
N = (N(2:end)+N(1:end-1))./2;
brate_actual = bfunc(N).*N;
drate_actual = dfunc(N).*N;
%Plot
fg4 = figure(4);
plot(N,brate_computed,'g+','MarkerSize',15)
hold on
plot(N,brate_actual,'b-','LineWidth',1.5)
set(gca,'FontSize',22)
hold on
plot(N,drate_computed,'go','MarkerSize',15)
hold on
plot(N,drate_actual,'r-','LineWidth',1.5)
xlabel('Number of Cells')
ylabel('Total Rates')
grid on
axis tight
AddLetters2Plots(fg4, {'(D)'},'FontSize',25)

%% Case 3: gamma = 1
gamma = 1;
% Functions
bfunc = @(n) max(0,b0-gamma*(r/K).*n);
dfunc = @(n) d0+(1-gamma)*(r/K).*n;
% Data Simulation
[Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
Xmat = Xmat';
tvec = tvec';
[nr,nc] = size(Xmat);
for j = 1:nc
    fg5 = figure(5);
    plot(tvec,Xmat(:,j))
    set(gca,'FontSize',22)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg5, {'(E)'},'FontSize',25)
end
% Estimation of Birth and Death Rates
[brate_computed,drate_computed,dNlengthvec,dNmeanvec,~,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);
N = (N(2:end)+N(1:end-1))./2;
brate_actual = bfunc(N).*N;
drate_actual = dfunc(N).*N;
% Plot
fg6 = figure(6);
plot(N,brate_computed,'m+','MarkerSize',15)
hold on
plot(N,brate_actual,'b-','LineWidth',1.5)
set(gca,'FontSize',22)
hold on
plot(N,drate_computed,'mo','MarkerSize',15)
hold on
plot(N,drate_actual,'r-','LineWidth',1.5)
xlabel('Number of Cells')
ylabel('Total Rates')
grid on
axis tight
AddLetters2Plots(fg6, {'(F)'},'FontSize',25)

        
