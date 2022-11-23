% Linh Huynh, September 2022
close all; clear all;
% Parameters
b0         = 2/120; % intrinsic birth rate
d0         = 0.1/120; % intrinsic death rate
r          = b0-d0; % intrinsic net growth rate
K          = 1e5; % carrying capacity
binsz_orig = 1e3; % bin size
nrun       = 100; % number of trajectories
% Functions
bfunc = @(n) max(b0,3*b0-4*(b0/K).*n);
dfunc = @(n) max(4*(b0/K).*n-3*b0,(8/14)*(b0/K).*n);
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
    %AddLetters2Plots(fg1, {'(A)'},'FontSize',25)
end
% Estimation of Birth and Death Rates
[brate_computed,drate_computed,~,~,~,N,~,dNvec,binvec_new,binvec] = separatebirthdeathrates(Xmat,dt,binsz_orig);
N = (N(2:end)+N(1:end-1))./2;
brate_actual = bfunc(N).*N;
drate_actual = dfunc(N).*N;
% Plot
fg2 = figure(2);
plot(N,brate_computed,'k+','MarkerSize',15)
hold on
plot(N,brate_actual,'b-','LineWidth',1.5)
set(gca,'FontSize',15)
hold on
plot(N,drate_computed,'ko','MarkerSize',15)
hold on
plot(N,drate_actual,'r-.','LineWidth',1.5)
xlabel('Number of Cells')
ylabel('Total Rates')
grid on
axis tight
legend('inferred birth', 'true birth', 'inferred death', 'true death')
%AddLetters2Plots(fg2, {'(B)'},'FontSize',25)
