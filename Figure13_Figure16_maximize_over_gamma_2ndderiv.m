% Linh Huynh, April 2022
close all; clear all;
deltagamma = 1e-3; %to approximate second derivative
%% Case gamma = 0
gamma_data1 = 0;
gamma_root_vec1 = [];
count_rootmax = 0;
% test and data parameters
K   = 1e5;
r   = 1/120;
b0  = 2/120;
d0  = 1/120;
% data
bfunc = @(n) max(b0-gamma_data1*(r/K).*n,0);
dfunc = @(n) d0+(1-gamma_data1)*(r/K).*n;
n_sample = 1e2;
n_skip = 0;
d2f_vec1 = [];
d2fexact_vec1 = [];
df_vec1 = [];
for id_run = 1:n_sample
    nrun  = 1;
    [Xvec,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    X_vec = Xvec';
    [lenrX,lencX] = size(X_vec);
    DeltaX_vec = X_vec(2:end,1) - X_vec(1:end-1,1);
    X_vec  = X_vec(1:end-1,1);
    Deltat = dt;   
    % functions
    bxj_nonneg = @(gamma) b0 - gamma*(r/K).*X_vec;
    bxj = @(gamma) max(bxj_nonneg(gamma),0);
    dxj = @(gamma) d0 + (1-gamma)*(r/K).*X_vec;
    % derivative of (bxj + dxj)
    d_bxj_p_dxj = @(gamma) -1.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma);
    % derivative of (bxj - dxj)
    d_bxj_m_dxj = @(gamma) 0.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma); 
    df  = @(gamma) sum((-0.5./(bxj(gamma) + dxj(gamma))).*d_bxj_p_dxj(gamma))...
                   + sum(((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat)./(bxj(gamma) + dxj(gamma))).*d_bxj_m_dxj(gamma))...
                   + sum(d_bxj_p_dxj(gamma).*(0.5/Deltat).*((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat).^2)./(X_vec.*(bxj(gamma) + dxj(gamma)).^2));          
    % root of df contrained on [0,1]
    gamma_initial = [-2 2];
    gamma_root = fzero(df,gamma_initial);
    gamma_root_vec1 = [gamma_root_vec1; gamma_root];
    d2f = (df(gamma_root + deltagamma) - df(gamma_root - deltagamma))/(2*deltagamma); %second derivative
    d2f_vec1 = [d2f_vec1; d2f];
    meandx = @ (gamma) X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat;
    vardx  = @ (gamma) X_vec.*(bxj(gamma)+dxj(gamma)).*Deltat;
    c4     = Deltat*(r/K).*X_vec.^2;
    c2     = Deltat*(r/K).*X_vec.^2;
    indx   = find(bxj(gamma_root) > 0);
    c4(indx) = 2*Deltat*(r/K).*X_vec.^2;
    c2(indx) = zeros(length(indx),1);
    d2fexact = @(gamma) sum((c4./(vardx(gamma)).^2).*(0.5*c4+2*c2.*(DeltaX_vec - meandx(gamma))-((DeltaX_vec - meandx(gamma)).^2).*c4./vardx(gamma)));
    d2fexact_root = d2fexact(gamma_root);
    d2fexact_vec1 = [d2fexact_vec1; d2fexact_root];
    df_root = df(gamma_root);
    df_vec1 = [df_vec1; df_root];
end

fg1 = figure(1);
hist(gamma_root_vec1)
set(gca,'FontSize',20)
xlim([0-0.1,1+0.1])
ylim([0,30])
xlabel('Estimated Optimal \gamma Value')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(gamma_root_vec)),', Var=',num2str(var(gamma_root_vec)),', Median=',num2str(median(gamma_root_vec))])
AddLetters2Plots(fg1, {'(A)'},'FontSize',25)

error_vec1 = gamma_data1 - gamma_root_vec1;
fg2 = figure(2);
hist(error_vec1)
set(gca,'FontSize',20)
xlim([-0.05,0.05])
ylim([0,30])
xlabel('\gamma Estimation Error')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(error_vec)),', Var=',num2str(var(error_vec)),', Median=',num2str(median(error_vec))])
AddLetters2Plots(fg2, {'(B)'},'FontSize',25)

fg7 = figure(7);
hist(d2fexact_vec1)
set(gca,'FontSize',18)
ylim([0,30])
xlabel('Second Derivative at Root')
ylabel('Numerical Frequency')
AddLetters2Plots(fg7, {'(A)'},'FontSize',25)

%% Case gamma = 0.5
gamma_data2 = 0.5;
gamma_root_vec2 = [];
count_rootmax = 0;
% test and data parameters
K   = 1e5;
r   = 1/120;
b0  = 2/120;
d0  = 1/120;
% data
bfunc = @(n) max(b0-gamma_data2*(r/K).*n,0);
dfunc = @(n) d0+(1-gamma_data2)*(r/K).*n;
n_sample = 1e2;
n_skip = 0;
d2f_vec2 = [];
d2fexact_vec2 = [];
for id_run = 1:n_sample
    nrun  = 1;
    [Xvec,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    X_vec = Xvec';
    [lenrX,lencX] = size(X_vec);
    DeltaX_vec = X_vec(2:end,1) - X_vec(1:end-1,1);
    X_vec  = X_vec(1:end-1,1);
    Deltat = dt;   
    % functions
    bxj_nonneg = @(gamma) b0 - gamma*(r/K).*X_vec;
    bxj = @(gamma) max(bxj_nonneg(gamma),0);
    dxj = @(gamma) d0 + (1-gamma)*(r/K).*X_vec;
    % derivative of (bxj + dxj)
    d_bxj_p_dxj = @(gamma) -1.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma);
    % derivative of (bxj - dxj)
    d_bxj_m_dxj = @(gamma) 0.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma); 
    df  = @(gamma) sum((-0.5./(bxj(gamma) + dxj(gamma))).*d_bxj_p_dxj(gamma))...
                   + sum(((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat)./(bxj(gamma) + dxj(gamma))).*d_bxj_m_dxj(gamma))...
                   + sum(d_bxj_p_dxj(gamma).*(0.5/Deltat).*((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat).^2)./(X_vec.*(bxj(gamma) + dxj(gamma)).^2));               
    % root of df contrained on [0,1]
    gamma_initial = [-2 2];
    gamma_root = fzero(df,gamma_initial);
    gamma_root_vec2 = [gamma_root_vec2; gamma_root];
    d2f = (df(gamma_root + deltagamma) - df(gamma_root))/deltagamma; %second derivative
    d2f_vec2 = [d2f_vec2; d2f];
    meandx = @ (gamma) X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat;
    vardx  = @ (gamma) X_vec.*(bxj(gamma)+dxj(gamma)).*Deltat;
    c4     = Deltat*(r/K).*X_vec.^2;
    c2     = Deltat*(r/K).*X_vec.^2;
    indx   = find(bxj(gamma_root) > 0);
    c4(indx) = 2*Deltat*(r/K).*X_vec.^2;
    c2(indx) = zeros(length(indx),1);
    d2fexact = @(gamma) sum((c4./(vardx(gamma)).^2).*(0.5*c4+2*c2.*(DeltaX_vec - meandx(gamma))-((DeltaX_vec - meandx(gamma)).^2).*c4./vardx(gamma)));
    d2fexact_root = d2fexact(gamma_root);
    d2fexact_vec2 = [d2fexact_vec2; d2fexact_root];
end

fg3 = figure(3);
hist(gamma_root_vec2)
set(gca,'FontSize',20)
xlim([0-0.1,1+0.1])
ylim([0,30])
xlabel('Estimated Optimal \gamma Value')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(gamma_root_vec)),', Var=',num2str(var(gamma_root_vec)),', Median=',num2str(median(gamma_root_vec))])
AddLetters2Plots(fg3, {'(C)'},'FontSize',25)

error_vec2 = gamma_data2 - gamma_root_vec2;
fg4 = figure(4);
hist(error_vec2)
set(gca,'FontSize',20)
xlim([-0.05,0.05])
ylim([0,30])
xlabel('\gamma Estimation Error')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(error_vec)),', Var=',num2str(var(error_vec)),', Median=',num2str(median(error_vec))])
AddLetters2Plots(fg4, {'(D)'},'FontSize',25)

fg8 = figure(8);
hist(d2fexact_vec2)
set(gca,'FontSize',18)
ylim([0,30])
xlabel('Second Derivative at Root')
ylabel('Numerical Frequency')
AddLetters2Plots(fg8, {'(B)'},'FontSize',25)

%% Case gamma = 1
gamma_data3 = 1;
gamma_root_vec3 = [];
count_rootmax = 0;
% test and data parameters
K   = 1e2;
r   = 1/120;
b0  = 2/120;
d0  = 1/120;
% data
bfunc = @(n) max(b0-gamma_data3*(r/K).*n,0);
dfunc = @(n) d0+(1-gamma_data3)*(r/K).*n;
n_sample = 1e2;
n_skip = 0;
d2f_vec3 = [];
d2fexact_vec3 = [];
for id_run = 1:n_sample
    nrun  = 1;
    [Xvec,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    X_vec = Xvec';
    [lenrX,lencX] = size(X_vec);
    DeltaX_vec = X_vec(2:end,1) - X_vec(1:end-1,1);
    X_vec  = X_vec(1:end-1,1);
    Deltat = dt;   
    % functions
    bxj_nonneg = @(gamma) b0 - gamma*(r/K).*X_vec;
    bxj = @(gamma) max(bxj_nonneg(gamma),0);
    dxj = @(gamma) d0 + (1-gamma)*(r/K).*X_vec;
    % derivative of (bxj + dxj)
    d_bxj_p_dxj = @(gamma) -1.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma);
    % derivative of (bxj - dxj)
    d_bxj_m_dxj = @(gamma) 0.5*(r/K).*X_vec...
                  - 0.5*(r/K).*X_vec.*abs(bxj_nonneg(gamma))./bxj_nonneg(gamma); 
    df  = @(gamma) sum((-0.5./(bxj(gamma) + dxj(gamma))).*d_bxj_p_dxj(gamma))...
                   + sum(((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat)./(bxj(gamma) + dxj(gamma))).*d_bxj_m_dxj(gamma))...
                   + sum(d_bxj_p_dxj(gamma).*(0.5/Deltat).*((DeltaX_vec - X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat).^2)./(X_vec.*(bxj(gamma) + dxj(gamma)).^2));          

    % root of df contrained on [0,1]
    gamma_initial = [-2 2];
    gamma_root = fzero(df,gamma_initial);
    gamma_root_vec3 = [gamma_root_vec3; gamma_root];
    d2f = (df(gamma_root + deltagamma) - df(gamma_root))/deltagamma; %second derivative
    d2f_vec3 = [d2f_vec3; d2f];
    meandx = @ (gamma) X_vec.*(bxj(gamma)-dxj(gamma)).*Deltat;
    vardx  = @ (gamma) X_vec.*(bxj(gamma)+dxj(gamma)).*Deltat;
    c4     = Deltat*(r/K).*X_vec.^2;
    c2     = Deltat*(r/K).*X_vec.^2;
    indx   = find(bxj(gamma_root) > 0);
    c4(indx) = 2*Deltat*(r/K).*X_vec.^2;
    c2(indx) = zeros(length(indx),1);
    d2fexact = @(gamma) sum((c4./(vardx(gamma)).^2).*(0.5*c4+2*c2.*(DeltaX_vec - meandx(gamma))-((DeltaX_vec - meandx(gamma)).^2).*c4./vardx(gamma)));
    d2fexact_root = d2fexact(gamma_root);
    d2fexact_vec3 = [d2fexact_vec3; d2fexact_root];
end

fg5 = figure(5);
hist(gamma_root_vec3)
set(gca,'FontSize',20)
xlim([0-0.1,1+0.1])
ylim([0,30])
xlabel('Estimated Optimal \gamma Value')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(gamma_root_vec)),', Var=',num2str(var(gamma_root_vec)),', Median=',num2str(median(gamma_root_vec))])
AddLetters2Plots(fg5, {'(E)'},'FontSize',25)

error_vec3 = gamma_data3 - gamma_root_vec3;
fg6 = figure(6);
hist(error_vec3)
set(gca,'FontSize',20)
xlim([-0.05,0.05])
ylim([0,30])
xlabel('\gamma Estimation Error')
ylabel('Numerical Frequency')
%title(['Mean=',num2str(mean(error_vec)),', Var=',num2str(var(error_vec)),', Median=',num2str(median(error_vec))])
AddLetters2Plots(fg6, {'(F)'},'FontSize',25)

fg9 = figure(9);
hist(d2fexact_vec3)
set(gca,'FontSize',18)
ylim([0,30])
xlabel('Second Derivative at Root')
ylabel('Numerical Frequency')
AddLetters2Plots(fg9, {'(C)'},'FontSize',25)

