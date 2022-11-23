
%% Case 1: gamma = 0
close all; clear all;
% Parameters
gamma_vec = [0; 0; 0.5]; 
b0_vec    = [1.1*(1/120); 1.1*(1/120); 1.1*(1/120)];
d0_vec    = [0.1*(1/120); 0.1*(1/120); 0.1*(1/120)];
K_vec     = [1e5; 0.5*1e5; 0.5*1e5];
binsz_orig = 1e3;

for i = 1:3
    switch i 
        case 1
            color = 'k';
        case 2
            color = 'r';
        case 3
            color = 'b';
    end
    gamma = gamma_vec(i);
    b0    = b0_vec(i);
    d0    = d0_vec(i);
    r     = b0 - d0;
    K     = K_vec(i);
    bfunc = @(n) max(0,b0-gamma*(r/K).*n);
    dfunc = @(n) d0+(1-gamma)*(r/K).*n;
    nrun  = 100;
    [Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    Xmat = Xmat';
    tvec = tvec';
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,~,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);
    brate_computed = max(brate_computed,0);
    Nendpoints = N;
    N = (N(2:end)+N(1:end-1))./2; %midpoints
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;
    [nrXmat,ncXmat] = size(Xmat);
    trajrand = Xmat(:,randi([1,ncXmat]));
    [binindex_vec,endpoints] = discretize(trajrand,Nendpoints);
    brate_time = brate_computed(binindex_vec);
    drate_time = drate_computed(binindex_vec);
    
    fg1 = figure(1);
    plot(N,brate_computed,[color,'+'],'MarkerSize',10)
    set(gca,'FontSize',23)
    hold on
    plot(N,drate_computed,[color,'o'],'MarkerSize',10)
    hold on
    xlabel('Number of Cells')
    ylabel('Estimated Total Rates')
    grid on
    AddLetters2Plots(fg1, {'(A)'},'FontSize',27)
    
    fg2 = figure(2);
    plot(tvec,trajrand,color,'LineWidth',2)
    set(gca,'FontSize',23)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg2, {'(D)'},'FontSize',27)
    
    fg3 = figure(3);
    plot(tvec,brate_time,[color,'+'],'MarkerSize',5)
    hold on
    plot(tvec,drate_time,[color,'o'],'MarkerSize',5)
    set(gca,'FontSize',23)
    xlabel('Time (Arbitrary Units)')
    ylabel('Estimated Total Rates')
    grid on
    AddLetters2Plots(fg3, {'(G)'},'FontSize',27)
end

%% Case 2: gamma = 0.5
% Parameters
gamma_vec = [0.5; 0.25; 0.75]; 
b0_vec    = [1.1*(1/120); 1.1*(1/120); 1.1*(1/120)];
d0_vec    = [0.1*(1/120); 0.1*(1/120); 0.1*(1/120)];
K_vec     = [1e5; 0.5*1e5; 0.5*1e5];
binsz_orig = 1e3;

for i = 1:3
    switch i 
        case 1
            color = 'g';
        case 2
            color = 'r';
        case 3
            color = 'b';
    end
    gamma = gamma_vec(i);
    b0    = b0_vec(i);
    d0    = d0_vec(i);
    r     = b0-d0;
    K     = K_vec(i);
    bfunc = @(n) max(0,b0-gamma*(r/K).*n);
    dfunc = @(n) d0+(1-gamma)*(r/K).*n;
    nrun  = 100;
    [Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    Xmat = Xmat';
    tvec = tvec';
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,~,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);
    brate_computed = max(brate_computed,0);
    Nendpoints = N;
    N = (N(2:end)+N(1:end-1))./2; %midpoints
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;
    [nrXmat,ncXmat] = size(Xmat);
    trajrand = Xmat(:,randi([1,ncXmat]));
    N_bin = [N; N(end)+ N(2)-N(1)];
    [binindex_vec,endpoints] = discretize(trajrand,Nendpoints);
    brate_time = brate_computed(binindex_vec);
    drate_time = drate_computed(binindex_vec);
    
    fg4 = figure(4);
    plot(N,brate_computed,[color,'+'],'MarkerSize',15)
    set(gca,'FontSize',23)
    hold on
    plot(N,drate_computed,[color,'o'],'MarkerSize',15)
    hold on
    xlabel('Number of Cells')
    ylabel('Estimated Total Rates')
    grid on
    AddLetters2Plots(fg4, {'(B)'},'FontSize',27)
    
    fg5 = figure(5);
    plot(tvec,trajrand,color,'LineWidth',2)
    set(gca,'FontSize',23)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg5, {'(E)'},'FontSize',27)
    
    fg6 = figure(6);
    plot(tvec,brate_time,[color,'+'],'MarkerSize',5)
    hold on
    plot(tvec,drate_time,[color,'o'],'MarkerSize',5)
    set(gca,'FontSize',23)
    xlabel('Time (Arbitrary Units)')
    ylabel('Estimated Total Rates')
    grid on 
    AddLetters2Plots(fg6, {'(H)'},'FontSize',27)
end
%% Case 3: gamma = 1
% Parameters
gamma_vec = [1; 0.5; 1]; 
b0_vec    = [1.1*(1/120); 1.1*(1/120); 1.1*(1/120)];
d0_vec    = [0.1*(1/120); 0.1*(1/120); 0.1*(1/120)];
K_vec     = [1e5; 0.5*1e5; 0.5*1e5];
binsz_orig = 1e3;

for i = 1:3
    switch i 
        case 1
            color = 'm';
        case 2
            color = 'r';
        case 3
            color = 'b';
    end
    gamma = gamma_vec(i);
    b0    = b0_vec(i);
    d0    = d0_vec(i);
    r     = b0-d0;
    K     = K_vec(i);
    bfunc = @(n) max(0,b0-gamma*(r/K).*n);
    dfunc = @(n) d0+(1-gamma)*(r/K).*n;
    nrun  = 100;
    [Xmat,tvec,dt] = datasimulation_Langevine(bfunc,dfunc,nrun);
    Xmat = Xmat';
    tvec = tvec';
    [brate_computed,drate_computed,dNlengthvec,dNmeanvec,~,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig);
    brate_computed = max(brate_computed,0);
    Nendpoints = N;
    N = (N(2:end)+N(1:end-1))./2; %midpoints
    brate_actual = bfunc(N).*N;
    drate_actual = dfunc(N).*N;
    [nrXmat,ncXmat] = size(Xmat);
    trajrand = Xmat(:,randi([1,ncXmat]));
    N_bin = [N; N(end)+ N(2)-N(1)];
    [binindex_vec,endpoints] = discretize(trajrand,Nendpoints);
    brate_time = brate_computed(binindex_vec);
    drate_time = drate_computed(binindex_vec);
    
    fg7 = figure(7);
    plot(N,brate_computed,[color,'+'],'MarkerSize',15)
    set(gca,'FontSize',23)
    hold on
    plot(N,drate_computed,[color,'o'],'MarkerSize',15)
    hold on
    xlabel('Number of Cells')
    ylabel('Estimated Total Rates')
    grid on
    AddLetters2Plots(fg7, {'(C)'},'FontSize',27)
    
    fg8 = figure(8);
    plot(tvec,trajrand,color,'LineWidth',2)
    set(gca,'FontSize',23)
    hold on
    xlabel('Time (Arbitrary Units)')
    ylabel('Number of Cells')
    grid on
    AddLetters2Plots(fg8, {'(F)'},'FontSize',27)
    
    fg9 = figure(9);
    plot(tvec,brate_time,[color,'+'],'MarkerSize',5)
    hold on
    plot(tvec,drate_time,[color,'o'],'MarkerSize',5)
    set(gca,'FontSize',23)
    xlabel('Time (Arbitrary Units)')
    ylabel('Estimated Total Rates')
    grid on  
    AddLetters2Plots(fg9, {'(I)'},'FontSize',27)
end
