% Linh Huynh, April 2022
function [brate_computed,drate_computed,dNlengthvec,dNmeanvec,dNvarvec,N,dt_method,CIbrupvec,CIbrlovec,CIdrupvec,CIdrlovec] = separatebirthdeathrates(Xmat,dt,binsz_orig)
    dt_data         = dt;
    dt_method       = dt_data;
    step            = floor(dt_method/dt_data);
    Nmat            = Xmat; %N(t)
    Nplmat          = Xmat((1+step):end,:); %N(t+dt)
    [nrNpl,ncNpl]   = size(Nplmat);
    [nrX,ncX]       = size(Xmat);
    rowdiffer       = nrX - nrNpl;
    Nendvec         = Xmat(end,:);
    Nappendmat      = repmat(Nendvec,[rowdiffer,1]);
    Nplmat          = [Nplmat;Nappendmat];
    Nvec            = reshape(Nmat,[nrX*ncX,1]);
    Nplvec          = reshape(Nplmat,[nrX*ncX,1]);
    [Nvecsrt,idsrt] = sort(Nvec);
    Nplvecsrt       = Nplvec(idsrt);
    Xmax            = max(Nvec);
    Xmin            = min(Nvec);
    %ninterval      = floor((Xmax-Xmin)/binsz_orig);
    ninterval       = ceil((Xmax-Xmin)/binsz_orig);
    %binsz_real    = (Xmax-Xmin)/ninterval;
    Xmax            = Xmin + ninterval*binsz_orig;
    edges           = (linspace(Xmin,Xmax,ninterval))';
    [binvec]        = discretize(Nvecsrt,edges);
    bin_unique = unique(binvec);
    maxbin  = max(binvec);
    minbin  = min(binvec);
    binnumb = (minbin:maxbin)';
    freq    = hist(binvec',binnumb');
    result  = [binnumb freq'];
    freq    = (freq(find(freq~=0)))';
    binvec_new = repelem((1:length(freq))',freq);
    save('binvec_new.mat','binvec_new')
    save('Nvecsrt.mat','Nvecsrt')
    dNvec       = Nplvecsrt - Nvecsrt;
    dNlengthvec = splitapply(@length,dNvec,binvec_new);
    dNmeanvec   = splitapply(@mean,dNvec,binvec_new);
    dNvarvec     = splitapply(@var,dNvec,binvec_new);
    N = edges;
    brate_computed = (dNvarvec+dNmeanvec)./(2*dt_method);
    drate_computed = (dNvarvec-dNmeanvec)./(2*dt_method);  
    freq_nz     = (freq(find(freq~=0)))'; %'nz' means no zeros
    degfreevec  = freq_nz-1;
    tupvec      = tinv(0.95,degfreevec);%t-Score
    SEMvec      = sqrt(dNvarvec)./sqrt(freq_nz); % Standard Error
    CImupvec    = dNmeanvec + tupvec.*SEMvec; 
    CImlovec    = dNmeanvec - tupvec.*SEMvec;
    chi2lo      = chi2inv(0.025,degfreevec); %chi square inverse
    chi2up      = chi2inv(0.975,degfreevec);
    CIvarupvec  = (degfreevec.*dNvarvec)./chi2lo;
    CIvarlovec  = (degfreevec.*dNvarvec)./chi2up;
    CIbrupvec  = (CIvarupvec+CImupvec)./(2*dt_method);
    CIbrlovec  = (CIvarlovec+CImlovec)./(2*dt_method);
    CIdrupvec  = (CIvarupvec-CImupvec)./(2*dt_method);
    CIdrlovec  = (CIvarlovec-CImlovec)./(2*dt_method);
  
end

