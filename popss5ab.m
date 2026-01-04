clear all
close all

%calculate fractional information loss for (locally optimal) 1D decoders of 2D
%population compared to 2D decoder.

% Select stimuli
nd = 16;
ns = 16;

%consider randomizing stimulus selection on a much finer grid and running 
%many more reps

dvals = [1:nd]/nd*2*pi;%[1:nd]/nd*2*pi; %radians
svals = linspace(1,7,ns);  %log2 speed

%priors on stim distributions
pdvals=ones(nd,1)/nd;  %conservative prior on dir vals
psvals=ones(ns,1)/ns;
pjoint = psvals*pdvals';
% note: choice of priors will impact the info loss from 1D decoders 

%Twin = 0.15;
%Twin = .15;

Twin_vec = 0.02:.02:.200;

rho = 0.4;   %mean peak correlation value  %loop on this eventually
rho_base = 0.05;
rmax = 120; %mean value
rbase = 5; %base
%bw_dir_mean = deg2rad(40); %+5deg
bw_dir_mean = deg2rad(35);
bw_spd_mean = 2; %log2 units %+.5 log2 units
 
nCells_vec = 2.^(4:11); 

%nCells_vec = 128;

nReps_vec = [40*ones(1,length(find(nCells_vec<=1024))) 10*ones(1,length(find(nCells_vec>1024)))];



%differential correlations:  
diffcorron=1;  %information limiting correlations on or off
 
behFIs=100;
behFId=1/(3/180*pi)^2;
 
maxFIs=2*behFIs;    % depends on s  or rather .1*2.^s so computed for each stim
maxFId=2*behFId;  %dir threshold fixed

gensamps=0; 

for nidx = 1:length(nCells_vec)
    ['nidx = ' num2str(nidx) ' of ' num2str(length(nCells_vec))]
    tic
    nn = nCells_vec(nidx);
    nReps = nReps_vec(nidx);
       
    for tidx = 1:length(Twin_vec)
       ['tidx = ' num2str(tidx) ' of ' num2str(length(Twin_vec))]
       Twin = Twin_vec(tidx);
         

       spriorvar = 4*sum(svals.^2.*psvals')-sum(svals.*psvals')^2;
       dpriorvar = sum(dvals.^2.*pdvals')-sum(dvals.*pdvals')^2;
       priorFIs = 1./spriorvar;
       priorFId = 1./dpriorvar;
       priorFI = diag([priorFIs,priorFId]);
       
       ENTs = 1/2*log2(2*pi*exp(1)*spriorvar);
       ENTd = 1/2*log2(2*pi*exp(1)*dpriorvar);


        for ridx = 1:nReps
            ['ridx = ' num2str(ridx) ' of ' num2str(nReps)] 
             %set up population, get rates for stims 

             [pop,Sigma,Q,overlap] = MTpoiss_rates2D(nn,dvals,svals,Twin,rho,rho_base,rmax,rbase,bw_dir_mean,bw_spd_mean);
             %note nPDirs, nPSpds = nCells
             mu = pop.mu;
             dmuds = pop.dmuds;
             dmudd = pop.dmudd;
             olap.dir(nidx,ridx) = mean(mean(overlap.dir));
             olap.spd(nidx,ridx) = mean(mean(overlap.spd));
        
             mubar=zeros(size(mu,1),1);
             Sigmabar=zeros(size(mu,1));
             
             %add differential correlations
              if(diffcorron)
                  ILCorr = 0.0*ones(ns,nd);
                 for s=1:ns
                     for d=1:nd
                         epsilon=diag(sqrt([1/maxFIs,1/maxFId]))*[1,ILCorr(s,d);ILCorr(s,d),1;]/(1-ILCorr(s,d)^2)*diag(sqrt([1/maxFIs,1/maxFId]));
%                          epsilon=epsilon/(Twin/max(Twin_vec));
                         Sigma(:,:,s,d) = Sigma(:,:,s,d) + [dmuds(:,s,d),dmudd(:,s,d)]*epsilon*[dmuds(:,s,d),dmudd(:,s,d)]';
                         Sigma(:,:,s,d) = (Sigma(:,:,s,d) + Sigma(:,:,s,d)')/2;
                         
                     end
                 end
              else
                  epsilon=diag([1/maxFIs,1/maxFId]);
              end
              
             mubar=zeros(size(mu,1),1);
             Sigmabar=zeros(size(mu,1));
             %construct 1D estimates of FI ignoring other stim dimension
             Sigmas=zeros(nn,nn,ns);
             Sigmad=zeros(nn,nn,nd);
             mus=zeros(nn,ns);
             dmusds=zeros(nn,ns);
             mud=zeros(nn,nd);
             dmuddd=zeros(nn,nd);
             Sigmans=zeros(nn,1);
             Sigmand=zeros(nn,1);
             smean=svals*psvals;
             dmean=dvals*pdvals;
             for s=1:ns
             for d=1:nd
                 mubar(:,1) = mubar(:,1) + mu(:,s,d)*psvals(s)*pdvals(d); 
             end
             end
             
             
             
             for s=1:ns
                 for d=1:nd
                     Sigmas(:,:,s) = Sigmas(:,:,s) + (Sigma(:,:,s,d) + mu(:,s,d)*mu(:,s,d)')*pdvals(d);
                     Sigmas(:,:,s) =  (Sigmas(:,:,s) +  Sigmas(:,:,s)')/2;
                     Sigmad(:,:,d) = Sigmad(:,:,d) + (Sigma(:,:,s,d) + mu(:,s,d)*mu(:,s,d)')*psvals(s);
                     Sigmad(:,:,d) =  (Sigmad(:,:,d) +  Sigmad(:,:,d)')/2;
                     mus(:,s) = mus(:,s) + mu(:,s,d)*pdvals(d);
                     mud(:,d) = mud(:,d) + mu(:,s,d)*psvals(s);
                     
                     Sigmans(:,1)=Sigmans(:,1) + (mu(:,s,d)-mubar)*(svals(s)-smean)*psvals(s)*pdvals(d);
                     Sigmand(:,1)=Sigmand(:,1) + (mu(:,s,d)-mubar)*(dvals(d)-dmean)*pdvals(s)*psvals(d);
                     
                     dmusds(:,s) = dmusds(:,s) + dmuds(:,s,d)*pdvals(d);
                     dmuddd(:,d) = dmuddd(:,d) + dmudd(:,s,d)*psvals(s);
                     
                       Sigmabar = Sigmabar + ...
                              (squeeze(Sigma(:,:,s,d)) + squeeze(mu(:,s,d))*squeeze(mu(:,s,d))')*psvals(s)*pdvals(d);
                 end
             end

             
             for s=1:ns
                 Sigmas(:,:,s) = Sigmas(:,:,s) - mus(:,s)*mus(:,s)';
             end
             for d=1:nd
                 Sigmad(:,:,d) = Sigmad(:,:,d) - mud(:,d)*mud(:,d)';
             end
             
             Sigmabar = Sigmabar - mubar*mubar';
             Sigmabar = (Sigmabar + Sigmabar')/2;
             %Compute total fisher information

             
             % At this point we have computed all the relevant quantities
             % needed to generate the parameters of the various decoders
             % and their biases and variances.  In the next section we will
             % compute the mean and variance of each decoder as a function
             % of the presented values of the stimuls as well as the
             % derivative of the mean with respect to s and d (for the
             % purposes of computing Fisher Information)
             
             % Global Decoder
             GL.ws = inv(Sigmabar)*Sigmans; %nn x 1
             GL.wd = inv(Sigmabar)*Sigmand;
             GL.w=[GL.ws,GL.wd]'; % 2 x nn
             
             GL.est.MI = ENTs + ENTd;
             GL.est.MIdiag = ENTs + ENTd;
             GL.s.mean=zeros(ns,1);
             GL.s.var=zeros(ns,1);
             GL.s.gradmean=zeros(ns,1);
             GL.d.mean=zeros(nd,1);
             GL.d.gradmean=zeros(nd,1);
             GL.d.var=zeros(nd,1);
             
             for s=1:ns
             for d=1:nd                 
                 GL.est.mean(:,s,d) = GL.w*(mu(:,s,d)-mubar) + [smean;dmean];
                 GL.est.bias(:,s,d) = GL.est.mean(:,s,d)-[svals(s);dvals(d)];
                 GL.est.gradmean(:,:,s,d) = [GL.w*dmuds(:,s,d),GL.w*dmudd(:,s,d)];
                 GL.est.var(:,:,s,d) = GL.w*Sigma(:,:,s,d)*GL.w';
                 GL.est.err(:,:,s,d) = GL.est.var(:,:,s,d) + (GL.est.bias(:,s,d)*GL.est.bias(:,s,d)');
                 GL.est.deterr(s,d) = det(GL.est.err(:,:,s,d));                 

                 GL.est.corr(s,d) = GL.est.var(1,2,s,d)/sqrt(GL.est.var(1,1,s,d)*GL.est.var(2,2,s,d));
                 
                 V=squeeze(GL.est.var(:,:,s,d));
                 B=squeeze(GL.est.gradmean(:,:,s,d));          
                 GL.est.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                 GL.est.invFI(:,:,s,d) = inv(GL.est.FI(:,:,s,d));  % variance and covariance of 
                                                                   % local bias removed estimates
                 GL.est.FIdiag(:,s,d)=diag(B).^2./diag(GL.est.var(:,:,s,d)) + diag(priorFI);         
                 GL.est.detFI(s,d) = det(GL.est.FI(:,:,s,d));
                 GL.est.detFIdiag(s,d) = prod(GL.est.FIdiag(:,s,d));
                 
                 GL.est.MI = GL.est.MI - 1/2*log2((2*pi*exp(1))^2./GL.est.detFI(s,d))*pjoint(s,d);
                 GL.est.MIdiag = GL.est.MIdiag - 1/2*log2((2*pi*exp(1))^2./GL.est.detFIdiag(s,d))*pjoint(s,d);
                     
                 GL.s.mean(s) = GL.s.mean(s) + GL.est.mean(1,s,d)*pdvals(d);
                 GL.s.gradmean(s) = GL.s.gradmean(s) + GL.est.gradmean(1,1,s,d)*pdvals(d);
                 GL.s.var(s) = GL.s.var(s) + (GL.est.var(1,1,s,d) + GL.est.mean(1,s,d).^2)*pdvals(d);
                 
                 GL.d.mean(d) = GL.d.mean(d) + GL.est.mean(2,s,d)*psvals(s);
                 GL.d.gradmean(d) = GL.d.gradmean(d) + GL.est.gradmean(2,2,s,d)*psvals(s);
                 GL.d.var(d) = GL.d.var(d) + (GL.est.var(2,2,s,d) + GL.est.mean(2,s,d).^2)*psvals(s);
             end
             end
             
             GL.est.MIerr = - 1/2*sum(sum(log2((2*pi*exp(1))^2.*GL.est.deterr).*pjoint)) + ENTs + ENTd;
             GL.s.var=GL.s.var-GL.s.mean.^2;
             GL.s.bias=GL.s.mean - svals';
             GL.s.err = GL.s.var + GL.s.bias.^2;
             GL.s.FI = GL.s.gradmean.^2./GL.s.var + priorFIs;
             GL.s.invFI = 1./GL.s.FI;
             GL.s.MI = -1/2*sum(log2(2*pi*exp(1)./GL.s.FI).*psvals) + ENTs;
             GL.s.MIerr = -1/2*sum(log2(2*pi*exp(1).*GL.s.err).*psvals) + ENTs;
             
             GL.d.var=GL.d.var-GL.d.mean.^2;
             GL.d.bias=GL.d.mean - dvals';
             GL.d.err = GL.d.var + GL.d.bias.^2;
             GL.d.FI = GL.d.gradmean.^2./GL.d.var + priorFId;
             GL.d.invFI = 1./GL.d.FI;
             GL.d.MI = -1/2*sum(log2(2*pi*exp(1)./GL.d.FI).*pdvals) + ENTd;
             GL.d.MIerr = -1/2*sum(log2(2*pi*exp(1).*GL.d.err).*pdvals) + ENTd;
             GL.MI = GL.est.MI;
             GL.Synergy = GL.est.MI - GL.s.MI - GL.d.MI;
             GL.StimSpecSynergy = -1/2*log2(1-GL.est.corr.^2);
             GL.meanSSS = sum(sum(GL.StimSpecSynergy.*pjoint));
             GL.ErrSynergy = GL.est.MIerr - GL.s.MIerr - GL.d.MIerr;

             % Population Vector Decoder (requires some simulation)
             
             PV.ws = [pop.PSpds;ones(1,nn);]; %nn x 2 %2xnn
             PV.wd = [cos(pop.PDirs);sin(pop.PDirs);]; %nn x 2 %2xnn
             
             PV.w = [PV.ws;PV.wd;];
             PV.MI= ENTs + ENTd;
             
            
             for s=1:ns
             for d=1:nd
                 PV.mu(:,s,d) = PV.w*mu(:,s,d);
                 PV.var(:,:,s,d) = PV.w*Sigma(:,:,s,d)*PV.w' +eps;
                 PV.gradmean(:,:,s,d) =  [PV.w*dmuds(:,s,d),PV.w*dmudd(:,s,d)];
                 
                 V=squeeze(PV.var(:,:,s,d))+eps;
                 B=squeeze(PV.gradmean(:,:,s,d));
                 
                 PV.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                 PV.invFI(:,:,s,d) = inv(PV.FI(:,:,s,d))+eps;  % variance and covariance of 
                                                                       % local bias removed estimates    
                 PV.detFI(s,d) = (det(PV.FI(:,:,s,d)));
                 PV.FIdiag(1,s,d)=B(1:2,1)'*inv(V(1:2,1:2))*B(1:2,1)+priorFI(1,1);
                 PV.FIdiag(2,s,d)=B(3:4,2)'*inv(V(3:4,3:4))*B(3:4,2)+priorFI(2,2);
                 PV.detFIdiag(s,d) = prod(PV.FIdiag(:,s,d));
                 
                 PV.MI = PV.MI - 1/2*log2((2*pi*exp(1))^2/PV.detFI(s,d))*pjoint(s,d);
                 
                 PV.A=sqrtm(squeeze(PV.var(:,:,s,d))); %4x4
                 
                 temp=squeeze(PV.mu(:,s,d)) + PV.A*randn(4,100000);
                 
                 temp=[temp(1,:)./(temp(2,:)+eps);angle((temp(3,:)+sqrt(-1)*temp(4,:))*exp(-sqrt(-1)*dvals(d)))+dvals(d)];
                 
                 PV.est.mean(:,s,d) = nanmean(temp')'; %2 x ns x nd
                 PV.est.var(:,:,s,d) = nancov(temp');
                                      
             end
             end
             
           

%end  bias correction
             
             for s=1:ns
             for d=1:nd
                 PV.est.bias(:,s,d) = PV.est.mean(:,s,d)-[svals(s);dvals(d)];  
             end
             end
             
             ds=svals(2)-svals(1);
             PV.est.gradmean=zeros(2,2,ns,nd);
             
%            Numerically estimate gradient of estimates mean (for FI calc)
             for s=2:ns-1
                 PV.est.gradmean(1:2,1,s,:) = (PV.est.mean(:,s+1,:)-PV.est.mean(:,s-1,:))/2/ds;
             end
             PV.est.gradmean(1:2,1,1,:) = (PV.est.mean(:,2,:)-PV.est.mean(:,1,:))/ds;
             PV.est.gradmean(1:2,1,ns,:) = (PV.est.mean(:,ns,:)-PV.est.mean(:,ns-1,:))/ds;
                          
             ds=dvals(2)-dvals(1);
             for s=1:ns
             for d=2:nd-1
                 PV.est.gradmean(:,2,s,d) = 1 + (PV.est.bias(:,s,d+1)-PV.est.bias(:,s,d-1))/ds/2;
             end
                 PV.est.gradmean(:,2,s,1) = 1 + (PV.est.bias(:,s,2)-PV.est.bias(:,s,nd))/ds/2;
                 PV.est.gradmean(:,2,s,nd) = 1 + (PV.est.bias(:,s,1)-PV.est.bias(:,s,nd-1))/ds/2;
             end
             
             PV.est.MI = ENTs + ENTd;
             PV.est.MIdiag = ENTs + ENTd;
             PV.s.mean=zeros(ns,1);
             PV.s.var=zeros(ns,1);
             PV.s.gradmean=zeros(ns,1);
             PV.d.mean=zeros(nd,1);
             PV.d.gradmean=zeros(nd,1);
             PV.d.var=zeros(nd,1);
             for s=1:ns
             for d=1:nd
                 
                 PV.est.err(:,:,s,d) = PV.est.var(:,:,s,d) + (PV.est.bias(:,s,d)*PV.est.bias(:,s,d)');
                 PV.est.deterr(s,d) = (det(PV.est.err(:,:,s,d)));
                 PV.est.corr(s,d) = PV.est.var(1,2,s,d)/sqrt(PV.est.var(1,1,s,d)*PV.est.var(2,2,s,d));
                 
                 V=squeeze(PV.est.var(:,:,s,d));
                 B=squeeze(PV.est.gradmean(:,:,s,d));
                 PV.est.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                 PV.est.invFI(:,:,s,d) = inv(PV.est.FI(:,:,s,d));  % variance and covariance of 
                                                                       % local bias removed estimates    
                 PV.est.detFI(s,d) = (det(PV.est.FI(:,:,s,d))+eps);
                 PV.est.FIdiag(:,s,d)=diag(B).^2./diag(PV.est.var(:,:,s,d))+diag(priorFI);
                 PV.est.detFIdiag(s,d) = prod(PV.est.FIdiag(:,s,d));
                 
                 PV.est.MI = PV.est.MI - 1/2*log2((2*pi*exp(1))^2/PV.est.detFI(s,d))*pjoint(s,d);
                 PV.est.MIdiag = PV.est.MIdiag - 1/2*log2((2*pi*exp(1))^2/PV.est.detFIdiag(s,d))*pjoint(s,d);

                 PV.s.mean(s) = PV.s.mean(s) + PV.est.mean(1,s,d)*pdvals(d);
                 PV.s.gradmean(s) = PV.s.gradmean(s) + PV.est.gradmean(1,1,s,d)*pdvals(d);
                 PV.s.var(s) = PV.s.var(s) + (PV.est.var(1,1,s,d) + PV.est.mean(1,s,d).^2)*pdvals(d);
                 
                 PV.d.mean(d) = PV.d.mean(d) + PV.est.mean(2,s,d)*psvals(s);
                 PV.d.gradmean(d) = PV.d.gradmean(d) + PV.est.gradmean(2,2,s,d)*psvals(s);
                 PV.d.var(d) = PV.d.var(d) + (PV.est.var(2,2,s,d) + PV.est.mean(2,s,d).^2)*psvals(s);
             end
             end
             
             PV.est.MIerr = - 1/2*sum(sum(log2((2*pi*exp(1))^2.*PV.est.deterr).*pjoint)) + ENTs + ENTd;
             
             PV.s.var=PV.s.var-PV.s.mean.^2;
             PV.s.bias=PV.s.mean - svals';
             PV.s.err = PV.s.var + PV.s.bias.^2;
             PV.s.FI = PV.s.gradmean.^2./PV.s.var + priorFIs;
             PV.s.invFI = 1./PV.s.FI;
             PV.s.MI = -1/2*sum(log2(2*pi*exp(1)./PV.s.FI).*psvals) + ENTs;
             PV.s.MIerr = -1/2*sum(log2(2*pi*exp(1).*PV.s.err).*psvals) + ENTs;
             PV.d.var=PV.d.var-PV.d.mean.^2;
             PV.d.bias=PV.d.mean - dvals';
             PV.d.err = PV.d.var + PV.d.bias.^2;
             PV.d.FI = PV.d.gradmean.^2./PV.d.var + priorFId;
             PV.d.invFI = 1./PV.d.FI;
             PV.d.MI = -1/2*sum(log2(2*pi*exp(1)./PV.d.FI).*pdvals) + ENTd;
             PV.d.MIerr = -1/2*sum(log2(2*pi*exp(1).*PV.d.err).*psvals) + ENTd;
             PV.MI = PV.est.MI;
             PV.Synergy = PV.MI - PV.s.MI - PV.d.MI;
             PV.StimSpecSynergy = -1/2*log2(1-PV.est.corr.^2);
             PV.meanSSS = sum(sum(PV.StimSpecSynergy.*pjoint));
             PV.ErrSynergy = PV.est.MIerr - PV.s.MIerr - PV.d.MIerr;
             
             % Marginal LOLE decoders (unbiased on average)
             Marg.ws=zeros(nn,ns);
             Marg.wd=zeros(nn,nd);
             
             NB.ws=zeros(nn,ns); % naieve bayesian deocder (no knowledge of correlations is neural response
                                 % or between s and d) 
             NB.wd=zeros(nn,ns); % formerally called diag.
                          
             for s=1:ns
                 idx=mus(:,s)>-Inf;
                 Marg.ws(idx,s) = Sigmas(idx,idx,s)\(dmusds(idx,s));  %marginal decoder
                 Marg.ws(idx,s) = Marg.ws(idx,s)/(Marg.ws(idx,s)'*dmusds(idx,s)); % 'unbiased' marginal decoder
                 NB.ws(idx,s) = diag(1./diag(Sigmas(idx,idx,s)))*dmusds(idx,s);
                 NB.ws(idx,s) = NB.ws(idx,s)/(NB.ws(idx,s)'*dmusds(idx,s));
                 %marginal decoder under assumption that conditional cov is diagonal
             end
             
             for d=1:nd
                 idx=mud(:,d)>-Inf;
                 Marg.wd(idx,d) = Sigmad(idx,idx,d)\(dmuddd(idx,d));
                 Marg.wd(idx,d) = Marg.wd(idx,d)/(Marg.wd(idx,d)'*dmuddd(idx,d)); % 'unbiased' marginal decoder
                 NB.wd(idx,d) = diag(1./diag(Sigmad(idx,idx,d)))*dmuddd(idx,d);
                 NB.wd(idx,d) = NB.wd(idx,d)/(NB.wd(idx,d)'*dmuddd(idx,d));
             end
             
             Marg.w=zeros(nn,2,ns,nd);
             Marg.est.FI=zeros(2,2,ns,nd);
             Marg.est.MI = ENTs + ENTd;
             Marg.est.MIdiag = ENTs + ENTd;
             Marg.s.mean=zeros(ns,1);
             Marg.s.var=zeros(ns,1);
             Marg.s.gradmean=zeros(ns,1);
             Marg.d.mean=zeros(nd,1);
             Marg.d.gradmean=zeros(nd,1);
             Marg.d.var=zeros(nd,1);
             
             for s=1:ns
             for d=1:nd
                 Marg.w(:,:,s,d) = [Marg.ws(:,s),Marg.wd(:,d)];
                 Marg.est.mean(1,s,d) = Marg.ws(:,s)'*(mu(:,s,d)-mus(:,s)) + svals(s);
                 Marg.est.mean(2,s,d) = Marg.wd(:,d)'*(mu(:,s,d)-mud(:,d)) + dvals(d);
                 Marg.est.bias(:,s,d) = Marg.est.mean(:,s,d)-[svals(s);dvals(d)];
                 
                 Marg.est.gradmean(:,:,s,d) = [Marg.w(:,:,s,d)'*dmuds(:,s,d),Marg.w(:,:,s,d)'*dmudd(:,s,d)];

                 Marg.est.var(:,:,s,d) = Marg.w(:,:,s,d)'*Sigma(:,:,s,d)*Marg.w(:,:,s,d);
                 Marg.est.err(:,:,s,d) = Marg.est.var(:,:,s,d) + (Marg.est.bias(:,s,d)*Marg.est.bias(:,s,d)');
                 Marg.est.deterr(s,d) = det(Marg.est.err(:,:,s,d));
                 Marg.est.corr(s,d) = Marg.est.var(1,2,s,d)/sqrt(Marg.est.var(1,1,s,d)*Marg.est.var(2,2,s,d));
                 
                 V=squeeze(Marg.est.var(:,:,s,d));
                 B=squeeze(Marg.est.gradmean(:,:,s,d));
%                 B=eye(2);
                 Marg.est.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                 Marg.est.invFI(:,:,s,d) = inv(Marg.est.FI(:,:,s,d));  % variance and covariance of 
                                                                       % local bias removed estimates    
                 Marg.est.detFI(s,d) = (det(Marg.est.FI(:,:,s,d)));
                 Marg.est.FIdiag(:,s,d)=diag(B).^2./diag(Marg.est.var(:,:,s,d))+diag(priorFI);
                 Marg.est.detFIdiag(s,d) = prod(Marg.est.FIdiag(:,s,d));
                 Marg.est.MI = Marg.est.MI - 1/2*log2((2*pi*exp(1))^2/Marg.est.detFI(s,d))*pjoint(s,d);
                 Marg.est.MIdiag = Marg.est.MIdiag - 1/2*log2((2*pi*exp(1))^2/Marg.est.detFIdiag(s,d))*pjoint(s,d);

                 Marg.s.mean(s) = Marg.s.mean(s) + Marg.est.mean(1,s,d)*pdvals(d);
                 Marg.s.gradmean(s) = Marg.s.gradmean(s) + Marg.est.gradmean(1,1,s,d)*pdvals(d);
%                 Marg.s.var(s) = Marg.s.var(s) + (Marg.est.var(1,1,s,d) + Marg.est.mean(1,s,d).^2)*pdvals(d);
                 Marg.s.var(s) = Marg.s.var(s) + (Marg.est.var(1,1,s,d) + Marg.est.mean(1,s,d).^2)*pdvals(d);
                 
                 Marg.d.mean(d) = Marg.d.mean(d) + Marg.est.mean(2,s,d)*psvals(s);
                 Marg.d.gradmean(d) = Marg.d.gradmean(d) + Marg.est.gradmean(2,2,s,d)*psvals(s);
                 Marg.d.var(d) = Marg.d.var(d) + (Marg.est.var(2,2,s,d) + Marg.est.mean(2,s,d).^2)*psvals(s);
%                 Marg.d.var(d) = Marg.d.var(d) + (Marg.est.var(2,2,s,d) + Marg.est.mean(2,s,d).^2)*psvals(s);
             end
             end
             
             Marg.est.MIerr = -1/2*sum(sum(log2((2*pi*exp(1))^2.*Marg.est.deterr).*pjoint)) + ENTs + ENTd;
             
%             Marg.s.m2 = Marg.s.var;
             Marg.s.var=Marg.s.var-Marg.s.mean.^2;
             Marg.s.bias=Marg.s.mean - svals';
             Marg.s.err = Marg.s.var + Marg.s.bias.^2;
             Marg.s.FI = Marg.s.gradmean.^2./Marg.s.var + priorFIs;
             Marg.s.invFI = 1./Marg.s.FI;
             Marg.s.MI = -1/2*sum(log2(2*pi*exp(1)./Marg.s.FI).*psvals) + ENTs;
             Marg.s.MIerr = -1/2*sum(log2(2*pi*exp(1).*Marg.s.err).*psvals) + ENTs;
                          
             Marg.d.var=Marg.d.var-Marg.d.mean.^2;
             Marg.d.bias=Marg.d.mean - dvals';
             Marg.d.err = Marg.d.var + Marg.d.bias.^2;
             Marg.d.FI = Marg.d.gradmean.^2./Marg.d.var + priorFId;
             Marg.d.invFI = 1./Marg.d.FI;
             Marg.d.MI = -1/2*sum(log2(2*pi*exp(1)./Marg.d.FI).*pdvals) + ENTd;
             Marg.d.MIerr = -1/2*sum(log2(2*pi*exp(1).*Marg.d.err).*pdvals) + ENTd;
             Marg.MI = Marg.est.MI;
             Marg.Synergy = Marg.est.MI - Marg.s.MI - Marg.d.MI;
             Marg.StimSpecSynergy = -1/2*log2(1-Marg.est.corr.^2);
             Marg.meanSSS = sum(sum(Marg.StimSpecSynergy.*pjoint));
             Marg.ErrSynergy = Marg.est.MIerr - Marg.s.MIerr - Marg.d.MIerr;
                         
             
             NB.w=zeros(nn,2,ns,nd);
             NB.est.FI=zeros(2,2,ns,nd);             
             NB.est.MI = ENTs + ENTd;
             NB.est.MIdiag = ENTs + ENTd;
             NB.s.mean=zeros(ns,1);
             NB.s.var=zeros(ns,1);
             NB.s.gradmean=zeros(ns,1);
             NB.d.mean=zeros(nd,1);
             NB.d.gradmean=zeros(nd,1);
             NB.d.var=zeros(nd,1);

             for s=1:ns
             for d=1:nd
                 NB.w(:,:,s,d) = [NB.ws(:,s),NB.wd(:,d)];
                 NB.est.mean(1,s,d) = NB.ws(:,s)'*(mu(:,s,d)-mus(:,s)) + svals(s);
                 NB.est.mean(2,s,d) = NB.wd(:,d)'*(mu(:,s,d)-mud(:,d)) + dvals(d);
                 NB.est.bias(:,s,d) = NB.est.mean(:,s,d)-[svals(s);dvals(d)];
                 
                 NB.est.gradmean(:,:,s,d) = [NB.w(:,:,s,d)'*dmuds(:,s,d),NB.w(:,:,s,d)'*dmudd(:,s,d)];
                 NB.est.var(:,:,s,d) = NB.w(:,:,s,d)'*Sigma(:,:,s,d)*NB.w(:,:,s,d);
                 NB.est.err(:,:,s,d) = NB.est.var(:,:,s,d) + (NB.est.bias(:,s,d)*NB.est.bias(:,s,d)');
                 NB.est.deterr(s,d) = det(NB.est.err(:,:,s,d));
                 
                 NB.est.corr(s,d) = NB.est.var(1,2,s,d)/sqrt(NB.est.var(1,1,s,d)*NB.est.var(2,2,s,d));
                 
                 V=squeeze(NB.est.var(:,:,s,d));
                 B=squeeze(NB.est.gradmean(:,:,s,d));
%                 B=eye(2);
                 NB.est.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                 NB.est.invFI(:,:,s,d) = inv(NB.est.FI(:,:,s,d));  % variance and covariance of 
                                                                       % local bias removed estimates    
                 NB.est.detFI(s,d) = (det(NB.est.FI(:,:,s,d)));
                 NB.est.FIdiag(:,s,d)=diag(B).^2./diag(NB.est.var(:,:,s,d)) + diag(priorFI);
                 NB.est.detFIdiag(s,d) = prod(NB.est.FIdiag(:,s,d));
                 
                 NB.est.MI = NB.est.MI - 1/2*log2((2*pi*exp(1))^2/NB.est.detFI(s,d))*pjoint(s,d);
                 NB.est.MIdiag = NB.est.MIdiag - 1/2*log2((2*pi*exp(1))^2/NB.est.detFIdiag(s,d))*pjoint(s,d);

                 NB.s.mean(s) = NB.s.mean(s) + NB.est.mean(1,s,d)*pdvals(d);
                 NB.s.gradmean(s) = NB.s.gradmean(s) + NB.est.gradmean(1,1,s,d)*pdvals(d);
                 NB.s.var(s) = NB.s.var(s) + (NB.est.var(1,1,s,d) + NB.est.mean(1,s,d).^2)*pdvals(d);
                 
                 NB.d.mean(d) = NB.d.mean(d) + NB.est.mean(2,s,d)*psvals(s);
                 NB.d.gradmean(d) = NB.d.gradmean(d) + NB.est.gradmean(2,2,s,d)*psvals(s);
                 NB.d.var(d) = NB.d.var(d) + (NB.est.var(2,2,s,d) + NB.est.mean(2,s,d).^2)*psvals(s);
                 
             end
             end
             
             NB.est.MIerr = -1/2*sum(sum(log2((2*pi*exp(1))^2.*NB.est.deterr).*pjoint)) + ENTs + ENTd;
             
             NB.s.var=NB.s.var-NB.s.mean.^2;
             NB.s.bias=NB.s.mean - svals';
             NB.s.err = NB.s.var + NB.s.bias.^2;
             NB.s.FI = NB.s.gradmean.^2./NB.s.var + priorFIs;
             NB.s.invFI = 1./NB.s.FI;
             NB.s.MI = -1/2*sum(log2(2*pi*exp(1)./NB.s.FI).*psvals) + ENTs;
             NB.s.MIerr = -1/2*sum(log2(2*pi*exp(1).*NB.s.err).*psvals) + ENTs;
             
             NB.d.var=NB.d.var-NB.d.mean.^2;
             NB.d.bias=NB.d.mean - dvals';
             NB.d.err = NB.d.var + NB.d.bias.^2;
             NB.d.FI = NB.d.gradmean.^2./NB.d.var + priorFId;
             NB.d.invFI = 1./NB.d.FI;
             NB.d.MI = -1/2*sum(log2(2*pi*exp(1)./NB.d.FI).*pdvals) + ENTd;
             NB.d.MIerr = -1/2*sum(log2(2*pi*exp(1).*NB.d.err).*pdvals) + ENTd;
             NB.MI = NB.est.MI;
             NB.Synergy = NB.est.MI - NB.s.MI - NB.d.MI;
             NB.StimSpecSynergy = -1/2*log2(1-NB.est.corr.^2);
             NB.meanSSS = sum(sum(NB.StimSpecSynergy.*pjoint));
             NB.ErrSynergy = NB.est.MIerr - NB.s.MIerr - NB.d.MIerr;
             
             
             % Naieve Bayes is joint completely unbiased but knows nothing
             % about correlations in neural responses...
             
             NBjoint.ws=zeros(nn,ns,nd); %knows about correlations in s and d but not in neural response
             NBjoint.wd=zeros(nn,ns,nd);
             NBjoint.w=zeros(nn,2,ns,nd);
             NBjoint.est.MI = ENTs + ENTd;
             NBjoint.est.MIdiag = ENTs + ENTd;
             NBjoint.s.mean=zeros(ns,1);
             NBjoint.s.var=zeros(ns,1);
             NBjoint.s.gradmean=zeros(ns,1);
             NBjoint.d.mean=zeros(nd,1);
             NBjoint.d.gradmean=zeros(nd,1);
             NBjoint.d.var=zeros(nd,1);

             for s=1:ns
                 for d=1:nd
 %                    idx=mu(:,s,d)>0;
                     
                     NBjoint.ws(idx,s,d)=dmuds(idx,s,d)./diag(Sigma(idx,idx,s,d)); %locally biased
                     NBjoint.wd(idx,s,d)=dmudd(idx,s,d)./diag(Sigma(idx,idx,s,d)); 
                     NBjoint.w(:,1,s,d)=NBjoint.ws(:,s,d);
                     NBjoint.w(:,2,s,d)=NBjoint.wd(:,s,d);

                     NBjoint.est.mean(:,s,d)=NBjoint.w(:,:,s,d)'*mu(:,s,d);
                     NBjoint.est.gradmean(:,:,s,d) = [NBjoint.w(:,:,s,d)'*dmuds(:,s,d),NBjoint.w(:,:,s,d)'*dmudd(:,s,d)];
                     NBjoint.est.var(:,:,s,d) = NBjoint.w(:,:,s,d)'*Sigma(:,:,s,d)*NBjoint.w(:,:,s,d);
                 
                     V=squeeze(NBjoint.est.var(:,:,s,d));
                     B=squeeze(NBjoint.est.gradmean(:,:,s,d));
                 
                     NBjoint.est.FI(:,:,s,d) = B'*inv(V)*B + priorFI;  %Assumes local bias removal applied to estimated s,d                 
                     NBjoint.est.invFI(:,:,s,d) = inv(NB.est.FI(:,:,s,d));  % variance and covariance of 
                                                                           % local bias removed estimates                         
                     NBjoint.est.detFI(s,d) = (det(NB.est.FI(:,:,s,d)));
                     NBjoint.est.FIdiag(:,s,d)=diag(B).^2./diag(NB.est.var(:,:,s,d)) + diag(priorFI);
                     NBjoint.est.detFIdiag(s,d) = prod(NBjoint.est.FIdiag(:,s,d));                     
                     
                     NBjoint.est.mean(:,s,d) = [svals(s);dvals(d);];
                     NBjoint.gradmean(:,:,s,d) = eye(2);
                     NBjoint.est.var(:,:,s,d) = NBjoint.est.invFI(:,:,s,d);
                     NBjoint.est.ERR(:,:,s,d) = NBjoint.est.invFI(:,:,s,d);
                     NBjoint.est.deterr(s,d) = det(NBjoint.est.ERR(:,:,s,d));
                     
                     NBjoint.s.mean(s) = svals(s);
                     NBjoint.s.gradmean(s) = 1;
                     NBjoint.s.var(s) = NBjoint.s.var(s) + NBjoint.est.var(1,1,s,d)*pdvals(d);

                     NBjoint.d.mean(d) = dvals(d);
                     NBjoint.d.gradmean(d) = 1;
                     NBjoint.d.var(d) = NBjoint.d.var(d) + NBjoint.est.var(2,2,s,d)*psvals(s);
                     
                     NBjoint.est.mean(1:2,s,d)=[svals(s);dvals(d)];
                     NBjoint.est.bias(1:2,s,d)=0;
                     NBjoint.est.gradmean(:,:,s,d)=1;
                     NBjoint.est.var(:,:,s,d)=NBjoint.est.invFI(:,:,s,d);
                     
                     NBjoint.est.corr(s,d)=NBjoint.est.var(1,2,s,d)/sqrt(NBjoint.est.var(1,1,s,d)*NBjoint.est.var(2,2,s,d));
                     NBjoint.est.MI = NBjoint.est.MI - 1/2*log2((2*pi*exp(1))^2/NBjoint.est.detFI(s,d))*pjoint(s,d);   
                     NBjoint.est.MIdiag = NBjoint.est.MIdiag - 1/2*log2((2*pi*exp(1))^2/NBjoint.est.detFIdiag(s,d))*pjoint(s,d);
                     
                 end
             end
             NBjoint.est.MIerr = -1/2*sum(sum(log2((2*pi*exp(1))^2.*NBjoint.est.deterr).*pjoint)) + ENTs + ENTd;
             
             NBjoint.s.bias = zeros(ns,1);
             NBjoint.s.err = NBjoint.s.var;
             NBjoint.s.FI = 1./NBjoint.s.var;
             NBjoint.s.invFI = NBjoint.s.var;
             NBjoint.s.MI = -1/2*sum(log2(2*pi*exp(1)./NBjoint.s.FI).*psvals) + ENTs;
             NBjoint.s.MIerr = -1/2*sum(log2(2*pi*exp(1).*NBjoint.s.err).*psvals) + ENTs;
             
             NBjoint.d.bias=zeros(nd,1);
             NBjoint.d.err = NBjoint.d.var;
             NBjoint.d.FI = 1./NBjoint.d.var;
             NBjoint.d.invFI = NBjoint.d.var;
             NBjoint.d.MI = -1/2*sum(log2(2*pi*exp(1)./NBjoint.d.FI).*pdvals) + ENTd;
             NBjoint.d.MIerr = -1/2*sum(log2(2*pi*exp(1).*NBjoint.d.err).*pdvals) + ENTd;
             
             NBjoint.MI = NB.est.MI;

             NBjoint.Synergy = NBjoint.est.MI - NBjoint.s.MI - NBjoint.d.MI;
             NBjoint.ErrSynergy = NBjoint.est.MIerr - NBjoint.s.MIerr - NBjoint.d.MIerr;
             
             NBNBjoint.ErrSynergy = NBjoint.est.MIerr - NB.s.MIerr - NB.d.MIerr;
           
             NBjoint.StimSpecSynergy = -1/2*log2(1-NBjoint.est.corr.^2);
             NBjoint.meanSSS = sum(sum(NBjoint.StimSpecSynergy.*pjoint));
             
             NBSynergy = NBjoint.est.MI - NB.s.MI - NB.d.MI;
             NBNBjoint.Synergy = NBjoint.est.MI - NB.s.MI - NB.d.MI;
             
             LOLE.dhds=zeros(nn,ns,nd);
             LOLE.dhdd=zeros(nn,ns,nd);
             LOLE.MI = ENTs + ENTd;
             LOLE.est.MIdiag = ENTs + ENTd;
             LOLE.s.mean=svals';
             LOLE.s.var=zeros(ns,1);
             LOLE.s.gradmean=zeros(ns,1);
             LOLE.d.mean=dvals';
             LOLE.d.gradmean=zeros(nd,1);
             LOLE.d.var=zeros(nd,1);
             LOLE.ws=zeros(nn,ns,nd);
             LOLE.wd=zeros(nn,ns,nd);
             LOLE.w=zeros(nn,2,ns,nd);

             for s=1:ns
                 for d=1:nd             
                     LOLE.dhds(idx,s,d)=(Sigma(idx,idx,s,d)\dmuds(idx,s,d)); %locally optimal joint decoder weights
                     LOLE.dhdd(idx,s,d)=(Sigma(idx,idx,s,d)\dmudd(idx,s,d));
                     
                     LOLE.FI(1,1,s,d) = dmuds(idx,s,d)'*LOLE.dhds(idx,s,d);
                     LOLE.FI(1,2,s,d) = dmuds(idx,s,d)'*LOLE.dhdd(idx,s,d);
                     LOLE.FI(2,1,s,d) = dmudd(idx,s,d)'*LOLE.dhds(idx,s,d);
                     LOLE.FI(2,2,s,d) = dmudd(idx,s,d)'*LOLE.dhdd(idx,s,d);
                     LOLE.detFI(s,d) = (det(LOLE.FI(:,:,s,d)));                     
                     LOLE.invFI(:,:,s,d)=inv(LOLE.FI(:,:,s,d));
                     
                     LOLE.ws(:,s,d)=LOLE.invFI(1,1,s,d)*LOLE.dhds(:,s,d)+LOLE.invFI(1,2,s,d)*LOLE.dhdd(:,s,d);
                     LOLE.wd(:,s,d)=LOLE.invFI(2,2,s,d)*LOLE.dhdd(:,s,d)+LOLE.invFI(1,2,s,d)*LOLE.dhds(:,s,d);
                     LOLE.w(:,:,s,d)=[LOLE.ws(:,s,d),LOLE.wd(:,s,d)];
                                          
                     LOLE.est.var(:,:,s,d)=LOLE.invFI(:,:,s,d);
                     LOLE.est.mean(1:2,s,d)=[svals(s);dvals(d)];
                     LOLE.est.bias(1:2,s,d)=0;               
                     LOLE.est.gradmean(:,:,s,d) = [LOLE.w(:,:,s,d)'*dmuds(:,s,d),LOLE.w(:,:,s,d)'*dmudd(:,s,d)];
                     LOLE.est.ERR(:,:,s,d)=LOLE.est.var(:,:,s,d);
                     LOLE.est.corr(s,d)=LOLE.est.var(1,2,s,d)/sqrt(LOLE.est.var(1,1,s,d)*LOLE.est.var(2,2,s,d));
                     LOLE.est.FIdiag(:,s,d)=1./diag(LOLE.invFI(:,:,s,d));
                     LOLE.est.detFIdiag(s,d) = prod(LOLE.est.FIdiag(:,s,d));
                     
                     LOLE.s.gradmean(s) = LOLE.s.gradmean(s) + LOLE.est.gradmean(1,1,s,d)*pdvals(d);
                     LOLE.s.var(s) = LOLE.s.var(s) + (LOLE.est.var(1,1,s,d))*pdvals(d);

                     LOLE.d.gradmean(d) = LOLE.d.gradmean(d) + LOLE.est.gradmean(2,2,s,d)*psvals(s);
                     LOLE.d.var(d) = LOLE.d.var(d) + (LOLE.est.var(2,2,s,d))*psvals(s);
                     
                     LOLE.FI(:,:,s,d) = LOLE.FI(:,:,s,d) + priorFI;
                     LOLE.detFI(s,d) = (det(LOLE.FI(:,:,s,d)));                     
                     LOLE.invFI(:,:,s,d)=inv(LOLE.FI(:,:,s,d));
                     LOLE.FI(:,:,s,d) = LOLE.FI(:,:,s,d) + priorFI;
                     
                     LOLE.MI = LOLE.MI - 1/2*log2((2*pi*exp(1))^2/LOLE.detFI(s,d))*pjoint(s,d);
                     LOLE.est.MIdiag = LOLE.est.MIdiag - 1/2*log2((2*pi*exp(1))^2/LOLE.est.detFIdiag(s,d))*pjoint(s,d);

                 end
             end
            
            
             LOLE.s.err = LOLE.s.var;
             LOLE.s.FI = LOLE.s.gradmean.^2./LOLE.s.var + priorFIs;
             LOLE.s.invFI = 1./LOLE.s.FI;
             LOLE.s.MI = -1/2*sum(log2(2*pi*exp(1)./LOLE.s.FI).*psvals) + ENTs;
             
             
             LOLE.d.err = LOLE.d.var;
             LOLE.d.FI = LOLE.d.gradmean.^2./LOLE.d.var + priorFId;
             LOLE.d.invFI = 1./LOLE.d.FI;
             LOLE.d.MI = (-1/2*sum(log2(2*pi*exp(1)./LOLE.d.FI).*pdvals)) + ENTd;
             LOLE.est.MI = LOLE.MI;
             
             LOLE.Synergy = LOLE.MI - LOLE.d.MI - LOLE.s.MI;
             LOLE.StimSpecSynergy = -1/2*log2(1-LOLE.est.corr.^2);
             LOLE.meanSSS = sum(sum(LOLE.StimSpecSynergy.*pjoint));
             LOLEMarg.Synergy = LOLE.MI - Marg.d.MI - Marg.s.MI;
             LOLEMarg.ErrSynergy = LOLE.MI - Marg.d.MIerr - Marg.s.MIerr;
             
%             Synergy = LOLE.MI - Marg.s.MI - Marg.d.MI;
             
             if(diffcorron==1)
                 MIsmax(tidx,nidx) = -1/2*log2(2*pi*exp(1)/behFIs) + ENTs;
                 MIdmax(tidx,nidx) = -1/2*log2(2*pi*exp(1)/behFId) + ENTd;
                 MIjointmax(tidx,nidx) = -1/2*log2((2*pi*exp(1))^2/det(inv(epsilon))) + ENTs + ENTd;
                 SynergyIn(tidx,nidx) = MIjointmax(tidx,nidx) - MIsmax(tidx,nidx) - MIdmax(tidx,nidx);
                 MImax(tidx,nidx)=MIsmax(tidx,nidx)+MIdmax(tidx,nidx);
             else
                 MIsmax(tidx,nidx)=LOLE.s.MI;
                 MIdmax(tidx,nidx)=LOLE.d.MI;
                 MIjointmax(tidx,nidx)=LOLE.MI;
                 SynergyIn(tidx,nidx) = MIjointmax(tidx,nidx) - MIsmax(tidx,nidx) - MIdmax(tidx,nidx);
                 MImax(tidx,nidx)=MIjointmax(tidx,nidx);
             end             
             
% Compute Summary Stats.             
             
             PV.s.meanFI = nanmean(PV.s.FI);
             PV.s.meanERR = nanmean(PV.s.err);
             PV.s.meanVar = nanmean(PV.s.var);
             PV.d.meanFI = nanmean(PV.d.FI);
             PV.d.meanERR = nanmean(PV.d.err);
             PV.d.meanVar = nanmean(PV.d.var);             
             
             GL.s.meanFI = mean(GL.s.FI);
             GL.s.meanERR = mean(GL.s.err);
             GL.s.meanVar = mean(GL.s.var);
             GL.d.meanFI = mean(GL.d.FI);
             GL.d.meanERR = mean(GL.d.err);
             GL.d.meanVar = mean(GL.d.var);             
             
             NB.s.meanFI = mean(NB.s.FI);
             NB.s.meanERR = mean(NB.s.err);
             NB.s.meanVar = mean(NB.s.var);
             NB.d.meanFI = mean(NB.d.FI);
             NB.d.meanERR = mean(NB.d.err);
             NB.d.meanVar = mean(NB.d.var);             
             
             NBjoint.s.meanFI = mean(NBjoint.s.FI);
             NBjoint.s.meanERR = mean(NBjoint.s.err);
             NBjoint.s.meanVar = mean(NBjoint.s.var);
             NBjoint.d.meanFI = mean(NBjoint.d.FI);
             NBjoint.d.meanERR = mean(NBjoint.d.err);
             NBjoint.d.meanVar = mean(NBjoint.d.var);             
             
             Marg.s.meanFI = mean(Marg.s.FI);
             Marg.s.meanERR = mean(Marg.s.err);
             Marg.s.meanVar = mean(Marg.s.var);
             Marg.d.meanFI = mean(Marg.d.FI);
             Marg.d.meanERR = mean(Marg.d.err);
             Marg.d.meanVar = mean(Marg.d.var);             
             
             LOLE.s.meanFI = mean(LOLE.s.FI);
             LOLE.s.meanERR = mean(LOLE.s.err);
             LOLE.s.meanVar = mean(LOLE.s.var);
             LOLE.d.meanFI = mean(LOLE.d.FI);
             LOLE.d.meanERR = mean(LOLE.d.err);
             LOLE.d.meanVar = mean(LOLE.d.var);             
             
             %pjoint = psvals*pdvals'; %pstim %%change this to generalize to correlated stims!!
             normd = maxFId; %dir thresh independent of dir, speed
             norms = maxFIs; %repmat(maxFIs*svals.^2,nd,1);  %ns x nd speed thresh is speed dependent

             MIerrstat.GL(ridx)=GL.est.MIerr;
             MIserrstat.GL(ridx)=GL.s.MIerr;
             MIderrstat.GL(ridx)=GL.d.MIerr;
             MIerrstat.PV(ridx)=PV.est.MIerr;
             MIserrstat.PV(ridx)=PV.s.MIerr;
             MIderrstat.PV(ridx)=PV.d.MIerr;
             MIerrstat.NB(ridx)=NB.est.MIerr;
             MIserrstat.NB(ridx)=NB.s.MIerr;
             MIderrstat.NB(ridx)=NB.d.MIerr;
             MIerrstat.NBjoint(ridx)=NBjoint.est.MIerr;
             MIserrstat.NBjoint(ridx)=NBjoint.s.MIerr;
             MIderrstat.NBjoint(ridx)=NBjoint.d.MIerr;
             MIerrstat.Marg(ridx)=Marg.est.MIerr;
             MIserrstat.Marg(ridx)=Marg.s.MIerr;
             MIderrstat.Marg(ridx)=Marg.d.MIerr;

             
             MIstat.GL(ridx)=GL.MI;
             MIs_stat.GL(ridx)=GL.s.MI;
             MId_stat.GL(ridx)=GL.d.MI;
             FIs_stat.GL(ridx)=GL.s.meanFI;
             FId_stat.GL(ridx)=GL.d.meanFI;             
             ERRs_stat.GL(ridx)=GL.s.meanERR;
             Vars_stat.GL(ridx)=GL.s.meanVar;
             ERRd_stat.GL(ridx)=GL.d.meanERR;
             Vard_stat.GL(ridx)=GL.d.meanVar;
             Synergy_stat.GL(ridx)=GL.Synergy;
             ErrSynergy_stat.GL(ridx)=GL.ErrSynergy;
             FracVars_stat.GL(ridx) = mean(2.^GL.s.var-1);
             degVard_stat.GL(ridx)=Vard_stat.GL(ridx)*(180/pi)^2;
             
             MIstat.PV(ridx)=PV.MI;
             MIs_stat.PV(ridx)=PV.s.MI;
             MId_stat.PV(ridx)=PV.d.MI;
             FIs_stat.PV(ridx)=PV.s.meanFI;
             FId_stat.PV(ridx)=PV.d.meanFI;             
             ERRs_stat.PV(ridx)=PV.s.meanERR;
             Vars_stat.PV(ridx)=PV.s.meanVar;
             ERRd_stat.PV(ridx)=PV.d.meanERR;
             Vard_stat.PV(ridx)=PV.d.meanVar;
             Synergy_stat.PV(ridx)=PV.Synergy;
             ErrSynergy_stat.PV(ridx)=PV.ErrSynergy;
             FracVars_stat.PV(ridx) = mean(2.^PV.s.var-1);
             degVard_stat.PV(ridx)=Vard_stat.PV(ridx)*(180/pi)^2;
             
             MIstat.NB(ridx)=NB.MI;
             MIs_stat.NB(ridx)=NB.s.MI;
             MId_stat.NB(ridx)=NB.d.MI;
             MIs_stat.NB(ridx)=NB.s.MI;
             MId_stat.NB(ridx)=NB.d.MI;
             FIs_stat.NB(ridx)=NB.s.meanFI;
             FId_stat.NB(ridx)=NB.d.meanFI;             
             ERRs_stat.NB(ridx)=NB.s.meanERR;
             Vars_stat.NB(ridx)=NB.s.meanVar;
             ERRd_stat.NB(ridx)=NB.d.meanERR;
             Vard_stat.NB(ridx)=NB.d.meanVar;
             Synergy_stat.NB(ridx)=NB.Synergy;
             ErrSynergy_stat.NB(ridx)=NB.ErrSynergy;
             FracVars_stat.NB(ridx) = mean(2.^NB.s.var-1);
             degVard_stat.NB(ridx)=Vard_stat.NB(ridx)*(180/pi)^2;
             
             MIstat.NBjoint(ridx)=NBjoint.MI;
             MIs_stat.NBjoint(ridx)=NBjoint.s.MI;
             MId_stat.NBjoint(ridx)=NBjoint.d.MI;
             FIs_stat.NBjoint(ridx)=NBjoint.s.meanFI;
             FId_stat.NBjoint(ridx)=NBjoint.d.meanFI;             
             ERRs_stat.NBjoint(ridx)=NBjoint.s.meanERR;
             Vars_stat.NBjoint(ridx)=NBjoint.s.meanVar;
             ERRd_stat.NBjoint(ridx)=NBjoint.d.meanERR;
             Vard_stat.NBjoint(ridx)=NBjoint.d.meanVar;
             Synergy_stat.NBjoint(ridx)=NBjoint.Synergy;
             ErrSynergy_stat.NBjoint(ridx)=NBjoint.ErrSynergy;
             FracVars_stat.NBjoint(ridx) = mean(2.^NBjoint.s.var-1);
             degVard_stat.NBjoint(ridx)=Vard_stat.NBjoint(ridx)*(180/pi)^2;
             
              ErrSynergy_stat.NBNBjoint(ridx)=NBNBjoint.ErrSynergy;
             
             MIstat.Marg(ridx)=Marg.MI;
             MIs_stat.Marg(ridx)=Marg.s.MI;
             MId_stat.Marg(ridx)=Marg.d.MI;
             FIs_stat.Marg(ridx)=Marg.s.meanFI;
             FId_stat.Marg(ridx)=Marg.d.meanFI;             
             ERRs_stat.Marg(ridx)=Marg.s.meanERR;
             Vars_stat.Marg(ridx)=Marg.s.meanVar;
             ERRd_stat.Marg(ridx)=Marg.d.meanERR;
             Vard_stat.Marg(ridx)=Marg.d.meanVar;
             Synergy_stat.Marg(ridx)=Marg.Synergy;
             ErrSynergy_stat.Marg(ridx)=Marg.ErrSynergy;
             FracVars_stat.Marg(ridx) = mean(2.^Marg.s.var-1);
             degVard_stat.Marg(ridx)=Vard_stat.Marg(ridx)*(180/pi)^2;
             
             MIstat.LOLE(ridx)=LOLE.MI;
             MIs_stat.LOLE(ridx)=LOLE.s.MI;
             MId_stat.LOLE(ridx)=LOLE.d.MI;
             FIs_stat.LOLE(ridx)=LOLE.s.meanFI;
             FId_stat.LOLE(ridx)=LOLE.d.meanFI;             
             ERRs_stat.LOLE(ridx)=LOLE.s.meanERR;
             Vars_stat.LOLE(ridx)=LOLE.s.meanVar;
             ERRd_stat.LOLE(ridx)=LOLE.d.meanERR;
             Vard_stat.LOLE(ridx)=LOLE.d.meanVar;
             Synergy_stat.LOLE(ridx)=LOLE.Synergy;
             ErrSynergy_stat.LOLE(ridx)=LOLE.Synergy;
             FracVars_stat.LOLE(ridx) = mean(2.^LOLE.s.var-1);
             degVard_stat.LOLE(ridx)=Vard_stat.LOLE(ridx)*(180/pi)^2;
             
             ErrSynergy_stat.LOLEMarg(ridx)=LOLEMarg.ErrSynergy;
             
             Synergy_stat.LOLEMarg(ridx) = LOLE.MI - Marg.s.MI - Marg.d.MI;
             Synergy_stat.NBNBjoint(ridx) = NBjoint.est.MI - NB.s.MI - NB.d.MI;
             
             SI_stat.GL(ridx)=1-mean(GL.est.corr(:).^2);
             SI_stat.PV(ridx)=1-mean(PV.est.corr(:).^2);
             SI_stat.NB(ridx)=1-mean(NB.est.corr(:).^2);
             SI_stat.NBjoint(ridx)=1-mean(NBjoint.est.corr(:).^2);
             SI_stat.Marg(ridx)=1-mean(Marg.est.corr(:).^2);
             SI_stat.LOLE(ridx)=1-mean(LOLE.est.corr(:).^2);

             FracSyn_stat.GL(ridx)=GL.Synergy/(GL.s.MI+GL.d.MI);
             FracSyn_stat.PV(ridx)=PV.Synergy/(PV.s.MI+PV.d.MI);
             FracSyn_stat.NB(ridx)=NB.Synergy/(NB.s.MI+NB.d.MI);
             FracSyn_stat.NBjoint(ridx)=NBjoint.Synergy/(NBjoint.s.MI+NBjoint.d.MI);
             FracSyn_stat.Marg(ridx)=Marg.Synergy/(Marg.s.MI+Marg.d.MI);
             FracSyn_stat.LOLE(ridx)=LOLE.Synergy/(LOLE.s.MI+LOLE.d.MI);

             %FracSyn_stat.LOLEMarg(ridx) = Synergy_stat.LOLEMarg(ridx)/(LOLE.s.MI+LOLE.d.MI);
             %FracSyn_stat.NBNBjoint(ridx) = Synergy_stat.NBNBjoint(ridx)/(NBjoint.s.MI+NBjoint.d.MI);
             FracSyn_stat.LOLEMarg(ridx) = Synergy_stat.LOLEMarg(ridx)/(Marg.s.MI+Marg.d.MI);
             FracSyn_stat.NBNBjoint(ridx) = Synergy_stat.NBNBjoint(ridx)/(NB.s.MI+NB.d.MI);
             
             StimSpecSynergy_stat.GL(ridx,:,:) = reshape(GL.StimSpecSynergy,1,ns,nd);
             StimSpecSynergy_stat.PV(ridx,:,:) = reshape(PV.StimSpecSynergy,1,ns,nd);
             StimSpecSynergy_stat.NB(ridx,:,:) = reshape(NB.StimSpecSynergy,1,ns,nd);
             StimSpecSynergy_stat.Marg(ridx,:,:) = reshape(Marg.StimSpecSynergy,1,ns,nd);
             StimSpecSynergy_stat.NBjoint(ridx,:,:) = reshape(NBjoint.StimSpecSynergy,1,ns,nd);
             StimSpecSynergy_stat.LOLE(ridx,:,:) = reshape(LOLE.StimSpecSynergy,1,ns,nd);
             
             %save pop count
             popcount_save(ridx)=nansum(mubar);
             
        end %ridx
        
       % Computes Relevant Means and stds for figures.
        
        Synergy.LOLEMarg(tidx,nidx) = mean(Synergy_stat.LOLEMarg);
        Synergy.NBNBjoint(tidx,nidx) = mean(Synergy_stat.NBNBjoint);
        Synergy_std.LOLEMarg(tidx,nidx) = std(Synergy_stat.LOLEMarg);
        Synergy_std.NBNBjoint(tidx,nidx) = std(Synergy_stat.NBNBjoint);

        FracVars.GL(tidx,nidx) = mean(FracVars_stat.GL);
        FracVars.PV(tidx,nidx) = mean(FracVars_stat.PV);
        FracVars.NB(tidx,nidx) = mean(FracVars_stat.NB);
        FracVars.NBjoint(tidx,nidx) = mean(FracVars_stat.NBjoint);
        FracVars.Marg(tidx,nidx) = mean(FracVars_stat.Marg);
        FracVars.LOLE(tidx,nidx) = mean(FracVars_stat.LOLE);
        
        degVard.GL(tidx,nidx) = mean(degVard_stat.GL);
        degVard.PV(tidx,nidx) = mean(degVard_stat.PV);
        degVard.NB(tidx,nidx) = mean(degVard_stat.NB);
        degVard.NBjoint(tidx,nidx) = mean(degVard_stat.NBjoint);
        degVard.Marg(tidx,nidx) = mean(degVard_stat.Marg);
        degVard.LOLE(tidx,nidx) = mean(degVard_stat.LOLE);
        
        FracVars_std.GL(tidx,nidx) = std(FracVars_stat.GL);
        FracVars_std.PV(tidx,nidx) = std(FracVars_stat.PV);
        FracVars_std.NB(tidx,nidx) = std(FracVars_stat.NB);
        FracVars_std.NBjoint(tidx,nidx) = std(FracVars_stat.NBjoint);
        FracVars_std.Marg(tidx,nidx) = std(FracVars_stat.Marg);
        FracVars_std.LOLE(tidx,nidx) = std(FracVars_stat.LOLE);
        
        degVard_std.GL(tidx,nidx) = std(degVard_stat.GL);
        degVard_std.PV(tidx,nidx) = std(degVard_stat.PV);
        degVard_std.NB(tidx,nidx) = std(degVard_stat.NB);
        degVard_std.NBjoint(tidx,nidx) = std(degVard_stat.NBjoint);
        degVard_std.Marg(tidx,nidx) = std(degVard_stat.Marg);
        degVard_std.LOLE(tidx,nidx) = std(degVard_stat.LOLE);
        
        
        StimSpecSynergy.GL(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.GL);
        StimSpecSynergy.PV(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.PV);
        StimSpecSynergy.NB(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.NB);
        StimSpecSynergy.NBjoint(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.NBjoint);
        StimSpecSynergy.Marg(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.Marg);
        StimSpecSynergy.LOLE(tidx,nidx,:,:) = mean(StimSpecSynergy_stat.LOLE);
        
        StimSpecSynergy_std.GL(tidx,nidx,:,:) = std(StimSpecSynergy_stat.GL);
        StimSpecSynergy_std.PV(tidx,nidx,:,:) = std(StimSpecSynergy_stat.PV);
        StimSpecSynergy_std.NB(tidx,nidx,:,:) = std(StimSpecSynergy_stat.NB);
        StimSpecSynergy_std.NBjoint(tidx,nidx,:,:) = std(StimSpecSynergy_stat.NBjoint);
        StimSpecSynergy_std.Marg(tidx,nidx,:,:) = std(StimSpecSynergy_stat.Marg);
        StimSpecSynergy_std.LOLE(tidx,nidx,:,:) = std(StimSpecSynergy_stat.LOLE);
                
        FracSyn.LOLEMarg(tidx,nidx) = mean(FracSyn_stat.LOLEMarg);
        FracSyn.NBNBjoint(tidx,nidx) = mean(FracSyn_stat.NBNBjoint);
        FracSyn_std.LOLEMarg(tidx,nidx) = std(FracSyn_stat.LOLEMarg);
        FracSyn_std.NBNBjoint(tidx,nidx) = std(FracSyn_stat.NBNBjoint);
        
       
        MIerr.GL(tidx,nidx) = mean(MIerrstat.GL);
        MIserr.GL(tidx,nidx) = mean(MIserrstat.GL);
        MIderr.GL(tidx,nidx) = mean(MIderrstat.GL);
        MIerr.PV(tidx,nidx) = mean(MIerrstat.PV);
        MIserr.PV(tidx,nidx) = mean(MIserrstat.PV);
        MIderr.PV(tidx,nidx) = mean(MIderrstat.PV);
        MIerr.Marg(tidx,nidx) = mean(MIerrstat.Marg);
        MIserr.Marg(tidx,nidx) = mean(MIserrstat.Marg);
        MIderr.Marg(tidx,nidx) = mean(MIderrstat.Marg);
        MIerr.NB(tidx,nidx) = mean(MIerrstat.NB);
        MIserr.NB(tidx,nidx) = mean(MIserrstat.NB);
        MIderr.NB(tidx,nidx) = mean(MIderrstat.NB);
        MIerr.NBjoint(tidx,nidx) = mean(MIerrstat.NBjoint);
        MIserr.NBjoint(tidx,nidx) = mean(MIserrstat.NBjoint);
        MIderr.NBjoint(tidx,nidx) = mean(MIderrstat.NBjoint);
        MIerr.LOLE(tidx,nidx) = mean(MIstat.LOLE);
        MIserr.LOLE(tidx,nidx) = mean(MIs_stat.LOLE);
        MIderr.LOLE(tidx,nidx) = mean(MId_stat.LOLE);

        MIerr_std.GL(tidx,nidx) = std(MIerrstat.GL);
        MIserr_std.GL(tidx,nidx) = std(MIserrstat.GL);
        MIderr_std.GL(tidx,nidx) = std(MIderrstat.GL);
        MIerr_std.PV(tidx,nidx) = std(MIerrstat.PV);
        MIserr_std.PV(tidx,nidx) = std(MIserrstat.PV);
        MIderr_std.PV(tidx,nidx) = std(MIderrstat.PV);
        MIerr_std.Marg(tidx,nidx) = std(MIerrstat.Marg);
        MIserr_std.Marg(tidx,nidx) = std(MIserrstat.Marg);
        MIderr_std.Marg(tidx,nidx) = std(MIderrstat.Marg);
        MIerr_std.NB(tidx,nidx) = std(MIerrstat.NB);
        MIserr_std.NB(tidx,nidx) = std(MIserrstat.NB);
        MIderr_std.NB(tidx,nidx) = std(MIderrstat.NB);
        MIerr_std.NBjoint(tidx,nidx) = std(MIerrstat.NBjoint);
        MIserr_std.NBjoint(tidx,nidx) = std(MIserrstat.NBjoint);
        MIderr_std.NBjoint(tidx,nidx) = std(MIderrstat.NBjoint);
        MIerr_std.LOLE(tidx,nidx) = std(MIstat.LOLE);
        MIserr_std.LOLE(tidx,nidx) = std(MIs_stat.LOLE);
        MIderr_std.LOLE(tidx,nidx) = std(MId_stat.LOLE);

        
        MI.GL(tidx,nidx) = mean(MIstat.GL);
        MI_std.GL(tidx,nidx) = std(MIstat.GL);
        MIs.GL(tidx,nidx) = mean(MIs_stat.GL);
        MIs_std.GL(tidx,nidx) = std(MIs_stat.GL);
        MId.GL(tidx,nidx) = mean(MId_stat.GL);
        MId_std.GL(tidx,nidx) = std(MId_stat.GL);
        FIs.GL(tidx,nidx) = mean(FIs_stat.GL);
        FIs_std.GL(tidx,nidx) = std(FIs_stat.GL);
        FId.GL(tidx,nidx) = mean(FId_stat.GL);
        FId_std.GL(tidx,nidx) = std(FId_stat.GL);
        ERRs.GL(tidx,nidx) = mean(ERRs_stat.GL);
        ERRs_std.GL(tidx,nidx) = std(ERRs_stat.GL);
        ERRd.GL(tidx,nidx) = mean(ERRd_stat.GL);
        ERRd_std.GL(tidx,nidx) = std(ERRd_stat.GL);
        Synergy.GL(tidx,nidx) = mean(Synergy_stat.GL);
        Synergy_std.GL(tidx,nidx) = std(Synergy_stat.GL);
        ErrSynergy.GL(tidx,nidx) = mean(ErrSynergy_stat.GL);
        ErrSynergy_std.GL(tidx,nidx) = std(ErrSynergy_stat.GL);
        SI.GL(tidx,nidx) = mean(SI_stat.GL);
        SI_std.GL(tidx,nidx) = std(SI_stat.GL);
        
        MI.PV(tidx,nidx) = mean(MIstat.PV);
        MI_std.PV(tidx,nidx) = std(MIstat.PV);
        MIs.PV(tidx,nidx) = mean(MIs_stat.PV);
        MIs_std.PV(tidx,nidx) = std(MIs_stat.PV);
        MId.PV(tidx,nidx) = mean(MId_stat.PV);
        MId_std.PV(tidx,nidx) = std(MId_stat.PV);
        FIs.PV(tidx,nidx) = mean(FIs_stat.PV);
        FIs_std.PV(tidx,nidx) = std(FIs_stat.PV);
        FId.PV(tidx,nidx) = mean(FId_stat.PV);
        FId_std.PV(tidx,nidx) = std(FId_stat.PV);
        ERRs.PV(tidx,nidx) = mean(ERRs_stat.PV);
        ERRs_std.PV(tidx,nidx) = std(ERRs_stat.PV);
        ERRd.PV(tidx,nidx) = mean(ERRd_stat.PV);
        ERRd_std.PV(tidx,nidx) = std(ERRd_stat.PV);
        Synergy.PV(tidx,nidx) = mean(Synergy_stat.PV);
        Synergy_std.PV(tidx,nidx) = std(Synergy_stat.PV);
        ErrSynergy.PV(tidx,nidx) = mean(ErrSynergy_stat.PV);
        ErrSynergy_std.PV(tidx,nidx) = std(ErrSynergy_stat.PV);
        SI.PV(tidx,nidx) = mean(SI_stat.PV);
        SI_std.PV(tidx,nidx) = std(SI_stat.PV);
        
        MI.NB(tidx,nidx) = mean(MIstat.NB);
        MI_std.NB(tidx,nidx) = std(MIstat.NB);
        MIs.NB(tidx,nidx) = mean(MIs_stat.NB);
        MIs_std.NB(tidx,nidx) = std(MIs_stat.NB);
        MId.NB(tidx,nidx) = mean(MId_stat.NB);
        MId_std.NB(tidx,nidx) = std(MId_stat.NB);
        FIs.NB(tidx,nidx) = mean(FIs_stat.NB);
        FIs_std.NB(tidx,nidx) = std(FIs_stat.NB);
        FId.NB(tidx,nidx) = mean(FId_stat.NB);
        FId_std.NB(tidx,nidx) = std(FId_stat.NB);
        ERRs.NB(tidx,nidx) = mean(ERRs_stat.NB);
        ERRs_std.NB(tidx,nidx) = std(ERRs_stat.NB);
        ERRd.NB(tidx,nidx) = mean(ERRd_stat.NB);
        ERRd_std.NB(tidx,nidx) = std(ERRd_stat.NB);
        Synergy.NB(tidx,nidx) = mean(Synergy_stat.NB);
        Synergy_std.NB(tidx,nidx) = std(Synergy_stat.NB);
        ErrSynergy.NB(tidx,nidx) = mean(ErrSynergy_stat.NB);
        ErrSynergy_std.NB(tidx,nidx) = std(ErrSynergy_stat.NB);
        SI.NB(tidx,nidx) = mean(SI_stat.NB);
        SI_std.NB(tidx,nidx) = std(SI_stat.NB);
        
        MI.NBjoint(tidx,nidx) = mean(MIstat.NBjoint);
        MI_std.NBjoint(tidx,nidx) = std(MIstat.NBjoint);
        MIs.NBjoint(tidx,nidx) = mean(MIs_stat.NBjoint);
        MIs_std.NBjoint(tidx,nidx) = std(MIs_stat.NBjoint);
        MId.NBjoint(tidx,nidx) = mean(MId_stat.NBjoint);
        MId_std.NBjoint(tidx,nidx) = std(MId_stat.NBjoint);
        FIs.NBjoint(tidx,nidx) = mean(FIs_stat.NBjoint);
        FIs_std.NBjoint(tidx,nidx) = std(FIs_stat.NBjoint);
        FId.NBjoint(tidx,nidx) = mean(FId_stat.NBjoint);
        FId_std.NBjoint(tidx,nidx) = std(FId_stat.NBjoint);
        ERRs.NBjoint(tidx,nidx) = mean(ERRs_stat.NBjoint);
        ERRs_std.NBjoint(tidx,nidx) = std(ERRs_stat.NBjoint);
        ERRd.NBjoint(tidx,nidx) = mean(ERRd_stat.NBjoint);
        ERRd_std.NBjoint(tidx,nidx) = std(ERRd_stat.NBjoint);
        Synergy.NBjoint(tidx,nidx) = mean(Synergy_stat.NBjoint);
        Synergy_std.NBjoint(tidx,nidx) = std(Synergy_stat.NBjoint);
        ErrSynergy.NBjoint(tidx,nidx) = mean(ErrSynergy_stat.NBjoint);
        ErrSynergy_std.NBjoint(tidx,nidx) = std(ErrSynergy_stat.NBjoint);
        SI.NBjoint(tidx,nidx) = mean(SI_stat.NBjoint);
        SI_std.NBjoint(tidx,nidx) = std(SI_stat.NBjoint);
        
        ErrSynergy.NBNBjoint(tidx,nidx) = mean(ErrSynergy_stat.NBNBjoint);
        ErrSynergy_std.NBNBjoint(tidx,nidx) = std(ErrSynergy_stat.NBNBjoint);
        
        MI.Marg(tidx,nidx) = mean(MIstat.Marg);
        MI_std.Marg(tidx,nidx) = std(MIstat.Marg);
        MIs.Marg(tidx,nidx) = mean(MIs_stat.Marg);
        MIs_std.Marg(tidx,nidx) = std(MIs_stat.Marg);
        MId.Marg(tidx,nidx) = mean(MId_stat.Marg);
        MId_std.Marg(tidx,nidx) = std(MId_stat.Marg);
        FIs.Marg(tidx,nidx) = mean(FIs_stat.Marg);
        FIs_std.Marg(tidx,nidx) = std(FIs_stat.Marg);
        FId.Marg(tidx,nidx) = mean(FId_stat.Marg);
        FId_std.Marg(tidx,nidx) = std(FId_stat.Marg);
        ERRs.Marg(tidx,nidx) = mean(ERRs_stat.Marg);
        ERRs_std.Marg(tidx,nidx) = std(ERRs_stat.Marg);
        ERRd.Marg(tidx,nidx) = mean(ERRd_stat.Marg);
        ERRd_std.Marg(tidx,nidx) = std(ERRd_stat.Marg);
        Synergy.Marg(tidx,nidx) = mean(Synergy_stat.Marg);
        Synergy_std.Marg(tidx,nidx) = std(Synergy_stat.Marg);
        ErrSynergy.Marg(tidx,nidx) = mean(ErrSynergy_stat.Marg);
        ErrSynergy_std.Marg(tidx,nidx) = std(ErrSynergy_stat.Marg);
        SI.Marg(tidx,nidx) = mean(SI_stat.Marg);
        SI_std.Marg(tidx,nidx) = std(SI_stat.Marg);
        
        MI.LOLE(tidx,nidx) = mean(MIstat.LOLE);
        MI_std.LOLE(tidx,nidx) = std(MIstat.LOLE);
        MIs.LOLE(tidx,nidx) = mean(MIs_stat.LOLE);
        MIs_std.LOLE(tidx,nidx) = std(MIs_stat.LOLE);
        MId.LOLE(tidx,nidx) = mean(MId_stat.LOLE);
        MId_std.LOLE(tidx,nidx) = std(MId_stat.LOLE);
        FIs.LOLE(tidx,nidx) = mean(FIs_stat.LOLE);
        FIs_std.LOLE(tidx,nidx) = std(FIs_stat.LOLE);
        FId.LOLE(tidx,nidx) = mean(FId_stat.LOLE);
        FId_std.LOLE(tidx,nidx) = std(FId_stat.LOLE);
        ERRs.LOLE(tidx,nidx) = mean(ERRs_stat.LOLE);
        ERRs_std.LOLE(tidx,nidx) = std(ERRs_stat.LOLE);
        ERRd.LOLE(tidx,nidx) = mean(ERRd_stat.LOLE);
        ERRd_std.LOLE(tidx,nidx) = std(ERRd_stat.LOLE);
        Synergy.LOLE(tidx,nidx) = mean(Synergy_stat.LOLE);
        Synergy_std.LOLE(tidx,nidx) = std(Synergy_stat.LOLE);
        ErrSynergy.LOLE(tidx,nidx) = mean(Synergy_stat.LOLE);
        ErrSynergy_std.LOLE(tidx,nidx) = std(Synergy_stat.LOLE);
        SI.LOLE(tidx,nidx) = mean(SI_stat.LOLE);
        SI_std.LOLE(tidx,nidx) = std(SI_stat.LOLE);
        
        ErrSynergy.LOLEMarg(tidx,nidx) = mean(ErrSynergy_stat.LOLEMarg);
        ErrSynergy_std.LOLEMarg(tidx,nidx) = std(ErrSynergy_stat.LOLEMarg);
        
        FracSyn.GL(tidx,nidx) = mean(FracSyn_stat.GL);
        FracSyn.PV(tidx,nidx)=mean(FracSyn_stat.PV);
        FracSyn.NB(tidx,nidx)=mean(FracSyn_stat.NB);
        FracSyn.NBjoint(tidx,nidx)=mean(FracSyn_stat.NBjoint);
        FracSyn.Marg(tidx,nidx)=mean(FracSyn_stat.Marg);
        FracSyn.LOLE(tidx,nidx)=mean(FracSyn_stat.LOLE);

        FracSyn_std.GL(tidx,nidx) = std(FracSyn_stat.GL);
        FracSyn_std.PV(tidx,nidx)=std(FracSyn_stat.PV);
        FracSyn_std.NB(tidx,nidx)=std(FracSyn_stat.NB);
        FracSyn_std.NBjoint(tidx,nidx)=std(FracSyn_stat.NBjoint);
        FracSyn_std.Marg(tidx,nidx)=std(FracSyn_stat.Marg);
        FracSyn_std.LOLE(tidx,nidx)=std(FracSyn_stat.LOLE);

        Popcount(tidx,nidx) = mean(popcount_save);
  toc  
    end %corr (tidx) loop
end %nCells (nidx) loop
mudir=squeeze(mean(mu,2));
musp=squeeze(mean(mu,3));
mumu=mean(mudir,2);

for n=1:nn
for s=1:ns
for d=1:nd
    musep(n,s,d)=mudir(n,d)*musp(n,s)/mumu(n);
end
end
end
muvec=reshape(mu,nn,ns*nd);
musepvec=reshape(musep,nn,ns*nd);
SIemp = 1-mean((muvec'-musepvec').^2)./var(muvec');


for n=1:nn
    temp=squeeze(mu(n,:,:));
    [U,S,V]=svd(temp);
    
SIsvd(n) = S(1,1)^2/sum(diag(S).^2);
end

figure
xvals = log2(nCells_vec);
yvals = [MIerr.GL(tidx,:)./MImax(tidx,:); MIerr.PV(tidx,:)./MImax(tidx,:); MIerr.NB(tidx,:)./MImax(tidx,:); MIerr.NBjoint(tidx,:)./MImax(tidx,:); MIerr.Marg(tidx,:)./MImax(tidx,:); MI.LOLE(tidx,:)./MImax(tidx,:)];
ebars = [MIerr_std.GL(tidx,:)./MImax(tidx,:); MIerr_std.PV(tidx,:)./MImax(tidx,:); MIerr_std.NB(tidx,:)./MImax(tidx,:); MIerr_std.NBjoint(tidx,:)./MImax(tidx,:); MIerr_std.Marg(tidx,:)./MImax(tidx,:); MIerr_std.LOLE(tidx,:)./MImax(:,nidx)];
for idx=1:6
    shadedErrorBar(xvals,yvals(idx,:),ebars(idx,:),'lineProps','o-','transparent',1);
    hold on
end
set(gca,'Box','Off','TickDir','Out','FontSize',20,'TickLength',[.025 , .01]);
legend('GL','PV','NB (1D)','NB (2D)','ML (1D)','ML (2D)');
xlabel('log2 N');
ylabel('MI/MIpurs');
ylim([0.4 1.1]);
xlim([min(log2(nCells_vec)) max(log2(nCells_vec))]);

figure
tidx=5; %100ms
set(gca,'Box','Off','TickDir','Out','FontSize',20,'TickLength',[.025 , .01]);
shadedErrorBar(log2(nCells_vec),ErrSynergy.GL(tidx,:)',ErrSynergy_std.GL(tidx,:)','lineProps','o-','transparent',1);
hold on
shadedErrorBar(log2(nCells_vec),ErrSynergy.PV(tidx,:)',ErrSynergy_std.PV(tidx,:)','lineProps','o-','transparent',1);

shadedErrorBar(log2(nCells_vec),FracSyn.NBNBjoint(tidx,:)',FracSyn_std.NBNBjoint(tidx,:)','lineProps','ko-','transparent',1);

shadedErrorBar(log2(nCells_vec),FracSyn.LOLEMarg(tidx,:)',FracSyn_std.LOLEMarg(tidx,:)','lineProps','ko-','transparent',1);
plot(log2(nCells_vec), 0.19*ones(1,length(nCells_vec)),'k--');

legend('GL','PV','NB','ML','pursuit')
xlabel('log2 N');
ylabel('frac synergy');
ylim([0 0.4]);
xlim([min(log2(nCells_vec)) max(log2(nCells_vec))]);