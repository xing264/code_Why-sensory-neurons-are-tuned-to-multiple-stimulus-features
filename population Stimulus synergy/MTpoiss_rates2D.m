
function [pop,Sigma,Q,overlap] = MTpoiss_rates2D(nCells,dirs,speeds,Twin,rho,rho_base,rmax,rbase,bw_dir_mean,bw_spd_mean)
%poiss_rates2D  generates firing rates for dir (radians),speed (log2 units) stimuli assuming
%Gaussian tuning functions.  Creates random draw of tuning curves params defined within the function
%Returns the rates, rate derivatives, pref dirs,
%pref speeds, tuning bw etc in the structure 'pop'

pop.nn = nCells;
pop.tWin = Twin;
pop.rho = rho;
pop.rho_base = rho_base;
%stims
nd = length(dirs);
ns = length(speeds);
nn=nCells;
dvals = dirs;
svals = speeds;

%mean cell params
rmax_mean = rmax;  %70; %spk/s %firing rate at preferred direction
            % lower rates will increase variability in the population
            % response but also make the FI calc unstable
%bw_theta_mean = deg2rad(40); %SD of tuning curve, radians
bw_theta_mean = bw_dir_mean;  %deg2rad(50); 
%bw_speed_mean = 1.5; %log2 units
bw_speed_mean = bw_spd_mean; %4;
r_base = rbase;  %5; %mean background firing rate, can set to zero

%tuning curve variability 
rmax_SD = 0.5*rmax_mean;
bw_theta_SD = 0.5*bw_theta_mean; % fraction of mean
bw_speed_SD = 0.5*bw_speed_mean; 
rbase_cells = r_base*(2*rand(1,nn));
rbase_cells(rbase_cells<0)=0; %create variability in bkgrnd rate


%SImin = 0.97;  % SETS MINIMIMUM VALUE OF SEPARABILITY INDEX
%SImin = 0.85;  %note that actual min SI will be much higher
SImin = 0.70;  %03/13/22 trying

nPDirs = nn;  
nPSpds = nn;

%select pref dirs from a distribution biased toward cardinal
         %directions
PDirs=[];
 %tvals = 0:.02:2*pi;
 tvals = 0:.02:2*pi;
 if length(tvals) < nn
     tvals = 0:2*pi/nn:2*pi;
 end
 density = round(10*(2*ones(1,length(tvals)) + 0.5*cos(4*tvals)));
 dsamps=[];
 for ii = 1:length(tvals)
     dsamps = [dsamps ones(1,density(ii))*tvals(ii)];
 end
 dsamps = dsamps(randperm(numel(dsamps)));
 PDirs(1:nPDirs) = dsamps(1:nPDirs); %radians
 
 % PDirs=[1:nPDirs]/nPDirs*2*pi-pi;
 
 PSpds = [];
% temp = .5 + gamrnd(2.7,1.1,1,nn);
 temp = -1.5+gamrnd(4,1.1,1,nn);
 idx=find(temp>7);
 
 temp(idx)=-1+8*rand(1,length(idx));
%  
%  temp=[1:nn]/nn*7-1;
 PSpds(1:nPSpds) = temp(randperm(nn));

 pop.PSpds = PSpds;
 pop.PDirs= PDirs;
 
 minbw_dir = (10*pi/180);
 
 bw_dir = minbw_dir + gamrnd(repmat((bw_theta_mean)^2/bw_theta_SD^2,1,nPDirs),repmat(bw_theta_SD^2/(bw_theta_mean),1,nPDirs));
% bw_dir = repmat(bw_theta_mean,1,nPDirs);
 
 display(['Pop avg. dir tuning SD ' num2str(mean(rad2deg(bw_dir))) ' deg']);
 %seeing some speed bw values that are much too big.  truncating large
 %values
 minbw_speed = 0.5;
 bw_speed = minbw_speed + gamrnd(repmat((bw_speed_mean)^2/bw_speed_SD^2,1,nPSpds),repmat(bw_speed_SD^2/(bw_speed_mean-minbw_speed),1,nPSpds));
 %bw_speed = repmat(bw_speed_mean,1,nPDirs);       
 %bw_speed = gamrnd(repmat(bw_speed_mean^2/bw_speed_SD^2,1,nPSpds),repmat(bw_speed_SD^2/bw_speed_mean,1,nPSpds));

% bw_speed = .5 + gamrnd(repmat(bw_speed_mean^2/bw_speed_SD^2,1,10*nPSpds),repmat(bw_speed_SD^2/bw_speed_mean,1,10*nPSpds));
  %bw_speed = bw_speed(bw_speed<8);
%  bw_speed = bw_speed(1:nPSpds);
 
display(['Pop avg. speed tuning SD ' num2str(mean(bw_speed)) ' log2 units']);
 %bw_speed = 0.5*PSpds;
 rmax_cells = gamrnd(repmat(rmax_mean^2/rmax_SD^2,1,nn),repmat(rmax_SD^2/rmax_mean,1,nn));
 %rmax_cells = repmat(rmax_mean,1,nn);
 
 rhoSI = sqrt(1-SImin)*(2*rand(1,nn)-1);
% rhoSI = 0.1*ones(1,nn);
 
 pop.bw_dir = bw_dir;
 pop.bw_spd = bw_speed;
 pop.rmax = rmax_cells;
 pop.rback = rbase_cells;
 pop.rhoSI = rhoSI;
 
 Sig_nn = zeros(nn);
 overlap.spd=zeros(nn);
 overlap.dir=zeros(nn);
 for i=1:nn
     for j=i:nn
         if i==j
             Sig_nn(i,j)=1;
             overlap.dir(i,j)=0;
             overlap.spd(i,j)=0;
         else
             bwd = mean([bw_dir(i) bw_dir(j)]);
             bws = mean([bw_speed(i) bw_speed(j)]);
             overlap.dir(i,j)= abs(PDirs(i)-PDirs(j))/bwd;
             overlap.spd(i,j) = abs(PSpds(i)-PSpds(j))/bws;
%             Sig_nn(i,j)=rho_base + (rho-rho_base)*exp(-((PSpds(i)-PSpds(j))/bws)^2/2)...
%                 *exp((cos(PDirs(i)-PDirs(j))-1)/bwd^2/2);
             Sig_nn(i,j)=rho_base + (rho-rho_base)*exp(-abs((PSpds(i)-PSpds(j))/bws).^2/2)...
                 *exp((cos(PDirs(i)-PDirs(j))-1)/bwd^2);
%              Sig_nn(i,j)=rho_base + rand(1)*(rho-rho_base)*exp(-abs((PSpds(i)-PSpds(j))/bws))...
%                  *exp(-sqrt(1-cos(PDirs(i)-PDirs(j)))/bwd);

             Sig_nn(j,i)=Sig_nn(i,j);  %initial correlation matrix
             overlap.dir(j,i)=overlap.dir(i,j);
             overlap.spd(j,i)=overlap.spd(i,j);
             
             %uniform pairwise correlation matrix
             %   rho = 0.1;
             %   Sig_nn = ones(nCells)*rho;
             %   Sig_nn(eye(nCells)~=0)=1;

         end
     end
 end
 
  % prevent singularities in FI calc
  
[V,D]=eig(Sig_nn);
e=diag(D);
 if(min(e)<0)
    display([ 'Sig_nn not positive definite: min(e) = ' num2str(min(e))]);

    D(D<0)=0;
    temp=V*D*V' + eye(nn)*0.001;
    diagtemp=1./sqrt(diag(temp));
    
    temp = diag(diagtemp)*temp*diag(diagtemp);
    
     if ~issymmetric(temp)
         temp = (temp + temp')/2;
     end
     
     display([ 'Sig_nn correction has RMS error: ' num2str(sqrt(mean((temp(:)-Sig_nn(:)).^2)))]);
    display([ '                  and max abs error: ' num2str(sqrt(max((temp(:)-Sig_nn(:)).^2)))]);
    %scatter(Sig_nn(:),temp(:))
    drawnow
     Sig_nn = temp;
 end
  Q = Sig_nn;  %save the correlation matrix for sample generation
  
 Sigma = repmat(Sig_nn,1,1,ns,nd); 
 
 rate_fn = zeros(nn,ns,nd);
 ang_diffs=zeros(1,nn);
 spd_diffs=zeros(1,nn);
 drateds = zeros(1,nn);
 dratedd = zeros(1,nn);
 mu = rate_fn;
 dmudd = mu;
 dmuds = mu;
         
 
 
        
 for didx = 1:nd
     for sidx = 1:ns
              
     ang_diffs = circ_dist(repmat(dvals(didx),1,nn),PDirs); %radians
     spd_diffs = repmat(svals(sidx),1,nn) - PSpds; %log2 speed

     rate_fn(1:nn,sidx,didx) = ...%
          rmax_cells.*exp((cos(dvals(didx)-PDirs)-1)./bw_dir.^2./(1-rhoSI.^2))...
         .*exp(-0.5*(spd_diffs./bw_speed).^2./(1-rhoSI.^2))...
         .*exp(-spd_diffs.*sin(dvals(didx)-PDirs)./bw_dir./bw_speed.*rhoSI./(1-rhoSI.^2));
     
     drateds = (-1./bw_speed.^2.*spd_diffs./(1-rhoSI.^2)-sin(dvals(didx)-PDirs)./bw_dir./bw_speed.*rhoSI./(1-rhoSI.^2))...
         .*rmax_cells.*exp((cos(dvals(didx)-PDirs)-1)./bw_dir.^2./(1-rhoSI.^2))...
         .*exp(-0.5*(spd_diffs./bw_speed).^2./(1-rhoSI.^2))...
         .*exp(-spd_diffs.*sin(dvals(didx)-PDirs)./bw_dir./bw_speed.*rhoSI./(1-rhoSI.^2));
     
     dratedd = (-sin(dvals(didx)-PDirs)./bw_dir.^2./(1-rhoSI.^2)-spd_diffs.*cos(dvals(didx)-PDirs)./bw_dir./bw_speed.*rhoSI./(1-rhoSI.^2))...
         .*rmax_cells.*exp((cos(dvals(didx)-PDirs)-1)./bw_dir.^2./(1-rhoSI.^2)) ...
         .*exp(-0.5*(spd_diffs./bw_speed).^2./(1-rhoSI.^2))...
         .*exp(-spd_diffs.*sin(dvals(didx)-PDirs)./bw_dir./bw_speed.*rhoSI./(1-rhoSI.^2));
      
     % mu(1:nn,sidx,didx) = round(rate_fn(1:nn,sidx,didx)*Twin); 
       mu(1:nn,sidx,didx) = rate_fn(1:nn,sidx,didx)*Twin + Twin*rbase_cells';   %trial avg count can be fractional 
            
      dmuds(1:nn,sidx,didx) =  drateds*Twin;
      dmudd(1:nn,sidx,didx) =  dratedd*Twin;
      
     %Turn Correlations matrix Sigma into Variance under marginal poisson assumption
     Sigma(:,:,sidx,didx) = diag(sqrt(mu(:,sidx,didx)))*Sigma(:,:,sidx,didx)*diag(sqrt(mu(:,sidx,didx)));
     Sigma(:,:,sidx,didx) = (Sigma(:,:,sidx,didx) + Sigma(:,:,sidx,didx)')/2;
     
     
 
     %%% Alternate Formula (simpler less likely to be singular)
%      sigmagsq = (rho-rho_base)/(1-rho+rho_base)/rmax/Twin;
%      
%      Sigma(:,:,sidx,didx) = sigmagsq*mu(:,sidx,didx)*mu(:,sidx,didx)'+diag(mu(:,sidx,didx)) ...
%                           + rho_base*sqrt(mu(:,sidx,didx)+sigmagsq*mu(:,sidx,didx).^2)*sqrt(mu(:,sidx,didx)+sigmagsq*mu(:,sidx,didx).^2)';
%                        
     if(min(eig(Sigma(:,:,sidx,didx)))<0)
        'still having min eig problems' 
     end
     %%% END BIG CHANGE

%       e = eig(Sigma(:,:,sidx,didx));
%     if min(e)<0
%          display([ 'Cov min(e) = ' num2str(min(e)) ': spd ' num2str(sidx) ', dir' num2str(didx)]);
%          Sigma(:,:,sidx,didx) = (Sigma(:,:,sidx,didx) - 2*min(e)*eye(nn))/(1-2*min(e));
%          Sigma(:,:,sidx,didx) = (Sigma(:,:,sidx,didx) + Sigma(:,:,sidx,didx)')/2;
%     end
     
     end %sidx
 end %didx
 
 
 
 %compute SI or each neuron
     for nndx = 1:nn
        Stemp = svd(reshape(mu(nndx,:,:),ns,nd));
        SI_nn(nndx) = Stemp(1).^2/sum(Stemp.^2);
     end
     
 %pop.rate = rate_fn;
 %pop.dratedd = dratedd;
 %pop.drateds = drateds;
 pop.mu = mu;
 pop.dmudd = dmudd;
 pop.dmuds = dmuds;
 pop.corr=Sig_nn;
 pop.SIvals = SI_nn;
 
end








