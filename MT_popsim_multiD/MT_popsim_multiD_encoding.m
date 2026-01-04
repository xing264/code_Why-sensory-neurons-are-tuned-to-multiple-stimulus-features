%MTpop_sim

clear 
close all

% if ispc
% rt='Z:/';
% elseif ismac
% rt='/Volumes/NewHub/';
% else
% rt='/NewHub/';
% end
% to run on linux machine: 
% matlab -r -nodisplay -nojvm
% cd /home/osbornelab/Desktop/popsim_combinatorial

experiment='MTpopsim_multiD';

nDims_use=1:7; %6-8 take a while: use up to 5 for basic checks
%nDims_use=1:4;

%change scale factors to look at effects of changing bandwidth and preferred stim
%uncomment 
scale_factors=[0.04 0.08 0.16 0.25 0.33 0.42 0.5]; %scale factors for bw_stim, and set prefstim below to 0
% scale_factors=[0:0.5:3]; %scale factors for prefstim, and set bw_stim below to 1
% scale_factors=1;
cols=colormap(jet(length(scale_factors)));

%interesting tradeoff:
%low BW = higher info/fractional synergy for low nDims, but high BW = higher info/FS for
%high nDims: because with high nDims and narrow tuning, higher likelihood
%of no spike rate change in a dim when other dims are at lowest point of
%TC, so wide tuning allows for better S/R relationship at range of stim
%values

%%
for fact_ind=1:length(scale_factors)

    fact_ind
    %reps to look at effect of randomly changing bandwidth or preferred stim
    repnum=1;
    scale_factor=scale_factors(fact_ind);
    bwvsinfo=zeros(length(nDims_use),repnum,2);
    
%for rep=1:repnum
    %rep
    %half bandwidth (ie SD)
    %uncomment if changing scale_factor to change bandwidth
    bw_stim = ones(max(nDims_use),1)*scale_factor; 
    %uncomment if not changing bandwidth
%     bw_stim = ones(max(nDims_use),1); 
    %randomly choose within range from 0 to 1, arbitrary units, random distribution
%     bw_stim=abs(randn(max(nDims_use),1))*scale_factor; 
%     bw_stim=sort(bw_stim,'descend');
    %if changing scale_factor to change preferred stim
%     prefstim=ones(max(nDims_use),1)*scale_factor;
    %if not changing preferred stim
    prefstim=ones(max(nDims_use),1)*0;
    %randomly choose within range from 0 to 1, arbitrary units, uniform distribution
%     prefstim=rand(max(nDims_use),1); 


for nDims=nDims_use
    clearvars -except experiment Sy Sy_given_x Sy_jt Sy_given_x_jt nDims rt nDims_use bw_stim repnum rep bwvsinfo ...
        prefstim cols scale_factor scale_factors fact_ind Ixy_1d Ixy_sum Ixy_multid meanct varct gaussent expent pjointsigma pjointsigma_count

    nDims
    tic
    %stimuli ranging from prefdir to 3*bandwidth
    %changing BW and stims together doesn't change info or fractional synergy,
    %so changing BW here is just changing the ratio between the stimulus range
    %and the dynamic range of the cell

    stims=-0.5:0.1:0.5;
    stims=reshape(stims,length(stims),1);
    nStimsx=length(stims);
    nStims = nStimsx^(nDims);

    %Create stim vectors
    %Matt's orig code used the same stim value for every dimension, instead of
    %all combinations of stim values.  Create an array of stimulus value
    %indices that gives all combinations
    sidx = zeros(nDims,nStims);
    for idx1=1:nDims
        temp=[];
        if idx1==1
            temp = repmat([1:nStimsx],1,nStimsx^(nDims-1));
            sidx(1,1:nStims) = temp;
        else
            for idx2=1:nStimsx
                temp = [temp repmat(idx2,1,nStimsx^(idx1-1))];
            end
            sidx(idx1,1:nStims) = repmat(temp,1,nStimsx^(nDims-idx1));
        end
    end %idx1
        
    
    %neurons

    rmax = 100; %spk/s %how does changing this affect info too
    k=40; %max spike count tested. MM used k=20 for rmax=(something small)
    %for PSTH
    % load('all_spk_tp_allcells.mat')
    % rmax_temp=cumsum(nanmean(nanmean(nanmean(all_spk_tp_allcells(:,6,:,:,:),4),5),3));

    times=0.2;
    % times=scale_factor;
    r_base = 0;
    bw_stim=bw_stim(1); %assume all BW are same

    %for t=1:length(times)
    Twin=times; %(t);
    stimdiff = stims - prefstim(1); %length nStimsx, sind(nDims,nStims)
    %note that all pref stims are the same
    
    if nDims==1
        rate_by_n_sqrd = stimdiff.^2;
    else
        sind_by_n = reshape(sidx,1,numel(sidx));
        temp = stimdiff(sidx);
        temp = reshape(temp,nDims,nStims);
        rate_by_n_sqrd = sum(temp.^2,1); %sum over dimensions to get the arg of the exp
    end
    rate_prod = r_base + rmax*Twin*exp(-0.5.*rate_by_n_sqrd(:)./bw_stim.^2);
    meanct(fact_ind,nDims)=mean(rate_prod(:));
    varct(fact_ind,nDims) = var(rate_prod(:));
       
    pjoint_all=arrayfun(@poissdistrib,rate_prod,ones(size(rate_prod))*k,'UniformOutput',0);
    pjoint_all=reshape([pjoint_all{:}],[k+1,size(rate_prod)]);
    pjoint_all=shiftdim(pjoint_all,1);
        
    pja_sum=sum(pjoint_all(:));
    Pjoint_stim = pjoint_all./pja_sum;
    nBins_x= nStims;
    nBins_y = k+1;
    Pcount = sum(Pjoint_stim,1);
    Pcount = Pcount./sum(Pcount);
    temp = sum(sum(Pjoint_stim.*log(Pjoint_stim./(sum(Pjoint_stim,1).*sum(Pjoint_stim,2)+eps) +eps)));
    Ixy_multid(fact_ind,nDims) = temp;
    

%loop over each dimension and then sum to get Ixy_stim
%We also need the joint distribution of counts and stim values for each
    %dimension individually
    %rate_by_dim = zeros(nStims,k+1);
    %nBins_x= nStimsx;
    nBins_y = k+1;
    Itemp=[];
    if nDims==1
        Ixy_sum(fact_ind,nDims) = Ixy_multid(fact_ind,nDims);
    else
         %this calculation should use the rate given by all dims
      % but the stim value only for the dim under consideration
       for dind=1:nDims 
           pjoint_by_dim=zeros(nStimsx,k+1);
           for sind=1:nStims
               ct = round(rate_prod(sind))+1;
               s_int = sidx(dind,sind);
               tempdata = poisspdf([0:k],rate_prod(sind));
               pjoint_by_dim(s_int,:) = pjoint_by_dim(s_int,:) + tempdata;
               %pjoint_by_dim(s_int,ct) = pjoint_by_dim(s_int,ct) +1/nStims;
            end
            pjoint_by_dim = pjoint_by_dim./sum(sum(pjoint_by_dim));
            Itemp(dind) = sum(sum(pjoint_by_dim.*log(pjoint_by_dim./(sum(pjoint_by_dim,1).*sum(pjoint_by_dim,2)+eps) +eps)));     
        end
        Ixy_sum(fact_ind,nDims) = sum(Itemp);
    end %if
    

toc
end %nDims



figure(1);hold all
% errorbar(nDims_use,mean(Ixy_stim,2),std(Ixy_stim,[],2),'x-','Color',cols(fact_ind,:),'LineWidth',2);hold all
h1=plot(nDims_use,squeeze(Ixy_multid(fact_ind,:)),'o-','Color',cols(fact_ind,:),'LineWidth',2);hold all
h1.MarkerSize=10;
h1.MarkerFaceColor=h1.Color;
xlim([0 max(nDims_use)+1]);
set(gca,'FontSize',18,'Box','Off','TickDir','Out');
ylabel('information multid (bits)');
xlabel('dimensions');

figure(6);hold all
% errorbar(nDims_use,mean(Ixy_stim,2),std(Ixy_stim,[],2),'x-','Color',cols(fact_ind,:),'LineWidth',2);hold all
h6=plot(nDims_use,squeeze(Ixy_multid(fact_ind,:)),'o-','Color',cols(fact_ind,:),'LineWidth',2);hold all
h6.MarkerSize=10;
h6.MarkerFaceColor=h6.Color;
xlim([0 max(nDims_use)+1]);
set(gca,'FontSize',18,'Box','Off','TickDir','Out');
ylabel('information sum (bits)')
xlabel('dimensions')

figure(2)
hold all;
fracsyn=(Ixy_multid(fact_ind,:)-Ixy_sum(fact_ind,:))./Ixy_sum(fact_ind,:);
fracsyn=squeeze(fracsyn);%LCO Ixy_1d must be the multiD?
%h4=errorbar(nDims_use,mean(fracsyn,2),std(fracsyn,[],2),'o-','Color',cols(fact_ind,:),'LineWidth',2,'MarkerSize',10);hold all
h2 = plot(nDims_use(:),fracsyn(:),'o-','Color',cols(fact_ind,:),'LineWidth',2,'MarkerSize',10);
hold all;
h2.MarkerSize=10;
h2.MarkerFaceColor=h2.Color;
xlim([0 max(nDims_use)+1]);
set(gca,'FontSize',18,'Box','Off','TickDir','Out');
ylabel('fractional synergy');
xlabel('dimensions');


% 
figure(3);
h3=plot(nDims_use,meanct(fact_ind,:),'o-','LineWidth',2,'Color',cols(fact_ind,:));hold all
h3.MarkerSize=10;
h3.MarkerFaceColor=h3.Color;
set(gca,'FontSize',18,'Box','Off','TickDir','Out');
ylabel('mean count');
xlabel('dimensions');
xlim([0 max(nDims_use)+1]);
pause(0.1)

figure(4);
h4=plot(nDims_use,varct(fact_ind,:),'o-','LineWidth',2,'Color',cols(fact_ind,:));hold all
h4.MarkerSize=10;
h4.MarkerFaceColor=h4.Color;
ylabel('var count');
xlabel('dimensions');
xlim([0 max(nDims_use)+1]);
set(gca,'FontSize',18,'Box','Off','TickDir','Out');
pause(0.1)

end %fact_ind
