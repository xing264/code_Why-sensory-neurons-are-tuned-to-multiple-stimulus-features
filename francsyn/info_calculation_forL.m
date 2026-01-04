

load neuron_data

orth_dir_cumI_inf=cell(1,30);
orth_spd_cumI_inf=cell(1,30);
orth_vec_cumI_inf=cell(1,30);
h = waitbar(0,'Please wait...');
for i=1:30
    
    [nts,ndir,nspd,ntri]=size(neuron_data{i});
    datause=permute(neuron_data{i},[2,3,4,1]);
    
    clear tpdirdata ;
    for j=1:ndir
        tpdirdata(j,:,:)=reshape(squeeze(datause(j,:,:,:)),nspd*ntri,nts);
    end
    
    clear tpspddata ;
    for j=1:nspd
        tpspddata(j,:,:)=reshape(squeeze(datause(:,j,:,:)),ndir*ntri,nts);
    end
    
    clear tpvecdata ;
    tpvecdata=reshape(datause,ndir*nspd,ntri,nts);
    
        cumdirdata=cumsum(tpdirdata,3);
        cumspddata=cumsum(tpspddata,3);
        cumvecdata=cumsum(tpvecdata,3);
    
    orth_dir_cumI_inf{i}=calculate_mutual_info(cumdirdata);
    orth_spd_cumI_inf{i}=calculate_mutual_info(cumspddata);
    orth_vec_cumI_inf{i}=calculate_mutual_info(cumvecdata);
    
    waitbar(i/30,h)
end




%%


Info_vec=[];
Info_dir=[];
Info_spd=[];

Info_vec_slp=[];
Info_dir_slp=[];
Info_spd_slp=[];

fsc_vec=[];
fsc_dir=[];
fsc_spd=[];

syn_ratio=[];

for i=1:30
    Info_vec(:,i)=orth_vec_cumI_inf{i}(2,:)';
    Info_dir(:,i)=orth_dir_cumI_inf{i}(2,:)';
    Info_spd(:,i)=orth_spd_cumI_inf{i}(2,:)';
    
    Info_vec_slp(:,i)=orth_vec_cumI_inf{i}(1,:)';
    Info_dir_slp(:,i)=orth_dir_cumI_inf{i}(1,:)';
    Info_spd_slp(:,i)=orth_spd_cumI_inf{i}(1,:)';
    
    fsc_vec(:,i)=orth_vec_cumI_inf{i}(1,:)./(orth_vec_cumI_inf{i}(1,:)+orth_vec_cumI_inf{i}(2,:));
    fsc_dir(:,i)=orth_dir_cumI_inf{i}(1,:)./(orth_dir_cumI_inf{i}(1,:)+orth_dir_cumI_inf{i}(2,:));
    fsc_spd(:,i)=orth_spd_cumI_inf{i}(1,:)./(orth_spd_cumI_inf{i}(1,:)+orth_spd_cumI_inf{i}(2,:));
end

% Info_vec(:,15)=Info_vec(:,15)-min(Info_vec(:,15));
% Info_dir(:,15)=Info_dir(:,15)-min(Info_dir(:,15));
% Info_spd(:,15)=Info_spd(:,15)-min(Info_spd(:,15));


for i=1:30
    syn_ratio1(:,i)=(Info_vec(:,i)-Info_dir(:,i)-Info_spd(:,i))./(Info_dir(:,i)+Info_spd(:,i)+eps);
    syn_ratio2(:,i)=(Info_vec(:,i)-Info_dir(:,i)-Info_spd(:,i))./(Info_vec(:,i)+eps);
end

stt_t1=1;
figure;
plot(stt_t1:200,syn_ratio1(stt_t1:end,:),'k')
hold on;
errorbar(stt_t1:200,mean(syn_ratio1(stt_t1:end,:),2),std(syn_ratio1(stt_t1:end,:),[],2),'r');
axis([0 200 -0.2 1])
axis square;
box off;
set(gca,'TickDir','out')
xlabel('time (ms)')
ylabel('fractional synergy')

mean(syn_ratio1(end,:))
std(syn_ratio1(end,:))

figure;
scatter(Info_dir(200,:)+Info_spd(200,:),Info_vec(200,:),'k')
hold on;plot(0:0.1:2.5,0:0.1:2.5,'k')

[p,s]=polyfit(Info_dir(200,:)+Info_spd(200,:),Info_vec(200,:),1)
[cfit,gds]=fit([Info_dir(200,:)+Info_spd(200,:)]',Info_vec(200,:)','poly1')
hold on;
plot(cfit,'r')
axis square
axis([0 2.5 0 2.5])
box off;
set(gca,'TickDir','out')
xlabel('Idir+Ispd (bits)')
ylabel('Ivec (bits)')



%%


load('toyvscell_fullresponse.mat')

pos_info_ratio=[];
I_vec_pos=[];
I_sum_pos=[];

spk_info_ratio=[];
spk_info_ratio2=[];

for i=1:29
    pos_info_ratio(:,i)=toyvscell{i,21}(1:200)./(toyvscell{i,18}(1:200)+toyvscell{i,19}(1:200));
    I_vec_pos(:,i)=toyvscell{i,20}(1:200);
    I_sum_pos(:,i)=toyvscell{i,18}(1:200)+toyvscell{i,19}(1:200);
    
    spk_info_ratio(:,i)=toyvscell{i,11}(1:200)./(toyvscell{i,8}(1:200)+toyvscell{i,9}(1:200));
    
    spk_info_ratio2(:,i)=toyvscell{i,12}(1:200);
end


pos_info_ratio(pos_info_ratio<-1)=-1;
pos_info_ratio(pos_info_ratio>5)=5;

stt_t=1;
figure;
plot(stt_t:200,pos_info_ratio(stt_t:end,:),'k');
hold on;
errorbar(stt_t:200,mean(pos_info_ratio(stt_t:end,:),2),std(pos_info_ratio(stt_t:end,:),[],2),'r')
axis square;
set(gca,'TickDir','out')
box off;
xlabel('time (ms)')
ylabel('fractional synergy')
axis([0 200 -0.5 2.5])



figure;
scatter(I_sum_pos(end,:),I_vec_pos(end,:),'k');
[ft,gf]=fit(I_sum_pos(end,:)',I_vec_pos(end,:)','poly1');
hold on;
plot(0:0.1:2,0:0.1:2,'k')
hold on;
xx=0:0.1:2;
plot(ft,'r')
axis square;
axis([0 2 0 2])
set(gca,'TickDir','out')
box off;
xlabel('Idir+Ispd (bits)')
ylabel('Ivec (bits)')

[cf,gd]=corrcoef(I_sum_pos(end,:),I_vec_pos(end,:))


mean(pos_info_ratio(end,:))
std(pos_info_ratio(end,:))














