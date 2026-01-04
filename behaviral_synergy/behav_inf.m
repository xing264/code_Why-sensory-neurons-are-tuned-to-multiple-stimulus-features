length_trial=500; %-smth_bin+1;

load('ga_allinfodata','I_vec','I_dir','I_spd','Ierr_vec');
I_vec_ga=I_vec;
I_dir_ga=I_dir;
I_spd_ga=I_spd;
Ierr_vec_ga=Ierr_vec;

load('xt_allinfodata','I_vec','I_dir','I_spd','Ierr_vec');
I_vec_xt=I_vec;
I_dir_xt=I_dir;
I_spd_xt=I_spd;
Ierr_vec_xt=Ierr_vec;

load('re_allinfodata','I_vec','I_dir','I_spd','Ierr_vec');
I_vec_re=I_vec;
I_dir_re=I_dir;
I_spd_re=I_spd;
Ierr_vec_re=Ierr_vec;

figure;
errorbar(1:length_trial,(I_vec_ga-I_dir_ga-I_spd_ga)./(I_dir_ga+I_spd_ga),Ierr_vec_ga ,'k')
hold on;
errorbar(1:length_trial,(I_vec_xt-I_dir_xt-I_spd_xt)./(I_dir_xt+I_spd_xt),Ierr_vec_xt ,'r')
hold on;
errorbar(1:length_trial,(I_vec_re-I_dir_re-I_spd_re)./(I_dir_re+I_spd_re),Ierr_vec_re ,'g')
axis square;
box off
set(gca,'TickDir','out')
xlabel('time (ms)')
ylabel('synergy ratio')
ylim([0 1])
xlim([100 300])


% % % % % r
%%

I_vec(I_vec<0)=0;
I_dir(I_dir<0)=0;
I_spd(I_spd<0)=0;

figure;
errorbar(1:length_trial,I_vec,Ierr_vec,'r')
hold on;
errorbar(1:length_trial,I_dir+I_spd,Ierr_dir+Ierr_spd,'c');
hold on;
errorbar(1:length_trial,I_dir,Ierr_dir,'b')
hold on;
errorbar(1:length_trial,I_spd,Ierr_spd,'g')
hold on;
errorbar(1:length_trial,I_dir_shuffle,Ierr_dir_shuffle,'k')
axis square;
box off
set(gca,'TickDir','out')
xlabel('time (ms)')
ylabel('information (bits)')
legend('vec','sum','dir','spd','shuffle')







