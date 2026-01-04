

function [mutual_I_inf]=calculate_mutual_info(mydata)

datafrac=[1 .95 .90 0.85 .80 .50];
numfracreps=30;

[condition,trialnum,ts]=size(mydata);

nbin=max(max(max(mydata)))+1;
xx=0:(nbin-1);
mutual_I=zeros(numfracreps,length(datafrac),ts);

for myrp=1:numfracreps
    for mypercent=1:length(datafrac)
        
        clear tpsn;
        tpsn=randperm(trialnum);
        mysn=tpsn(1:floor(datafrac(mypercent)*length(tpsn)));

        P_joint=zeros(condition,nbin,ts);
        
        for tt=1:ts
            for i=1:condition
                clear N;
                N=hist(mydata(i,mysn,tt),xx)./length(mysn);
                P_joint(i,:,tt)=N./condition;
            end
        end

        
        for i=1:ts
            mutual_I(myrp,mypercent,i)=sum(sum(P_joint(:,:,i).*log2(P_joint(:,:,i)./(repmat(sum(P_joint(:,:,i),1),size(P_joint(:,:,i),1),1).*repmat(sum(P_joint(:,:,i),2),1,size(P_joint(:,:,i),2))+eps)+eps)));
        end
    end
end

nstop=length(datafrac)-1;
for i=1:ts
    clear tempdata;
    tempdata=squeeze(mutual_I(:,:,i));
    mutual_I_inf(1:2,i)=polyfit(1./datafrac(1:nstop),squeeze(mean(mutual_I(:,1:nstop,i),1)),1);
    mutual_I_inf(3,i)=std(mutual_I(:,length(datafrac),i),[],1);
end

















