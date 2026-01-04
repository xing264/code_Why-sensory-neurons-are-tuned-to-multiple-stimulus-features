function pkn=poissdistrib(rate,k)

pkn=zeros(k+1,1);
%using Poisson distribution
    for kInd=0:k %0:k spikes
        pkn(kInd+1,1)=(rate.^kInd.*exp(-rate))./factorial(kInd);
    end
    %nStims^nDims x k
    pkn=pkn./sum(pkn(:));
end
