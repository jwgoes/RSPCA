function [data,SBSP] =  datagen2(n,k,nin,nout,CV,TM)


%%%%%%%%% outliers

dCV=diag(CV);
mdCV=mean(dCV);
CV=TM*(1/n)*CV/mdCV;
CV_sqrt = sqrtm(CV);
outs=CV_sqrt*random('normal',0,1,n,nout);




%%%%%%%%% inliers

I=randperm(n);
CVin=zeros(n,1);
CVin(I(1:k))=1/k;
CVin=diag(CVin);
SBSP=CVin(:,I(1:k));
CVin_sqrt = sqrt(CVin);
ins=CVin_sqrt*random('normal',0,1,n,nin);

data=[outs ins]';
data=data(randperm(length(data)),:);
