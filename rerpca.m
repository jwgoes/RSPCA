function [ out ] = rerpca(data,d)
%RERPCA Summary of this function goes here
%   Detailed explanation goes here
%   Data input as NxD
%   rho(t) = 2/pi*arctan(pi/2*t/c^2)
%   rho'(t) = 1/(1+(pi/2*t/c^2)^2)*1/c^2

niter=5;
D=size(data,2);

alpha = 1-1/size(data,1);
c=.787;
delta=.5;

mean_init=1/(size(data,1))*sum(data,1);

initdata=data(1:200,:);
[u,s,v]=svd(initdata);
evect_init=v(:,1:d);
evals_init=diag(s);evals_init=evals_init(1:d);

sumr2=zeros(200,1);
for i=1:200
    sumr2(i)=norm(initdata(i,:)'-evect_init*evect_init'*initdata(i,:)')^2;
end

sig2=mean(sumr2);
hpi=pi/2;
for i=1:d
    sig2=sig2*mean(((1/hpi)*ones(size(sumr2,1),size(sumr2,2))).*atan(hpi*sumr2/sig2))/delta;
end

vk=mean(1./(ones(size(sumr2,1),size(sumr2,2))+(hpi*sumr2/sig2).^2));
qk=mean(sumr2./(ones(size(sumr2,1),size(sumr2,2))+(hpi*sumr2/sig2).^2));
uk=1;




C_prev=zeros(size(data,2));
C_prev=evect_init*diag(evals_init)*evect_init';

% for i=1:200
%     wn=rhoprime(sumr2(i)/sig2);
%    C_prev=C_prev+sig2*wn*(data(i,:)-mean_init)'*(data(i,:)-mean_init);
% end
% den=0;
% for i=1:200
%     den=den+rhoprime(sumr2(i)/sig2)*sumr2(i);
% end
% C_prev=C_prev*sig2/den;

mean_prev=mean_init;
evect_prev=evect_init;
evals_prev=evals_init';

vk_prev=vk;
qk_prev=qk;
uk_prev=uk;

for j=1:1
    data=data(randperm(size(data,1)),:);
for i=1:size(data,1);
    
    y=data(i,:)-mean_prev;
    r=y'-evect_prev*evect_prev'*y';
    sumr2=r'*r;
    
    weight1 = 1 / (1.0+ (hpi * sumr2/sig2)^2);
    weight2 = (1.0/hpi) * atan(hpi * sumr2/sig2) / (sumr2/sig2);
    
    vk = alpha*vk_prev + weight1;          
    qk = alpha*qk_prev + weight1*sumr2;
    uk = alpha*uk_prev + 1;                

    gamma1 = (alpha*vk_prev) / vk   ;
    gamma2 = (alpha*qk_prev) / qk   ;
    gamma3 = (alpha*uk_prev) / uk   ;
    
    vk_prev=vk;
    qk_prev=qk;
    uk_prev=uk; 
   
    mean_prev=gamma1*mean_prev+(1-gamma1)*data(i,:);
    
    sig2 = gamma3*sig2 + (1-gamma3)*weight2*sumr2 / delta;
    
    A=zeros(D,d+1);
    A(:,1:d)=repmat((gamma2*evals_prev).^.5,D,1).*evect_prev;
    A(:,d+1)=(1-gamma2)^.5*y*(sig2/sumr2)^.5;
    
    
    
%   [ui,si,vi]=svd(A);
%    evect_prev=ui(:,1:d);
%    evals_prev=diag(si)';
%    evals_prev=evals_prev(1:d);

     B=A'*A;
     
     [ui,si,vi]=svd(B);
     evect_prev=A*ui';
     evect_prev=evect_prev(:,1:d);
     evals_prev=diag(si)';
     evals_prev=evals_prev(1:d);
   
     C_prev=C_prev*gamma2+(1-gamma2)*sig2*y'*y/sumr2;
     
     [u,v]=eig(C_prev);
     evect_prev=u(:,1:d);
    
end
end
out=evect_prev;



end

function [out] = rho(t)
c=.787;
out = 2/pi*atan(pi/2*t/c^2);
end

function [out] = rhoprime(t)
c=.787;
out = 1/(1+(pi/2*t/c^2)^2)*1/c^2;
end
