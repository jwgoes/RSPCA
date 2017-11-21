function [p_angles,iterations,U1]=first(n,k,data,initvecs,SBSP,swtch)

d=n-k;

L=length(data);
p_angles=[];
iterations=[];


U1=initvecs(:,1:k);



L=length(data);

outputindex = linspace(1,L,25);
outputindex = floor(outputindex);

for i=1:L
    if any(outputindex==i)
        perc=100*i/L;
        perc=round(perc);
        %string = sprintf('%d percent of algorithm 1, with switch %d, complete', perc,swtch);
        %disp(string);
        
    end
    nu=1/sqrt(i);

    newx=data(i,:)';
    if swtch==1
        newx/norm(newx);
    elseif swtch==0
    else
        error('Invalid switch')
    end
    
    
    
    U1=U1+(nu*newx*newx'*U1);
    [U1 S V]=svd(U1,0);
    %U1 = U1(:,1:k);
    
%     if ~isequal(size(U1),size(SBSP))
%         size(U1)
%         error('sizes not equal!')
%     end
    D1=subspace(U1,SBSP);
    if D1<eps*100
        D1=0;
    end
    
    
    
    D1=single(D1);
    D1=D1/(pi/2);
    p_angles=[p_angles, D1];
    iterations=[iterations,i];
    
    
end



return










