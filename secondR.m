function [p_angles,iterations,U,S]=secondR(n,k,data,initvecs,initvals,SBSP,swtch)

iterations=[];
p_angles=[];

M=initvecs*diag(initvals)*initvecs';  %
M=M/trace(M);

size(initvecs);
epsilon=10^(-22);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=length(data);

outputindex = linspace(1,L,25);
outputindex = floor(outputindex);

[U S V] = svd(M);

U=U(:,1:k);S=S(1:k,1:k);


for i=1:L

    newx=data(i,:)';
    
    if any(outputindex==i)
        perc=100*i/(L);
        perc=round(perc);
        string = sprintf('%d percent of algorithm 2R, with switch %d, complete', perc,swtch);
        %disp(string);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if swtch==1
        newx=newx/norm(newx);
    elseif swtch==0
    else
        error('Invalid switch')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    C = U*S*U'; 
    
    
    y=U'*newx;
    z=U*S*y;
    

    R1guy = -(z*z')/(1+newx'*z);
    [x d v] = svd(R1guy);
    x=x(:,1);
    
    xhat = U'*x;
    xperp = x - U*U'*x;
    
    
    
    Q=[[S+xhat*xhat',norm(xperp)*xhat];[norm(xperp)*xhat',norm(xperp)*norm(xperp)]];
    

    
    [Uprime,SS,Vprime]=svd(Q);
    

    UU=[U,xperp/(norm(xperp))]*Uprime;
    UU=UU(:,1:k);
    SS=SS(1:k,1:k);
    
    IA = UU*SS*UU';
    denom = norm(IA*newx);

    
    update=C + newx*newx'/denom;

    [U S V] = svd(update);
    
    U=U(:,1:k);
    S=S(1:k,1:k);
    
    
    




    D2=subspace(U,SBSP)/(pi/2);
    
    if size(U)~=size(SBSP)
        error('sizes dont match')
    end
    
    iterations=[iterations,i];
    p_angles=[p_angles,D2];



    
end


return