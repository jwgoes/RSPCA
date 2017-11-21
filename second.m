function [p_angles,iterations,U]=second(n,k,data,initvecs,initvals,SBSP,swtch)

L=length(data);
p_angles=[];
iterations=[];

U=initvecs(:,1:k);
S=initvals(1:k);S=diag(S);

for i=1:L
    if mod(i,100)==0
        perc=100*i/L;
        perc=round(perc);
        string = sprintf('%d percent of algorithm 2, with switch %d, complete', perc,swtch);
        %disp(string);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    newx=data(i,:)';
    if swtch==1
        newx/norm(newx);
    elseif swtch==0
    else
        error('Invalid switch')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
    

    
    xhat=U'*newx;
    xperp=newx-U*U'*newx;
    Q=[[S+xhat*xhat',norm(xperp)*xhat];[norm(xperp)*xhat',norm(xperp)*norm(xperp)]];
    
    [Uprime,S,Vprime]=svd(Q);
    II=[1:n];
    
    S=S(1:k,1:k);
    

    U=[U,xperp/(norm(xperp))]*Uprime;
    U=U(:,1:k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    D2=subspace(U,SBSP)/(pi/2);  
    p_angles=[p_angles, D2];
    iterations=[iterations,i];
    
end





return