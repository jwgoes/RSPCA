function [p_angles,iterations,U1]=firstR(n,k,data,initvecs,SBSP,p,swtch)

p_angles=[];       % normalized subspace angles
iterations=[];     % keep track of which iterations we compute the p_angles


U1=initvecs;      % initialize subspace randomly


L=length(data);   

outputindex = linspace(1,L,25);
outputindex = floor(outputindex);
%%
for i=1:(L)
    if any(outputindex==i)
        perc=100*i/(L);                                                        %output to display the percentage of the algorithm completed
        perc=round(perc);
        string = sprintf('%d percent of algorithm 1R (first version), with switch %d, complete', perc,swtch);
        %disp(string);
        
    end
    
    nu=1/sqrt(i);  % learning rate
    newx=data(i,:)';
    
    %%%%
    if swtch==1
       newx = newx/norm(newx);               %if swtch==1, normalize to the sphere
    elseif swtch==0
    else
        error('Invalid switch')
    end
    %%%%
    %%
    
    wt=norm(U1'*newx);
    wt=wt^(2-p);

    U1=U1+(nu*newx*newx'*U1)/(wt);
    
    [U1 S V]=svd(U1,0);
    U1=U1(:,1:k);
 
    
    %%
    if size(U1)~=size(SBSP)
        error('sizes dont match')
    end
    D1=subspace(U1,SBSP);
    D1=D1/(pi/2);
    p_angles=[p_angles, D1];
    iterations=[iterations,i];

end



return








