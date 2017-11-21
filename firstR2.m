function [p_angles,iterations,U1Perp]=firstR2(n,k,data,initvecs,INDsize,SBSP,p,swtch)

p_angles=[];       % normalized subspace angles
iterations=[];     % keep track of which iterations we compute the p_angles


U1=orthcomp(initvecs);      % initialize subspace randomly


L=length(data);   


outputindex = linspace(1,L,INDsize);
outputindex = floor(outputindex);


for i=1:(L)
    if any(outputindex==i)
        perc=100*i/(L);                 %output to display the percentage of the algorithm completed
        perc=round(perc);
        string = sprintf('%d percent of algorithm 1R (second version), with switch %d, complete', perc,swtch);
        %disp(string);
        
    end
    
    nu=1/sqrt(i);            % learning rate

    newx=data(i,:)';
    
    %%%%
    if swtch==1
       newx = newx/norm(newx);               %if swtch==1, normalize to the sphere
    elseif swtch==0
    else
        error('Invalid switch')
    end
    %%%
    
    wt=norm(U1'*newx);
    wt=wt^(2-p);
    U1=(U1-(nu*newx*newx'*U1)/(wt));
   
    [U1 S V]=svd(U1);
    size(U1);
    %U1=U1(:,k+1:end);
    U1=U1(:,1:end-k);

    if any(outputindex==i)
        U1Perp=orthcomp(U1);
        
        if ~isequal(size(U1Perp),size(SBSP))
            size(U1Perp)
            error('wrong sizes')
            %pause
        end
        D1=subspace(U1Perp,SBSP);
        
        D1=D1/(pi/2);
        p_angles=[p_angles, D1];
        iterations=[iterations,i];
        
    end

end

U1Perp=orthcomp(U1);
% plot(p_angles)
% pause


return










