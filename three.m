function [p_angles,iterations,bvecs,values]=three(n,k,data,initvecs,SBSP,swtch)
L = length(data);

d=n-k;
vectors = random( 'normal', 0, 1, n, n- k );
vectors = vectors * pinv( sqrtm( vectors' * vectors ) );
vectors = real( vectors );


values = ones( d,1 ) / ( d);
M = vectors * diag(values) * vectors';

p_angles=[];iterations=[];

for i=1:L
    if mod(i,100)==0
        perc=100*i/(L);
        perc=round(perc);
        string = sprintf('%d percent of algorithm 3, with switch %d, complete', perc,swtch);
        %disp(string);
    end
    eta=1/sqrt(i);
    newx=data(i,:)';
    if swtch==1
        newx=newx/norm(newx);
    elseif swtch==0
    else
        error('Invalid switch')
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Mprev=vectors*diag(values)*vectors';
inexpM = vectors * diag(log(values)) * vectors'- eta*newx*newx' ;







inexpM = .5*(inexpM + inexpM');

[vectors,values] = eig(inexpM);
%[vectors,values,bs] = svd(inexpM);

values = diag(values);
values = exp(values);

if length(values)~=n
    error('values wrong size')
end



values = qre1(values,d);




bvecs = vectors(:,1:k);




if size(SBSP)~=size(bvecs)
    error('subspace is the wrong size')
end
D3=subspace(SBSP,bvecs)/(pi/2);

p_angles=[p_angles D3];
iterations=[iterations i];


end


return










