function heatmap_gen(D,d,iter,rho_ins,rho_outs)
%
%function heatmap_gen(D,d,iter,rho_ins,rho_outs)
%
%D - ambient dimension
%d - subspace dimension
%iter - number of iterations for each square
%rho_ins - vector of rho_in values to loop over
%rho_outs - vector of rho_out values to loop over
%
%Example: heatmap_gen(100,3,10,[2:10],[2:8])



results_mat = nan(length(rho_ins), length(rho_outs), iter);

for rin = 1:length(rho_ins)
    in = rho_ins(rin);
    for rout = 1:length(rho_outs)
        
        out = rho_outs(rout);
        
        N_in = in * d;
        N_out = out * D;
        
        for i = 1:iter
            rin_rout_i = [rin,rout,i]
            
            %Generate data
            
            cov=[1:D];
            cov=cov/sum(cov);
            
            [data,SBSP] = datagen2(D,d,N_in,N_out,cov,10);
            [initvecs, initvals ] = pca_initialize_random_orthogonal(d,D);
            
            [p_angles,I]=threeR(D,d,data,initvecs,SBSP,1);
            
            
            
            
            
            results_mat(rin,end-rout+1,i) =  p_angles(end);
            
            
        end
    end
end

final_mat = mean(results_mat,3)



imshow(final_mat,'XData',rho_ins,'YData', flip(rho_outs),'InitialMagnification', 1800)
set(gca,'YDir','normal')
axis on;
xlabel('$\rho_{in}$','Interpreter','latex','FontSize',22)
ylabel('$\rho_{out}$','Interpreter','latex','FontSize',22)
