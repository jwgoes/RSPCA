function astrodata_final()
% astrodata_final()
%
% Run to re-create experiments on astronomical data



clear all

k=4;
a=fitsread('specarr_robpca_short.fits');
data=fitsread('specarr_robpca_short.fits','Image',1);
data = data(randperm(3485),:);


mean=zeros(1,134);
for i=1:length(data)
    data(i,:)=data(i,:)/median(data(i,:));
    mean = mean + data(i,:);
end
mean = mean/3485;

data=data-repmat(mean,3485,1);


[numpts,n]=size(data);

[initvecs, initvals ] = pca_initialize_random_orthogonal(k,n);

SBSP=initvecs;




p=1;
  d=4;
  re_pca=rerpca(data,d);

 [p_angles1R,iterations,U3]=threeR(n,k,data,initvecs,SBSP,0);
 [p_angles1R,iterations,U3_reg]=three(n,k,data,initvecs,SBSP,0);


 




p = panel();
p.pack(3,4); 
p.de.margin = 3;

for i=1:4
    p(1,i).select()
    %subplot(3,4,i);
    plot(a,re_pca(:,i));
    xlim([min(a), max(a)]);
    ylim([-.5,.5]);
     set(findall(gca, 'Type', 'Line'),'LineWidth',2,'Color','blue');
    if i > 1
        set(gca, 'xtick', [], 'ytick', []);
    else
        set(gca, 'xtick', []);
    end
end
%set(gca,'xaxisLocation','top')
p(1,1).title('\fontsize{15} First eigenspectrum');
p(1,2).title('\fontsize{15} Second eigenspectrum');
p(1,3).title('\fontsize{15} Third eigenspectrum');
p(1,4).title('\fontsize{15} Fourth eigenspectrum');
p(1,1).ylabel('RE-PCA');
%p(1,1).margintop = 1122;
p.margintop = 10;

for i=1:4
    p(2,i).select()
    %subplot(3,4,i+4)
    plot(a,U3_reg(:,i))
    xlim([min(a), max(a)])
    ylim([-.5,.5]);
     set(findall(gca, 'Type', 'Line'),'LineWidth',2,'Color','blue');
    if i > 1
        set(gca, 'xtick', [], 'ytick', []);
    else
        set(gca, 'xtick', []);
    end
end

p(2,1).ylabel('Online PCA');

for i=1:4
    p(3,i).select()
    %subplot(3,4,i+8)
    plot(a,U3(:,i))
     set(findall(gca, 'Type', 'Line'),'LineWidth',2,'Color','blue');
    xlim([min(a), max(a)])
    ylim([-.5,.5]);
    if i > 1
        set(gca, 'xtick', [], 'ytick', []);
    else
        set(gca, 'xtick', []);
    end
end

p(3,1).ylabel('Robust Online PCA');


