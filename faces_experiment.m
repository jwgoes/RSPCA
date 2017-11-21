function faces_experiment()
% faces_experiment()
%
% Run experiment with faces


load data_background
load yale_subsample
data=[J(find(index<2),:);C];  % J(find(index<2),:) are inliers and C are outliers. By index<2 we choose the face from the first person. 
data=[J(find(index<8),:);C];
k=9;%  the dimension of the faces from a same person is at most 9 (there are theory for it..)
swtch =1 ;
p=1;
n=size(data,2);
[U0,S0,V0]=svd(J(find(index<2),:));
SBSP=V0(:,1:k);% this is the true subspace, obtained by PCA on the set of inliers
initvecs= pca_initialize_random_orthogonal(k, size(data,2));

%data=repmat(data,2,1); possibly repeat the inputs
data=data(randperm(size(data,1)),:);% make sure the order of input data points is random
[p_angles{1},iterations,U{1}]=first(n,k,data,initvecs,SBSP,swtch);
note{1}='first';
[p_angles{2},iterations,U{2}]=firstR(n,k,data,initvecs,SBSP,2,swtch);
note{2}='firstR';
[p_angles{3},iterations,U{3}]=firstR2(n,k,data,initvecs,1,SBSP,p,swtch);
U{3}=orth(U{3});
note{3}='firstR2';

initvals = ones( 1, k ) / ( ( n - k ) * ( k) );
[p_angles{4},iterations,U{4}]=second(n,k,data,initvecs,initvals,SBSP,swtch);
note{4}='second';
[p_angles{5},iterations,U{5},S]=secondR(n,k,data,initvecs,initvals,SBSP,swtch);
note{5}='secondR';
[p_angles{6},iterations,U{6},values]=three(n,k,data,initvecs,SBSP,swtch);
note{6}='three';
[p_angles{7},iterations,U{7}]=threeR(n,k,data,initvecs,SBSP,swtch);
note{7}='threeR';
for i=1:7
finalerror(1,i)=p_angles{i}(end); %quantitive measure 1
finalerror(2,i)=norm(U{i}*U{i}'-SBSP*SBSP'); %quantitive measure 2
finalerror(3,i)=norm(U{i}*U{i}'-SBSP*SBSP','fro'); %quantitive measure 3
end


facedata=zeros(20,20*8);facedata2=facedata;facedata3=facedata;
facedata(:,1:20)=reshape(J(1,:),20,20);
facedata2(:,1:20)=reshape(J(50,:),20,20);
facedata3(:,1:20)=reshape(J(20,:),20,20);
for i=1:7
facedata(:,i*20+1:(i+1)*20)=reshape(J(1,:)*U{i}*U{i}',20,20);% project the first inlier to the fitted subspace
facedata2(:,i*20+1:(i+1)*20)=reshape(J(50,:)*U{i}*U{i}',20,20);% project the 50-th inlier to the fitted subspace
facedata3(:,i*20+1:(i+1)*20)=reshape(J(20,:)*U{i}*U{i}',20,20);% project the 50-th inlier to the fitted subspace
end

figure;
imshow(facedata/255);
set(gca,'position',[0 0 1 1]); % plot the images as previous attachment...
figure;
imshow(facedata2/255);
set(gca,'position',[0 0 1 1]);


%save facedata_rep2 facedata facedata2   
