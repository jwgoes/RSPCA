function diagS1=qre1(diagS,para)
para=1/para;
D=length(diagS);
for i=1:D
r(i)=diagS(i)/sum(min(diagS(i),diagS));
end
[sort_r sort_r_ind]=sort(r,'ascend');
index=find(sort_r>para,1);
if length(index)==0
diagS1=diagS/(sum(diagS));
else
threshold=sum(diagS(sort_r_ind(1:index-1)))/(1/para-(1-index+D));
diagS1=(min(diagS,threshold));
diagS1=diagS1/sum((diagS1));
end
diagS2(sort_r_ind)=diagS1;