function [Matrix_Sub,E]=diagonal_LZSub(a,b)
%%% Constructing Subspace Matrix
dim_sub=length(a)
aM=diag(a)
b(1)=[]
b(end)=[]
bM=diag(b)
bM1=[zeros(dim_sub-1,1) bM]
bM1=[bM1;zeros(1,dim_sub)]
bM2=[bM zeros(dim_sub-1,1)]
bM2=[zeros(1,dim_sub);bM2]
Matrix_Sub=aM+bM1+bM2

%%% Eignecalue
E=eig(Matrix_Sub)
end