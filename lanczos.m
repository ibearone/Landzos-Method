clc;clear;close all;
%%% Setting Parameters%%%
dim=1000 %%%dimention of matrix
loops=100
%%%%%
A=randi(10,dim)
A=(A+transpose(A))/2
b(1)=0
v_0=zeros(dim,1)
v{1}=randi(100,dim,1)
v{1}=v{1}./norm(v{1})

w{1}=A*v{1}
for k=1
    w{k}=A*v{k}-b(k).*v_0
    a(k)=w{k}'*v{k}
    w{k}=w{k}-a(k)*v{k}
    b(k+1)=sqrt(w{k}'*w{k})
    v{k+1}=w{k}./b(k+1)
end
for k=2:loops
    w{k}=A*v{k}-b(k).*v{k-1}
    a(k)=w{k}'*v{k}
    w{k}=w{k}-a(k)*v{k}
    b(k+1)=sqrt(w{k}'*w{k})
    v{k+1}=w{k}./b(k+1)
end
aM=diag(a)
b(1)=[]
b(end)=[]
bM=diag(b)
bM1=[zeros(loops-1,1) bM]
bM1=[bM1;zeros(1,loops)]
bM2=[bM zeros(loops-1,1)]
bM2=[zeros(1,loops);bM2]
H_M=aM+bM1+bM2
plot(a)
hold on;
plot(b)
[V,D,W] = eig(A)

[x_1,y_1,z_1]=eig(H_M)
A_eigen=sort(diag(D))
H_M_eigen=sort(diag(y_1))