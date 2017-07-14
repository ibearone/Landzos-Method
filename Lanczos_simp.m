function [a,b,Matrix_Sub,E_ground,loops]=Lanczos_simp(Matrix_target,tolerance)
%%% parameter
dim=length(Matrix_target) %%%Dimension of Matrix
loops=3
%%%initial parameter
b(1,1)=0
v_0=zeros(dim,1)
v{1}=randi(100,dim,1)
v{1}=v{1}./norm(v{1})
w{1}=Matrix_target*v{1}

%%% Lanczos loops
for k=1
    w{k}=Matrix_target*v{k}-b(k).*v_0
    a(k)=w{k}'*v{k}
    w{k}=w{k}-a(k)*v{k}
    b(k+1)=sqrt(w{k}'*w{k})
    v{k+1}=w{k}./b(k+1)
    k=k+1
end
while k<loops
    w{k}=Matrix_target*v{k}-b(k).*v{k-1}
    a(k)=w{k}'*v{k}
    w{k}=w{k}-a(k)*v{k}
    b(k+1)=sqrt(w{k}'*w{k})
    v{k+1}=w{k}./b(k+1)
    [Matrix_Sub,E]=diagonal_LZSub(a,b)
    E_ground(k)=min(E)
    D_E_ground=abs(E_ground(k)-E_ground(k-1))
    k=k+1
    if D_E_ground>tolerance
        loops=loops+1
    else
        loops=loops
    end
end

end