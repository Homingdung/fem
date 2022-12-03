function U=fem1d(x)
M=length(x);
for i=1:M-1
    h(i)=x(i+1)-x(i);
end
%Initialization of the assembling
A=sparse(M,M);
F=zeros(M,1);
A(1,1) = 1;
F(1)=0;
A(M,M)=1;
F(M)=0;
A(2,2)=1./h(1);
F(2)=simpson_hat1(x(1),x(2));
for i=2:M-2
    A(i,i)=A(i,i)+1/h(i);
    A(i,i+1)=A(i,i+1)-1/h(i);
    A(i+1,i)=A(i+1,i)-1/h(i);
    A(i+1,i+1)=A(i+1,i+1)+1/h(i);
    F(i)=F(i)+simpson_hat2(x(i),x(i+1));
    F(i+1)=F(i+1)+simpson_hat1(x(i),x(i+1));
end

A(M-1,M-1)=A(M-1,M-1)+1/h(M-1);
F(M-1)=F(M-1)+simpson_hat2(x(M-1),x(M));

U=A\F;
return
end


