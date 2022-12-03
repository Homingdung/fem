%Problem statement 
%Boundary value problem -u''=f, x(0)=x(1)=0;


%Mesh generation
x0=0;
xN=1;
N=10;
h=(xN-x0)/N;
x=x0:h:xN;
x_m=x0:h/2:xN;

%Basis function: qudratic interpolation and derivatives
phi0=@(x)2.*(x-1).*(x-1/2);
phi1=@(x)4.*x.*(1-x);
phi2=@(x)2.*x.*(x-1/2);

dphi0 = @(x)4.*x-3;
dphi1 = @(x)4-8.*x;
dphi2 = @(x)4.*x-1;
%function f
f=@(x)1;
M=2.*N+1;
A=sparse(M,M);
F=zeros(M,1);
%Contribution from a paticular element
K=sparse([7,-8,1;-8,16,-8;1,-8,7]).*(1/(3.*h));
%Assemble stiffness matrix A
for i = 1:N
    A(2.*i-1:2.*i+1,2.*i-1:2.*i+1)=...
        A(2.*i-1:2.*i+1,2.*i-1:2.*i+1)+K;
end
A = A(2:end-1,2:end-1); 
%Assemble load vector F

for i = 1:N
    f0 = @(m) f(h.*m + x(i)).*phi0(m).*h;
    f1 = @(m) f(h.*m + x(i)).*phi1(m).*h;
    f2 = @(m) f(h.*m + x(i)).*phi2(m).*h;
    b0 = integral(f0,0,1);
    b1 = integral(f1,0,1);
    b2 = integral(f2,0,1);
    F(2*i-1:2*i+1)=F(2*i-1:2*i+1)+[b0;b1;b2];
end
F=F(2:end-1);
%Solve for the linear system
U=A\F;
U=[0;U;0];
figure
plot(x_m,U)
title('1D FEM qudratic interpolation')

