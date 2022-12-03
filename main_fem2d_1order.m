xa=0;
xN=1;
ya=0;
yN=1;
Nx=50;
Ny=50;
hx=(xN-xa)/Nx;
hy=(yN-ya)/Ny;
[x,y] = meshgrid(linspace(xa,xN,Nx + 1),linspace(ya,yN,Ny + 1));
x = reshape(x,(Nx + 1)*(Ny + 1),1);
y = reshape(y,(Nx + 1)*(Ny + 1),1);
DT = delaunayTriangulation(x,y);
%
f=@(x,y)1;
%Assemble stiffness matrix A
%Number of nodes
[N,~]=size(DT.Points);
%Number of elements
[Ne,~]=size(DT.ConnectivityList);
A = zeros(N,N);
F = zeros(N,1);
%DT(1,:) can access the information of the first triangle
%DT(1,1) can access the information of the first node


for index=1:Ne
    
    i=DT.ConnectivityList(index,1);
    j=DT.ConnectivityList(index,2);
    k=DT.ConnectivityList(index,3);
    
    xi=DT.Points(i,1); yi=DT.Points(i,2);
    xj=DT.Points(j,1); yj=DT.Points(j,2);
    xk=DT.Points(k,1); yk=DT.Points(k,2);
    
    K = det([xi,yi,1;xj,yj,1;xk,yk,1])/2;
    
    A(i,i)=A(i,i) + ...
        ((yj-yk)^2+(xk-xj)^2)/(4*K);
    A(i,j)=A(i,j) + ...
        ((yj-yk)*(yk-yi)+(xk-xj)*(xi-xk))/(4*K);
    A(i,k)=A(i,k)+ ...
        ((yj-yk)*(yi-yj)+(xk-xj)*(xj-xi))/(4*K);
    
    A(j,i)=A(j,i)+ ...
        ((yj-yk)*(yk-yi)+(xk-xj)*(xi-xk))/(4*K);
    A(j,j)=A(j,j)+ ...
        ((yk-yi)^2+(xi-xk)^2)/(4*K);
    A(j,k)=A(j,k)+ ...
        ((yk-yi)*(yi-yj)+(xi-xk)*(xj-xi))/(4*K);
    
    A(k,i)=A(k,i)+ ...
        ((yj-yk)*(yi-yj)+(xk-xj)*(xj-xi))/(4*K);
    A(k,j)=A(k,j)+ ...
        ((yk-yi)*(yi-yj)+(xi-xk)*(xj-xi))/(4*K);
    A(k,k) = A(k,k) + ...
        ((yi-yj)^2+(xj-xi)^2)/(4*K);
    
    F(i)=F(i)+(K/3)*f(xi,yi);
    F(j)=F(j)+(K/3)*f(xj,yj);
    F(k)=F(k)+(K/3)*f(xk,yk);
end


for index = 1:N 
    xi = DT.Points(index,1); yi = DT.Points(index,2);
    if (xi == 0) || (xi == 1) || (yi == 0) || (yi == 1)
        A(index,:) = 0;
        A(:,index) = 0;
        A(index,index) = 1;
        F(index) = 0;
    end
end

U = A\F;
%plot
figure
trimesh(DT.ConnectivityList,x,y,U);


