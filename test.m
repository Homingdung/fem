xa = 0;
xb = 1;
ya = 0;
yb = 1;
Nx = 40;
Ny = 40;
hx = (xb - xa)/Nx;
hy = (yb - ya)/Ny;
[x,y] = meshgrid(linspace(xa,xb,Nx + 1),linspace(ya,yb,Ny + 1));
x = reshape(x,(Nx + 1)*(Ny + 1),1);
y = reshape(y,(Nx + 1)*(Ny + 1),1);
DT = delaunayTriangulation(x,y);

draw_mesh = 0; % 绘制网格
draw_solution = 1; % 绘制数值解图像

if draw_mesh == 1
    triplot(DT)
    hold on
    vxlabels = arrayfun(@(n) {sprintf('P%d', n)}, (1:(Nx + 1)*(Ny + 1))');
    Hpl = text(x, y, vxlabels, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'BackgroundColor', 'none');
    ic = incenter(DT);
    numtri = size(DT,1);
    trilabels = arrayfun(@(x) {sprintf('T%d', x)}, (1:numtri)');
    Htl = text(ic(:,1), ic(:,2), trilabels, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'blue');
    axis([xa - 0.2*hx,xb + 0.2*hx,ya - 0.2*hy,yb + 0.2*hy])
    hold off
end

f = @(x,y) ((x - 1).*sin(x) - 2.*cos(x)).*(y - 1).*sin(y) + ((y - 1).*sin(y) - 2.*cos(y)).*(x - 1).*sin(x);
u = @(x,y) (x - 1).*(y - 1).*sin(x).*sin(y);
ux = @(x,y) (sin(x) + (x - 1).*cos(x)).*(y - 1).*sin(y);
uy = @(x,y) (sin(y) + (y - 1).*cos(y)).*(x - 1).*sin(x);

% 构造局部刚度矩阵和右端项
[N,~] = size(DT.Points); % 节点总个数
[Nk,~] = size(DT.ConnectivityList); % 单元总个数
A = zeros(N,N);
F = zeros(N,1);

for count = 1:Nk % 按单元循环
    
    i = DT.ConnectivityList(count,1);
    j = DT.ConnectivityList(count,2);
    k = DT.ConnectivityList(count,3);
    
    xi = DT.Points(i,1); yi = DT.Points(i,2);
    xj = DT.Points(j,1); yj = DT.Points(j,2);
    xk = DT.Points(k,1); yk = DT.Points(k,2);
    
    Delta = det([xi,yi,1;xj,yj,1;xk,yk,1])/2;
    
    A(i,i) = A(i,i) + ( (yj - yk)^2 + (xk - xj)^2 )/(4*Delta);
    A(i,j) = A(i,j) + ( (yj - yk)*(yk - yi) + (xk - xj)*(xi - xk) )/(4*Delta);
    A(i,k) = A(i,k) + ( (yj - yk)*(yi - yj) + (xk - xj)*(xj - xi) )/(4*Delta);
    
    A(j,i) = A(j,i) + ( (yj - yk)*(yk - yi) + (xk - xj)*(xi - xk) )/(4*Delta);
    A(j,j) = A(j,j) + ( (yk - yi)^2 + (xi - xk)^2 )/(4*Delta);
    A(j,k) = A(j,k) + ( (yk - yi)*(yi - yj) + (xi - xk)*(xj - xi) )/(4*Delta);
    
    A(k,i) = A(k,i) + ( (yj - yk)*(yi - yj) + (xk - xj)*(xj - xi) )/(4*Delta);
    A(k,j) = A(k,j) + ( (yk - yi)*(yi - yj) + (xi - xk)*(xj - xi) )/(4*Delta);
    A(k,k) = A(k,k) + ( (yi - yj)^2 + (xj - xi)^2 )/(4*Delta);
    
    F(i) = F(i) + (Delta/3)*f(xi,yi);
    F(j) = F(j) + (Delta/3)*f(xj,yj);
    F(k) = F(k) + (Delta/3)*f(xk,yk);

end

for count = 1:N
    xi = DT.Points(count,1); yi = DT.Points(count,2);
    if (xi == 0) || (xi == 1) || (yi == 0) || (yi == 1)
        A(count,:) = 0;
        A(:,count) = 0;
        A(count,count) = 1;
        F(count) = 0;
    end
end

U = A\F;
figure
trimesh(DT.ConnectivityList,x,y,U);
colormap(cool)
