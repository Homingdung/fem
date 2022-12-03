function meshgeneration(Nx,Ny,L)
%mesh generation
%Nx, Ny: number of nodes
%L: size of the domain
dx=L/(Nx-1);
dy=L/(Ny-1);
% mesh node
p=zeros(2,Nx*Ny);
store=zeros(3,2*(Nx-1)*(Ny-1));

k=0;
for i=1:Ny
    y=(i-1)*dy;
    for j = 1:Nx
        x=(j-1)*dx;
        k = k +1;
        p(1,k)=x;
        p(2,k)=y;
    end
end

% connectivity list
k=0;
for i=1:Ny-1
    for j=1:Nx-1
        k1=j+(i-1)*Nx;
        k2=k1+1;
        k3=k2+Nx;
        k4=k1+Nx;
        k=k+1;
        store(:,k)=[k1;k3;k4];
        k = k + 1;
        store(:,k)=[k1,k2,k3];
    end
end
figure
hold all 
patch('faces',store','vertices',p','facecolor','r','edgecolor','k');
plot(p(2,:),p(2,:),'o','color','k');

end
