x=linspace(0,1,20);
U=fem1d(x);
sol=[];
for i=1:length(x)
    sol(i)=fem_soln(x,U,x(i));
end
figure
plot(x,sol);
title('1D FEM linear interpolation')