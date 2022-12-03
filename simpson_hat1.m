
function y = simpson_hat1(x1,x2)

  xm = (x1+x2)*0.5;
  y = (x2-x1)*(f(x1)*hat1_1order(x1,x1,x2) + 4*f(xm)*hat1_1order(xm,x1,x2)...
               + f(x2)*hat1_1order(x2,x1,x2) )/6;

  return

  




