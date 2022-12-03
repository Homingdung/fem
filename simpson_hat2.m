
function y = simpson_hat2(x1,x2)

% This function evaluate \int_{x1}^{x2} f(x) \phi(x) dx,
% where   \phi(x) = (x-x1)/(x2-x1), using the Simpson rule

  xm = (x1+x2)*0.5;
  y = (x2-x1)*(f(x1)*hat2_1order(x1,x1,x2) + 4*f(xm)*hat2_1order(xm,x1,x2)...
               + f(x2)*hat2_1order(x2,x1,x2) )/6;

  return

  




