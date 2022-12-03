# fem
Finite element method for 1D and 2D problems

# Problem Statement 
$$
-\Delta u=f, x\in \Omega=(0,1), u(0)=u(1)=0
$$
1D problem
```main_fem1d_1order.m```: Finite element method for 1D problem  with linear interpolation.

```main_fem1d_2order.m``` : Finite element method for 1D problem with qudratic interpolation. 

2D problem
```meshgeneration.m```: Generate the mesh. 
However, we have more advanced tool: working with Delaunay triangulations in MATLAB.
[Delaunay Triangulations](https://uk.mathworks.com/help/matlab/math/delaunay-triangulation.html)

```main_fem2d_1order.m``` : Finite element method for 2D problem with **Delaunay Triangulations**. 

# References
[1] Li, Zhilin, Zhonghua Qiao, and Tao Tang. _Numerical solution of differential equations: introduction to finite difference and finite element methods_. Cambridge University Press, 2017.
