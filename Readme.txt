
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 This program is designed to solve the Laplaican problem by the weak Galerkin method
 -\Delta u=f in \Omega
         u=g on \partial\Omega.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

User Guide:

The main program is WGElliptic.m. The syntax is 

  err=WGElliptic(n,d0,db,dnw);

The input arguments:

 n is the mesh grid parameter
   The domain is [0,1]*[0,1] with uniform triangle partition,
   and the length of right angle side is set to be 1/n.
 d0 is the degree of u_0, i.e. u_0\in P_{d0}(T)
   d0 should be an integer, and 1<=d0<=5.
 db is the degree of u_b, i.e. u_b\in P_{db}(e)
   db should be an integer, and 1<=db<=5.
 dnw is the degree of weak gradient, i.e. \nabla_w v\in [P_{dnw}(T)]^2
   dnw should be an integer, and 0<=dnw<=5.

The output arguments:

 err is a 1*2 vector containing the triple-bar norm (which is equivalent to H^1 norm) and the L^2 norm of the error.

Moreover, if a figure for numerical solution is necessary, please set the variable draw=1 in presetting.m. The Dirichlet boundary condition and the source term can also be alternated in presetting.m. The detailed description can be found in presetting.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Algorithm introduction:

In this program, we consider a uniform triangle mesh for the unit square domain [0,1]*[0,1]. The descripitions about partition can be found in partition.m.

Denote T_h the triangle partition. T is the element in T_h, and e is the edge in T_h. For each element T, denote h_T its diameter.

Now, we introduce the weak Galerkin method for Laplacian problem briefly. First, define the following weak functions space 

V_h=\{ v=(v_0,v_b);v_0\in P_k(T), v_b\in P_j(e)\},

V_h=\{ v=(v_0,v_b);v_0\in P_k(T), v_b\in P_j(e), v_b=0 on \partial\Omega\},

where j,k are given non-negative integers.

Define a vector valued function space G=[P_l(T)]^2, where l is an integer. Define the weak gradient operator \nabla_w: V_h->G by the following equation.
Suppose v\in V_h, for any \phi\in G, define \nabla_w v\in G, such that on each element T\in T_h

(\nabla_w v,\phi)_T = -(v_0,\nabla\cdot\phi)_T+( v_b,\phi\cdot n)_\partial T,

where (\cdot,\cdot) denotes the L^2 inner-product.

We also need to define the L^2 projection operator Q_m on each edge e. Denote m=max{j,l}, and on each edge e, Q_m is L^2 projection from L^2(e) onto P_m(e).

Now, define the following two bilinear forms:

For any v,w\in V_h,

    s(v,w)=\sum_{T\in T_h} h_T^{-1} (Q_m v_0-v_b, Q_m w_0-w_b)_\partial T,

    a(v,w)=\sum_{T\in T_h} (\nabla_w v,\nabla_w w)+s(v,w).

The weak Galerkin algorithm for Laplacian problem is:

Find u_h\in V_h, such that u_h=Q_j g on \partial\Omega, and for any v\in V_h, 

    a(u_h,v)=(f,v_0).

If there is any questions regarding the code, please contact:
Qilong Zhai, 
Ruishu Wang,
Lin Mu, mul1@ornl.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

How to cite this code:

Zhai, Q, Wang, R, and Mu, L. (2016) Implementation of Weak Galerkin Finite Element Methods [Computer program]. Available at https://sites.google.com/site/weakgalerkin/home/codes/WGFEMPoisson.zip.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


The theoretical analysis refers to 
Mu, L., Wang, J., & Ye, X. (2015). A weak Galerkin finite element method with polynomial reduction. Journal of Computational and Applied Mathematics, 285, 45¨C58. http://doi.org/10.1016/j.cam.2015.02.001.

Mu, L., Wang, J., & Ye, X. (2015). Weak Galerkin finite element methods on polytopal meshes. International Journal of Numerical Analysis and Modeling, 12(1), 31¨C53. Retrieved from http://arxiv.org/abs/1204.3655