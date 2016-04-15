function z=WGElliptic(n,d0,db,dnw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is designed to solve the Laplaican problem
% -\Delta u=f in \Omega
%         u=g on \partial\Omega.
%
%----The input arguments:
% n is the mesh grid parameter
%   The domain is [0,1]*[0,1] with uniform triangle partition,
%   and the length of right angle side is set to be 1/n.
% d0 is the degree of u_0, i.e. u_0\in P_{d0}(T)
%   d0 should be an integer, and 1<=d0<=5.
% db is the degree of u_b, i.e. u_b\in P_{db}(e)
%   db should be an integer, and 1<=db<=5.
% dnw is the degree of weak gradient, i.e. \nabla_w v\in [P_{dnw}(T)]^2
%   dnw should be an integer, and 0<=dnw<=5.
% There are several choices for the WG finite element (u0,ub,grad_w u)
% eg. PkPkPk element
%     PkPkPk-1 element
%     PkPk-1Pk-1 element
%     and others
%
%----The output arguments:
% z(1) is the error corresponding to the energy norm
% z(2) is the error corresponding to the L2 norm
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear global

global err dof
global d1 d2 d3


d1=d0;
d2=db;
d3=dnw;
% load the inputs

pre_setting;
% Some parameter settings before calculation

tic
% Start timing

WG_partition(n);
% Mesh partition

stiffness_matrix;
% Reform the stiffness matrix

stiffnessmat_C;
% Reform the mass matrix

right_side_pb;
% Reform the right-hand side

bdcond_mat_pb;
% Adjust the stiffness matrix according to the boundary condition

bdcond_rhs_pb;
% Adjust the right side according to the boundary condition

solve_WG;
% Solving the matrix

toc
% End timing

error_est;
% Calculate the errors in triple-bar norm and L^2 norm


z=err;
% Output

clear global