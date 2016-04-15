function bdcond_rhs_pb
% This program is used to adjust the right-side for the Dirichlet boundary
% condition
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global p s
global F
global nb 
global bound_idx bound_loc
global solu 
global Q_b

p1=p(:,s(1,bound_idx));
% p1 is the coordinate of the starting points of all edges on
% \partial\Omega
p2=p(:,s(2,bound_idx));
% p2 is the coordinate of the ending points of all edges on \partial\Omega
N=length(bound_idx);
% N is the number of edges on \partial\Omega

F_bound=zeros(size(bound_loc));

fun=@(x,y)solu(x,y);
z=Qb_global(nb,1,p1,p2,fun,@ub,Q_b);
% Calculate the L^2 projection of exact solution on boundary edges

F_bound(1:nb*N)=z;
F(bound_loc)=F_bound;
% The corresponding position is assigned by the projection values