function Qhu=Qh1(fun)
% This function is used to calculate the Q_h projection 
% 
% Input:
% fun is the continuous function for projection
% 
% Output:
% Qhu is the vector of coefficients of u_0 and u_b 
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dof dou0 doub  nb
global p s  NT NE
global gs2_pt gs2_wt
global n0 
global Q_0 Q_b

Qhu=zeros(dof,1);

[phy_x,phy_y]=ref2phy(gs2_pt);
% Transform Gauss integration points on physical elements
idx=reshape(ones(n0,1)*(1:NT),dou0,1);
% The index of each basis functions

u_global=fun(phy_x,phy_y);
u_global=u_global(idx,:);
% Generate the value of function fun at each Gauss point on all physical
% elements

u0_global=u0_pb(gs2_pt(1,:)',gs2_pt(2,:)');
% Produce the value of u_0 at each Gauss point on all physical elements

z=reshape(u_global.*u0_global*gs2_wt',n0,NT);
% Multiply the value of rhs and u_0, and then multiply the weight and area
% to get the integration value
Qhu(1:dou0)=Q_0\z;
% Divide the mass matrix Q_0 to get the coefficients of u_0


p1=p(:,s(1,:));
% p1 is the coordinate of the starting points of all edges 
p2=p(:,s(2,:));
% p2 is the coordinate of the ending points of all edges 


Qhu_bound=zeros(doub,1);

z=Qb_global(nb,1,p1,p2,fun,@ub,Q_b);
% Calculate the L^2 projection of exact solution on all edges
Qhu_bound(1:nb*NE)=z;
% The corresponding position is assigned by the projection values

Qhu(dou0+1:end)=Qhu_bound;