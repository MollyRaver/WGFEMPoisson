function right_side_pb
% This program is used to reform the right-hand side (f,v_0)
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global gs2_pt gs2_wt
global n0 NT dou0 dof
global F rhs
global tri_area

F=zeros(dof,1);
% F is the right side

[phy_x,phy_y]=ref2phy(gs2_pt);
% Transform Gauss integration points on physical elements
idx=reshape(ones(n0,1)*(1:NT),dou0,1);
% The index of each basis functions

F_global=rhs(phy_x,phy_y);
F_global=F_global(idx,:);
% Generate the value of function rhs at each Gauss point on all physical
% elements

u0_global=u0_pb(gs2_pt(1,:)',gs2_pt(2,:)');
% Produce the value of u_0 at each Gauss point on all physical elements

F(1:dou0)=F_global.*u0_global*gs2_wt'.*tri_area(idx,:);
% Multiply the value of rhs and u_0, and then multiply the weight and area
% to get the integration value