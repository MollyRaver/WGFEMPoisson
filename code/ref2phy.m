function [phy_x,phy_y]=ref2phy(ref_xy)
% This program is used to transform reference coordinates into physical
% corrdinates
%
% Inputs:
% ref_xy is the reference coordinates, x-axis for the first row and y-axis
% for the second row
%
% Outputs:
% phy_x is x-axis coordinates on all physical elements
% phy_y is y-axis coordinates on all physical elements
%
% if ref_xy is a 2*m matrix, and there are NT triangle elements, then phy_x
% and phy_y are NT*m matrice.
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global p t

px=p(1,:);
py=p(2,:);
% px, py are the coordiantes of all vertice
area_xyz=[1-sum(ref_xy);ref_xy];
% Transform the reference coordinates into area coordiantes
phy_x=px(t(1:3,:))'*area_xyz;
phy_y=py(t(1:3,:))'*area_xyz;
% The physical coordiantes are the linear combination of area coordinates