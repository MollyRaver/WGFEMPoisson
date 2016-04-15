function stiffnessmat_C
% This program is used to reform the mass matrix 
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global p t 
global n0 dof
global NT C
global A00_pre

Q_0=A00_pre;
% Q_0 is the mass matrix on reference element

dof_loc=n0;
% The size of local mass matrix is n0
ia=zeros(n0,n0*NT)';
% ia is the row number of each entry 
ja=zeros(n0,n0*NT)';
% ja is the column number of each entry
va=zeros(n0,n0*NT)';
% va is the value of each entry

for i=1:NT
    vtx=p(:,t(1:3,i));
    % vtx is the coordinate of three vertice
    K=0.5*abs(det([ones(1,3);vtx]));
    % K is the measure of element
    loc=n0*(i-1)+1:n0*i;
    % loc is the index of the interior basis functions
    
    [x,y]=meshgrid(loc,loc);
    ia(dof_loc*(i-1)+(1:dof_loc),:)=x;
    ja(dof_loc*(i-1)+(1:dof_loc),:)=y;
    va(dof_loc*(i-1)+(1:dof_loc),:)=Q_0*K;
    % Transform the local mass matrix into ia, ja, and va
end
C=sparse(ia,ja,va,dof,dof);
% Generate the global mass matrix