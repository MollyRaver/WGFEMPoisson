function bdcond_mat_pb
% This program is used to adjust the stiffness matrix for the Dirichlet
% boundary condition
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment
global s 
global nb dou0 
global bound_idx bound_loc
global A

bound_idx=find(s(4,:)==0);
% Find all the edge on \partial\Omega


bound_loc=zeros(1,length(bound_idx)*(nb));
[X,Y]=meshgrid(nb*(bound_idx-1),1:nb);
loc_b=X+Y+dou0;
% Find the index of these edges
bound_loc(:)=loc_b;

A(bound_loc,:)=0;
A(bound_loc,bound_loc)=eye(length(bound_loc));
% Fix the stiffness for these edges. The value of edges on \partial\Omega
% are given by the right-side directly