function loc_b=loc_ub(i)
% This program is used to calculate the index of edge basis functions
% u_b
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% i is the number of element 
% Output:
% loc_b is the vector of index of edge basis functions
global nb dou0 t

loc_b=reshape(nb*ones(nb,1)*t(4:6,i)'+(1:nb)'*ones(1,3)-nb,1,3*nb)+dou0;