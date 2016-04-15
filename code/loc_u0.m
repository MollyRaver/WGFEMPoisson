function loc_0=loc_u0(i)
% This program is used to calculate the index of interior basis functions
% u_0
%
% Input:
% i is the number of element 
% Output:
% loc_0 is the vector of index of interior basis functions
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global n0 

loc_0=(1:n0)+n0*(i-1);