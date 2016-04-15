function z=u0_pb(x,y)
% This program is used to produce the value of u_0 on all elements
%
% Inputs:
% x,y are coordinates on reference element
%
% Outputs:
% z is the vector of u_0 at (x,y) for each element
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NT

z=u01(x,y)';
% The value of u_0 at (x,y) on reference element
z=repmat(z,NT,1);
% Repeat for NT times