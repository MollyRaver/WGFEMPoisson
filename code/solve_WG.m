function solve_WG
% This function is used to solve the numerical solution from the stiffness
% matrix A and right side F
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global A F 
global uh

uh=A\F;
% Use Gauss elimination method to solve the problem directly. If you want
% to use iteration methods, just fix here

