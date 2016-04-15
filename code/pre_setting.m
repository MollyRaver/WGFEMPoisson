function pre_setting
% This program is used to set some parameters, including user parameters
% and some data for the calculation
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global n0 nb
global S_pre A00_pre
global nabla_w_pre DwDw_pre 
global draw 
global solu rhs
global alpha beta gamma para
global gs1_pt gs1_wt gs1_num
global Q_0 Q_b
global d1 d2 d3
% All the global variables are set in this program

format long e
% All the outputs are double format


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters below are user parameters, please set them properly for
% your problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

draw=1;
% This variable determines if there will be a picture for the numerical
% solution

solu=@u_real;
% This variable is the function for the exact solution. It shall 
% affect the boundary condition and the error estimate
rhs=@f;
% This variable is the function for the right side. 

alpha=1;
% This variable is the coefficent of \Delta u in the equation
gamma=0;
% This variable is the coefficent of u in the equation
% alpha and gamma should be constant

para=1;
% This variable is the coefficent of the stabilizer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters below are program parameters. Do not change these parameters
% unless you changed the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n0=length(u01(0,0));
% n0 is the degree of freedom on each T, and is affected by d0
nb=length(ub(0));
% nb is the degree of freedom on each e, and is affected by db

[A00_pre,Q_b,DwDw_pre]=mass_ref;
Q_0=A00_pre;
% Calculate the mass matrix on reference element

nabla_w_pre=nablaw_ref;
% Calculate the weak gradient on reference element
db=d2;
dnw=d3;



if db>dnw
    S_pre=stable_ref1;
else
    S_pre=stable_ref2;
end
% Calculate the stabilizer on reference element


gs1_wt=[0.118463442528095,0.239314335249683,...
    0.284444444444445,0.239314335249683,...
    0.118463442528095];
gs1_pt=[0.046910077030668 0.230765344947158 0.5 ...
    0.769234655052842 0.953089922969332];
gs1_num=5;
% Set the 1d Gauss integration points and their weights
integral_transform(1);
% Set the 2d Gauss integration points and their weights
