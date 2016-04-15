function error_est
% This function is used to draw the solution and calculate the errors
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dou0 doub 
global uh A C
global solu 
global draw
global err


Qhu=Qh1(solu);
% Qhu is the Q_h projection of exact solution solu
eh=Qhu-uh;
% eh is the error between Q_h u and u_h

if draw~=0
    close all
    % Close all figures
    WG_draw(uh);
    % Draw the numerical solution
end

% Calculate the triple-bar error |||e_h|||^2=a(e_h,e_h)
% The triple-bar norm can be calculated by the quadratic form of A that
% |||e_h|||^2=eh'*A*eh;
loc_u=1:dou0+doub;
err_trb=eh(loc_u)'*A(loc_u,loc_u)*eh(loc_u);
err_trb=sqrt(err_trb);
% err_trb is the triple-bar norm of e_h

% Similarly, The L^2 norm can be calculated by the quadratic form of C that
% ||e_0||^2=eh'*C*eh;
err_L2=sqrt(eh'*C*eh);
% err_L2 is the L^2 norm of e_h

err=[err_trb err_L2];