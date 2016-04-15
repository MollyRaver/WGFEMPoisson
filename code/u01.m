function z=u01(x,y)
% This function is the basis of u_0
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global d1

switch d1
    case 1
        z=[ones(size(x)) x y];
    case 2
        z=[ones(size(x)) x y x.^2 x.*y y.^2];
    case 3
        z=[ones(size(x)) x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3];
    case 4
        z=[ones(size(x)) x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4];
    case 5
        z=[ones(size(x)) x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4 ...
           x.^5 x.^4.*y x.^3.*y.^2 x.^2.*y.^3 x.*y.^4 y.^5];
end