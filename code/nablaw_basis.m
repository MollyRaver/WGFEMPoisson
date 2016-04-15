function z=nablaw_basis(x,y)
% This function is the basis of weak gradient space
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global d3

switch d3
    case 0
        z=[1 0;
            0 1];
    case 1
        z=[1 x y 0 0 0;
            0 0 0 1 x y];
    case 2
        z=[1 x y x^2 x*y y^2 0 0 0 0   0   0;
            0 0 0 0   0   0   1 x y x^2 x*y y^2];
    case 3
        z=[1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 0 0 0 0   0   0   0   0     0     0   ;
            0 0 0 0   0   0   0   0     0     0   1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3];
    case 4
        z=[1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4];
    case 5
        z=[1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4 x.^5 x.^4.*y x.^3.*y.^2 x.^2.*y.^3 x.*y.^4 y.^5 ...
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  ...
          1 x y x.^2 x.*y y.^2 x.^3 x.^2.*y x.*y.^2 y.^3 x.^4 x.^3.*y x.^2.*y.^2 x.*y.^3 y.^4 x.^5 x.^4.*y x.^3.*y.^2 x.^2.*y.^3 x.*y.^4 y.^5];
end