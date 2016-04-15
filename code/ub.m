function z=ub(t)
% This function is the basis of u_b
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global d2

switch d2
    case 0
        z=ones(size(t));
    case 1
        z=[t 1-t];
    case 2
        z=[t t.*(1-t) 1-t ];
    case 3
        z=[t.^3 t.^2.*(1-t) t.*(1-t).^2 (1-t).^3 ];
    case 4
        z=[t.^4 t.^3.*(1-t) t.^2.*(1-t).^2 t.*(1-t).^3 (1-t).^4 ];
    case 5
        z=[t.^5 t.^4.*(1-t) t.^3.*(1-t).^2 t.^2.*(1-t).^3 t.*(1-t).^4 (1-t).^5];
end