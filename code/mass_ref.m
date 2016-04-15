function [v0v0,vbvb,qqidx]=mass_ref
% This program is used to calculate the mass matrix on the reference
% element
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The calculation is based on symbolic operations

syms x y t real

v0=u01(x,y);
vb=ub(t);
q=nablaw_basis(x,y);

qq=int(int(q'*q,x,0,1-y),y,0,1)*2;
v0v0=eval(int(int(v0'*v0,x,0,1-y),y,0,1)*2);
vbvb=eval(int(vb'*vb,t,0,1));

for i=1:2
    for j=1:2
        qqidx(:,:,2*i+j-2)=eval(int(int(q(i,:)'*q(j,:),x,0,1-y),y,0,1)*2);
    end
end