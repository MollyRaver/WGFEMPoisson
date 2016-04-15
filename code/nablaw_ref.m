function z=nablaw_ref
% This program is used to calculate the weak gradient on the reference
% element
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The calculation is based on symbolic operations

syms x y t real

v0=u01(x,y);
vb=ub(t);
q=nablaw_basis(x,y);

n1=[sqrt(2)/2 sqrt(2)/2];
n2=[-1 0];
n3=[0 -1];
l1=sqrt(2);
l2=1;
l3=1;

divq=diff(q(1,:),x)+diff(q(2,:),y);
qn1=subs(n1*q,{x,y},{1-t,t});
qn2=subs(n2*q,{x,y},{0,1-t});
qn3=subs(n3*q,{x,y},{t,0});

qq=int(int(q'*q,x,0,1-y),y,0,1);
divqv0=-int(int(divq'*v0,x,0,1-y),y,0,1);
qn1vb=l1*int(qn1'*vb,t,0,1);
qn2vb=l2*int(qn2'*vb,t,0,1);
qn3vb=l3*int(qn3'*vb,t,0,1);

z=eval(qq\[divqv0 qn1vb qn2vb qn3vb]);
