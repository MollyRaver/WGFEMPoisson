function z=stable_ref1
% This program is used to calculate the stabilizer on the reference
% element, and Q_b is used in the stabilizer
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The calculation is based on symbolic operations

syms x y t real

v0=u01(x,y);
vb=ub(t);

v01=subs(v0,{x,y},{1-t,t});
v02=subs(v0,{x,y},{0,1-t});
v03=subs(v0,{x,y},{t,0});

vbvb=int(vb'*vb,t,0,1);

Qbv01=vbvb\int(vb'*v01,t,0,1);
Qbv02=vbvb\int(vb'*v02,t,0,1);
Qbv03=vbvb\int(vb'*v03,t,0,1);

zr=zeros(length(vb),length(vb));
v1=[Qbv01 -eye(length(vb)) zr zr];
v2=[Qbv02 zr -eye(length(vb)) zr];
v3=[Qbv03 zr zr -eye(length(vb))];

z(:,:,1)=v1'*vbvb*v1;
z(:,:,2)=v2'*vbvb*v2;
z(:,:,3)=v3'*vbvb*v3;

z=double(z);
