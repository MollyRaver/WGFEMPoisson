function z=stable_ref2
% This program is used to calculate the stabilizer on the reference
% element, and Q_w (Projection) is used in the stabilizer
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The calculation is based on symbolic operations

syms x y t real

v0=u01(x,y);
vb=ub(t);
vl=ul(t);

v01=subs(v0,{x,y},{1-t,t});
v02=subs(v0,{x,y},{0,1-t});
v03=subs(v0,{x,y},{t,0});

vlvl=int(vl'*vl,t,0,1);

Qbv01=vlvl\int(vl'*v01,t,0,1);
Qbv02=vlvl\int(vl'*v02,t,0,1);
Qbv03=vlvl\int(vl'*v03,t,0,1);

zr=zeros(length(vl),length(vb));
vbvl=vlvl\int(vl'*vb,t,0,1);

v1=[Qbv01 -vbvl zr zr];
v2=[Qbv02 zr -vbvl zr];
v3=[Qbv03 zr zr -vbvl];

z(:,:,1)=v1'*vlvl*v1;
z(:,:,2)=v2'*vlvl*v2;
z(:,:,3)=v3'*vlvl*v3;

z=double(z);

