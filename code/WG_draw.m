function WG_draw(uh)
% This function is used to plot the figure of uh
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global p s t NT NE
global dou0 nb

figure
hold on
for i=1:NT
    vtx=p(:,t(1:3,i));
    loc_0=loc_u0(i);
    trimesh([1 2 3],p(1,t(1:3,i)),p(2,t(1:3,i)),...
        u01([0;1;0],[0;0;1])*uh(loc_0));
      
end
% Draw u_0
for i=1:NE
    p1=p(:,s(1,i));
    p2=p(:,s(2,i));
    loc_b=(1:nb)+nb*(i-1)+dou0;
    plot3([p1(1) p2(1)],[p1(2) p2(2)],ub([0;1])*uh(loc_b));
    
    
end
% axis([0 1 0 1 0 1])
hold off
view(-38,30)
% Draw u_b

figure
hold on
hold on
for i=1:NT
    vtx=p(:,t(1:3,i));
    loc_0=loc_u0(i);

    h = trisurf([1 2 3],p(1,t(1:3,i)),p(2,t(1:3,i)),...
        u01([0;1;0],[0;0;1])*uh(loc_0),...
                'FaceColor', 'interp', 'EdgeColor', 'interp');
      
end
% Draw u_0 as surfaces