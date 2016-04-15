function [quad_xy_2d quad_w_2d quad_num_2d]=integral_transform(rule)
% Set the 2d Gauss integral points on triangle elements
% There are two types of rule 
% Rule=1 represents the 27 points Gauss integration
% Rule=2 represents the 16 points Gauss integration
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gs2_wt gs2_pt gs2_n

[point_type sub_xyz sub_w]=integral_rule(rule);
quad_num_2d=point_type*[1;3;6];
quad_xy_2d=zeros(2,quad_num_2d);
quad_w_2d=zeros(1,quad_num_2d);
pt_order=[1  2  3  2  3  1;
          2  3  1  1  2  3];
m=0;
n=0;
for i=1:point_type(1)
    location=[1/3;1/3;1/3];
    quad_xy_2d(:,i+m)=location(1:2);
    quad_w_2d(i+m)=sub_w(i+n);
end
m=point_type(1);
n=point_type(1);
for i=1:point_type(2)
    a=sub_xyz(1,i+n);
    location=[a;a;1-2*a];
    for j=1:3
        quad_xy_2d(:,3*(i-1)+j+m)=location(pt_order(:,j));
        quad_w_2d(3*(i-1)+j+m)=sub_w(i+n);
    end
end
m=point_type(1)+3*point_type(2);
n=point_type(1)+point_type(2);
for i=1:point_type(3)
    a=sub_xyz(1,i+n);
    b=sub_xyz(2,i+n);
    location=[a;b;1-a-b];
    for j=1:6
        quad_xy_2d(:,6*(i-1)+j+m)=location(pt_order(:,j));
        quad_w_2d(6*(i-1)+j+m)=sub_w(i+n);
    end
end

gs2_pt=quad_xy_2d;
gs2_wt=quad_w_2d;
gs2_n=quad_num_2d;