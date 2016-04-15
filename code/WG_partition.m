function WG_partition(n)
% This function is for the uniform partition
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global p s t
global NT NE
global n0 nb
global dou0 doub dof

[p s t]=uniform_mesh1(n,0,1);
% Partition the unit square [0,1]*[0,1] into 2*n^2 triangles
%
% p is the vertice of partiton, the first row represents x-axis coordinate
% and the second row represents y-axis coordinate
%
% s is the edges of partition. The first two rows are the numbers of
% endpoins, and the next two rows are numbers of elements by this edge
%
% t is the elements of partition. The first three rows are numbers of
% vertice counterclockwise, and the next three rows are numbers of edges

NT=size(t,2);      
% NT is the number of triangle elements
NE=size(s,2);            
% NE is the number of edges

dou0=n0*NT;
% Degree of freedoms of interior basis functions
doub=nb*NE;
% Degree of freedoms of edge basis function
dof=dou0+doub;
% Total degree of freedoms