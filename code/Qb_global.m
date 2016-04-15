function z=Qb_global(nb,n,p1,p2,fun,basis,Q_b)
% This program is used to project some continuous functions on edges
%
% Inputs:
% nb is the degree of freedom of polynomial u_b
% fun is function for projection
% n is the dimension of fun. If fun is a scalar function, then n=1
% basis is a basis function of the polynomial space
% p1 is the starting points of edges
% p2 is the endging points of edges
% Q_b is the mass matrix of basis
% 
% Outputs:
% z is the coefficients on all edges
% The column of z represents different basis functions
% The row of z represents the number of edge
% If fun is vector-valued function, then the level of z represents the
% components of fun
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global gs1_pt gs1_wt gs1_num

N=size(p1,2);
% N is the number of edges
z=zeros(nb,N,n);
for i=1:n
    z1=zeros(nb,N);
    % z1 is the projection coefficients for the i-th component of fun
    for j=1:gs1_num
        t=gs1_pt(j);
        % t represents the j-th integration point
        [px,py]=trs_1d(t,p1,p2);
        % px and py are the coordinates of j-th integration point on all
        % edges
        ug_global=repmat(basis(t)',1,N);
        % The value of projection polynomials are same for all edges, so
        % repeat N times
        fun_global=fun(px,py);
        % fun_global is the value of all edges 
        fun_global=repmat(fun_global(i,:),nb,1);
        % Take the i-th component of fun_global
        z1=z1+ug_global.*fun_global*gs1_wt(j);
    end
    z(:,:,i)=(Q_b\z1);
end


function [px,py]=trs_1d(t,p1,p2)
% This function is used to calculate the coordinate of the point with
% proportion t
% The outputs are the coordinate on all edges with starting points p1 and
% ending points p2
px=(1-t)*p1(1,:)+t*p2(1,:);
py=(1-t)*p1(2,:)+t*p2(2,:);