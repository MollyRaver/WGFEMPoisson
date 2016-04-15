function stiffness_matrix
% This program is used to reform the stiffnes matrix
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global p t
global NT 
global n0 nb  
global dof
global S_pre A00_pre nabla_w_pre DwDw_pre 
global A
global h
global alpha gamma para
global tri_area

h=0;
% h is the mesh size

dof_loc=n0+3*nb;
% dof_loc is the size of local stiffness matrix
ia=zeros(dof_loc,dof_loc*NT)';
% ia is the row number of each entry 
ja=zeros(dof_loc,dof_loc*NT)';
% ja is the column number of each entry
va=zeros(dof_loc,dof_loc*NT)';
% va is the value of each entry
tri_area=zeros(NT,1);
% tri_area is the measure of all the elements

for i=1:NT
% Loop by elements, and generate the local stiffness matrix on each
% element
    P=p(1:2,t(1:3,i));
    % The coordinate of three vertice 
    % P=[x1 x2 x3;...
    %    y1 y2 y3];
    E=P(:,[3 1 2])-P(:,[2 3 1]);
    % The vector form of three edges
    L=sqrt(sum(E.^2));
    % The length of three edges
    K=0.5*abs(det([ones(1,3);P]));
    tri_area(i)=K;
    % The measure of element
    h_T=0.5*prod(L)/K;
    % The diameter of element
    h=max([h_T,h]);
    
    D=inv([E(:,3) -E(:,2)])';
    % D is the affine transformation from physical element to reference
    % element
    nabla_w=nabla_w_pre';
    % The coefficient of weak gradient
    T=(D'*D);
    DwDw=(T(1,1)*DwDw_pre(:,:,1)+T(1,2)*DwDw_pre(:,:,2)+T(2,1)*DwDw_pre(:,:,3)+T(2,2)*DwDw_pre(:,:,4))*K;
    % The inner-product of weak gradient basis functions on physical element
   
    A00=A00_pre*K*gamma;
    % The inner-product of u_0 on physical element
    Abb=nabla_w*DwDw*nabla_w'*alpha;
    % The inner-product of weak gradient on physical element
    
    S_loc=zeros(n0+3*nb);
    % S_loc is the local stabilizer
    for j=1:3
        S_loc=S_loc+S_pre(:,:,j)*L(j)/h_T;
    end
    % S_loc is transformed from the stabilizer on reference element
    
    A_loc=blkdiag(A00,zeros(3*nb))+Abb;
    A_loc=A_loc+S_loc*para;
    % Combining A_loc and S_loc together as the local stiffness matrix
    
    loc_0=loc_u0(i);
    % loc_0 is the index of interior basis functions
    loc_b=loc_ub(i);
    % loc_b is the index of edge basis functions
    loc=[loc_0 loc_b];
    % loc is the index of all basis functions
    
    for j=1:3
        if t(6+j,i)==-1
            loc([n0+nb*j-nb+1:n0+nb*j])=loc([n0+nb*j:-1:n0+nb*j-nb+1]);
        end
    end
    % Adjust the order of basis functions to make it anti-clock wise
    
    [x,y]=meshgrid(loc,loc);
    ia(dof_loc*(i-1)+(1:dof_loc),:)=x;
    ja(dof_loc*(i-1)+(1:dof_loc),:)=y;
    va(dof_loc*(i-1)+(1:dof_loc),:)=A_loc;
    % Transform the data in A_loc into ia, ja, and va

end
A=sparse(ia,ja,va,dof,dof);
% Generate the global stiffness matrix
