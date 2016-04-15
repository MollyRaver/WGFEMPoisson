function [p1 s1 t1]=uniform_mesh1(n,left,right,draw)
% This function is used to produce a uniform triangle partition of a square
% domain
%
% Inputs:
% n is the number of partition on square edge
% left and right are the coordinates of square
% draw is the flag for plotting. If draw=1, a figure for partition shall be
% plotted
%
% Outputs:
% p1 is the information for vertice
% s1 is the information for edges
% t1 is the information for elements
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write=0;
% If write=1, the corresponding .ele, .edge, and .node files shall be
% producted. They are used for C++ programs.
if nargin<4
    draw=0;
end
if nargin<3
    right=1;
end
if nargin<2
    left=0;
end
if nargin<1
    n=10;
end
h=(right-left)/n;
p=zeros(3,(n+1)^2);
t=zeros(3,2*n^2);
for i=1:n+1
    for j=1:n+1
        m=tr(i,j,n);
        p(1,m)=right-(j-1)*h;
        p(2,m)=left+(i-1)*h;
        p(3,m)=0;
        if (i-1)*(j-1)*(i-n-1)*(j-n-1)==0
            p(3,m)=1;
        end
    end
end
for i=1:n
    for j=1:n
        m=trs(i,j,n);
        k=2*n^2-trs1(n+1-i,n+1-j,n)+1;
        t(:,m)=[tr(i,j,n);tr(i+1,j,n);tr(i,j+1,n)];
        t(:,k)=[tr(i,j+1,n);tr(i+1,j,n);tr(i+1,j+1,n)];
    end
end
e=zeros(2,4*n);
for i=1:n
    e(:,i)=[tr(i,1,n);tr(i+1,1,n)];
    e(:,n+i)=[tr(n+1,i,n);tr(n+1,i+1,n)];
    e(:,2*n+i)=[tr(n+2-i,n+1,n);tr(n+1-i,n+1,n)];
    e(:,3*n+i)=[tr(1,n+2-i,n);tr(1,n+1-i,n)];
end
[s t]=sort_edge(e,t);
for i=1:size(t,2)
    [~,order]=sort(-t(1:3,i));
    order_inv(order)=1:3;
    t(4:6,i)=t(3+order_inv,i);
end
ptorder=[1 2 3 1 2 3];
for i=1:size(t,2)
    for j=1:3
        if t(ptorder(j+2),i)<t(ptorder(j+1),i)
            t(j+6,i)=-1;
        else
            t(j+6,i)=1;
        end
    end
end


if draw==1;
    save
    hold on
    for i=1:size(p,2)
        temp=num2str(i);
        text(p(1,i)+0.1*h,p(2,i)-0.1*h,temp,'Color','blue');
        plot(p(1,i),p(2,i),'r*');
    end
    for i=1:size(s,2)
        temp=num2str(i);
        text(sum(p(1,s(1:2,i)))/2,sum(p(2,s(1:2,i)))/2,temp,'Color','green');
    end
    for i=1:size(t,2)
        temp=num2str(i);
        text(sum(p(1,t(1:3,i)))/3,sum(p(2,t(1:3,i)))/3,temp,'Color','red');
        plot(p(1,t([1 2],i)),p(2,t([1 2],i)));
        plot(p(1,t([2 3],i)),p(2,t([2 3],i)));
        plot(p(1,t([1 3],i)),p(2,t([1 3],i)));
    end
    hold off
end

p1=p(1:2,:);
s1=s;
t1=t;
p=[size(p,2) 0 0 0;(1:size(p,2))' p'];
t=[size(t,2) 0 0 0;(1:size(t,2))' t(1:3,:)'];
s=[size(s,2) 0 0 0;(1:size(s,2))' s(1:2,:)' ~s(4,:)'];

if write==1
    fid=fopen('fem.1.node','w');
    for i=1:p(1,1)+1
        fprintf(fid,'%d ',p(i,1));
        for j=2:3
            fprintf(fid,'%.15f ',p(i,j));
        end
        fprintf(fid,'%d ',p(i,4));
        fprintf(fid,'\r\n');
    end
    fclose(fid);
    fid=fopen('fem.1.ele','w');
    for j=1:3
        fprintf(fid,'%d ',t(1,j));
    end
    fprintf(fid,'\r\n');
    for i=2:t(1,1)+1
        for j=1:4
            fprintf(fid,'%d ',t(i,j));
        end
        fprintf(fid,'\r\n');
    end
    fclose(fid);
    fid=fopen('fem.1.edge','w');
    for i=1:s(1,1)+1
        for j=1:4
            fprintf(fid,'%d ',s(i,j));
        end
        fprintf(fid,'\r\n');
    end
    fclose(fid);
end

function m=tr(i,j,n)
if i+j<n+3
    m=(i+j-2)*(i+j-1)/2+j;
    if mod(i+j,2)==1
        m=(2*i+2*j-2)*(i+j-1)/2+1-m;
    end
else
    m=(n+1)^2-tr(n+2-i,n+2-j,n)+1;
end

function m=trs(i,j,n)
if i+j<n+2
    m=(i+j-2)^2+2*j-1;
    if mod(i+j,2)==1
        m=(i+j-1)^2+(i+j-2)^2+1-m;
    end
else
    m=2*n^2-((2*n+1-i-j)^2+2*(n+1-j))+1;
    if mod(i+j,2)==1
        m=4*n^2-((2*n+2-i-j)^2+(2*n+1-i-j)^2-1)-m;
    end
end

function m=trs1(i,j,n)
if i+j<n+2
    m=(i+j-2)^2+2*j-1;
    if mod(i+j,2)==0
        m=(i+j-1)^2+(i+j-2)^2+1-m;
    end
else
    m=2*n^2-((2*n+1-i-j)^2+2*(n+1-j))+1;
    if mod(i+j,2)==0
        m=4*n^2-((2*n+2-i-j)^2+(2*n+1-i-j)^2-1)-m;
    end
end