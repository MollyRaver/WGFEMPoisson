function [s t]=sort_edge(e,t)
%      by Qilong Zhai, Ruishu Wang, and Lin Mu
%                 04/14/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=size(e,2);
N=size(t,2);
alledges=zeros(3,3*N);
s=zeros(4,(3*N+M)/2);
e1=zeros(2,M+1);
for m=1:N
    alledges(:,3*m-2)=[t([1 2],m);m];
    alledges(:,3*m-1)=[t([1 3],m);m];
    alledges(:,3*m)=[t([2 3],m);m];
end
e1(:,1:M)=[min(e(1:2,:));max(e(1:2,:))];
alledges(1:2,:)=[min(alledges(1:2,:));max(alledges(1:2,:))];
e1(:,1:M)=sortrows(e1(:,1:M)')';
alledges=sortrows(alledges')';
j=1;
k=1;

t(4:6,:)=zeros(3,N);
ind=4*ones(1,N);
for i=1:(3*N+M)/2
    if any(alledges(1:2,j)~=e1(:,k))
        s(:,i)=[alledges(:,j);alledges(3,j+1)];
        t(ind(alledges(3,j)),alledges(3,j))=i;
        t(ind(alledges(3,j+1)),alledges(3,j+1))=i;
        ind([alledges(3,j) alledges(3,j+1)])=ind([alledges(3,j) alledges(3,j+1)])+1;
        j=j+2;
    else
        s(:,i)=[alledges(:,j);0];
        t(ind(alledges(3,j)),alledges(3,j))=i;
        ind(alledges(3,j))=ind(alledges(3,j))+1;
        j=j+1;
        k=k+1;
    end
end