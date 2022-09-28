function [v]=V_residue(t)
% t=7199;
v=zeros(1,6);
m=zeros(t+1,6);
x=zeros();
x(1)=0;
load('Vi.mat')  % 供油曲线数据
load('INIT.mat') % 飞行器参数
m(1,:)=INIT(:,7)'*850;
A=[1 -1 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 -1 1];
for t1=1:t
    x(t1+1)=t1;
    m(t1+1,:)=m(t1,:)-Vi(t1,:)*A;

end
% v(1,:)=m(t+1,:)/850;
 v(1,:)=m(t+1,:);
