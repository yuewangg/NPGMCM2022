clear all;clc;close all ;
load("../data/data_q1.mat");
%%阵列天线参数
M=86;% 阵列天线数量
K=256;%快拍数
%theta=[-3 1];%信号源与天线轴方向夹角
%dis=[60 61];%信号源
N=1;%信号源数目
%f0=24e9;%射频
c=2.99792458e8;
%B=150e6;%带宽
%lambda=c/f0;%波长
%FT=2e-3;%一个锯齿波持续时间
%Fs=3.2e4;%采样频率
d=0.0815;%阵元间距
%%接收信号
X=zeros(M,K);
A=zeros(M,N);
s=zeros(N,K);

%X=A*s
for i=1:M
    for j=1:K
        X(i,j)=Z(i,j);
    end
end
snr=20; %信噪比
%X1=awgn(X,snr,'measured');
X1=X;
%%构造伪谱
mm=((s*s')/K);
R=X1*X1'/K;
[V,D]=eig(R);%eig 求特征值 D，特征向量 V
Eig=diag(D);%特征值
Eig=sort(Eig,'descend'); %从大到小排特征值
V=fliplr(V);%特征向量对应排
%%信源数估计 SORTE 方法估计出相干信源组一个独立信源看作是一个相干信源组 in[2]
D=diag(D);
D=sort(D,'descend');
lambda_1=zeros(1,M-1);

for i=1:M-1
lambda_1(i)=D(i)-D(i+1);
end
var_k=zeros(1,M-1);
for i=1:M-1
var_k(i)=var(lambda_1(i:M-1));
end
for i=1:M-2
if(var_k(i)~=0)
SORTE(i)=var_k(i+1)/var_k(i);
else
SORTE(i)=1e10;
end
end
[aa num_1]=sort(SORTE(1:M-3));
N=num_1(1);%信源数
Us=V(:,1:N);
Un=V(:,N+1:M);
%%峰值搜索
P_mus=zeros(1,360);
a=zeros(M,1);
angle=[-90:0.5:90-0.5];


for theta_1=1:360%遍历每个角度
    a=exp(-1i*2*pi*d*sin(angle(theta_1)/180*pi)*(0:M-1)');
    %a=exp(-1i*2*pi*d/lambda*sin(angle(theta_1)/180*pi)*(0:M-1)');
    P_mus(theta_1)=1/(a'*Un*Un'*a);
end
figure;
h=plot(angle,10*log10(abs(P_mus)/abs(max(P_mus))));

[num,loc] = findpeaks(10*log10(abs(P_mus)/abs(max(P_mus)))) ;                         %寻找出全部峰值
[a_num,a_loc] = max(num)            ;              %在全部峰值里面找出最大的一个a_num，包含其位置a_loc
location_in_x_1 = loc(a_loc)  ;
fprintf('The angle is %8.5f\n',angle(location_in_x_1))
set(h,'Linewidth',2);
xlabel('入射角/(degree)');
ylabel('归一化空间谱/(dB)');
title('music 算法估计 DOA');
grid on
