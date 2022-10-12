
%% music 算法估计 DOA
clear all;clc;close all ;
load("data_q1.mat");
load("data_q2.mat");
load("data_q3.mat");
%%阵列天线参数
M=86;% 阵列天线数量
K=256;%快拍数
F0=78.8e9;%射频
c=299792458;

Lamda=c/F0;%波长
T=1.25e-7;
Fs=1/T;%采样频率
d=0.0815/85;%阵元间距
%%接收信号
X1=Z_noisy;
%X1=Z;
R=X1*X1'/K;
[Vettor,D]=eig(R);%eig 求特征值 D，特征向量 V
eEig=diag(D);%特征值
eEig=sort(eEig,'descend'); %从大到小排特征值
Vettor=fliplr(Vettor);%特征向量对应排
N=2;
Us=Vettor(:,1:N);
UN=Vettor(:,N+1:M);
%%峰值搜索
Pmuss=zeros(1,2000);
a=zeros(M,1);
angle=[-50:0.05:50-0.05];
e1=zeros(M,1);
e1(1,1)=1;
for theta_1=1:2000%遍历每个角度
a=exp(-1i*2*pi*d/Lamda*sin(angle(theta_1)/180*pi)*(0:M-1)');
Pmuss(theta_1)=1/(a'*UN*UN'*a);
%P_mus(theta_1)=1/abs(a'*Un*Un'*e1)^2;
end
figure;
hhh=plot(angle,10*log10(abs(Pmuss)/abs(max(Pmuss))));
[num,loc] = findpeakm( 10*log10(abs(Pmuss)/abs(max(Pmuss)))) ;                      %寻找出全部峰值
[a_num,a_loc]=sort(loc,1,'descend')

maxPink=a_num(1);
scMaxPink=a_num(2);
[pag]=ind2sub(size(10*log10(abs(Pmuss)/abs(max(Pmuss)))),find(10*log10(abs(Pmuss)/abs(max(Pmuss)))==maxPink));
[pag1]=ind2sub(size(10*log10(abs(Pmuss)/abs(max(Pmuss)))),find(10*log10(abs(Pmuss)/abs(max(Pmuss)))==scMaxPink));

location_in_x_1 = pag; 
location_in_x_2 = pag1; 

fprintf('The angle is %8.5f and %8.5f\n',angle(location_in_x_1),angle(location_in_x_2));

set(hhh,'Linewidth',2);
xlabel('入射角/(degree)');
ylabel('归一化空间谱/(dB)');
title('music 算法估计 DOA');
%axis([-50 50 -4 0]);
grid on;

