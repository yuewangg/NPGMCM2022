%% music 算法估计 DOA
clear all;clc;close all ;
load("data_q1.mat");
load("data_q2.mat");
load("data_q3.mat");
load("data_q4.mat");
%%阵列天线参数
M=86;
Knum=256;
fF0=78.8e9;
c=299792458;
slope=78.986e12;
lambda=c/fF0;
T=1.25e-7;
Fs=1/T;
d=0.0815/85;
X1=Z_noisy;
X2=X1.';
R=X2*X2'/Knum;
[vetor,D]=eig(R);
eigG=diag(D);%特征值
eigG=sort(eigG,'descend'); 
vetor=fliplr(vetor);
N=2;
Us=vetor(:,1:N);
Un=vetor(:,N+1:Knum);
taoMax=1/(T*slope);
rmax=c*taoMax/2;
P_mus=zeros(1,10000);
a=zeros(Knum,1);
%distance=[0:rmax/1500:rmax-rmax/1500];
tao=[0:taoMax/10000:taoMax-taoMax/10000];
%tao=[0:taoMax/1000:taoMax-1*taoMax/1000];
for dis_1=1:10000
a=exp(1i*2*pi*slope*tao(dis_1)*T*(0:Knum-1)');
P_mus(dis_1)=1/(a'*Un*Un'*a);
end
figure;
dis=c*tao/2;
h=plot(dis,10*log10(abs(P_mus)/abs(max(P_mus))));
[num,loc] = findpeakm( 10*log10(abs(P_mus)/abs(max(P_mus))));                       %寻找出全部峰值
[a_num,a_loc]=sort(loc,1,'descend');
maxPink=a_num(1);
scMaxPink=a_num(2);
[PAGN]=ind2sub(size(10*log10(abs(P_mus)/abs(max(P_mus)))),find(10*log10(abs(P_mus)/abs(max(P_mus)))==maxPink));
[PAGN1]=ind2sub(size(10*log10(abs(P_mus)/abs(max(P_mus)))),find(10*log10(abs(P_mus)/abs(max(P_mus)))==scMaxPink));
L1 = PAGN; 
L2 = PAGN1; 
dis1=c*tao(L1)/2;
dis2=c*tao(L2)/2;
fprintf('The dis is %8.5f and %8.5f\n',dis1,dis2);
set(h,'Linewidth',2);
xlabel('距离/(m)');
ylabel('归一化空间谱/(dB)');
title('music 算法估计距离');
axis([7.5 9 -30 0]);
grid on;