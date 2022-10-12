
%% music 算法估计 DOA
clear all;clc;close all ;
load("data_q1.mat");
load("data_q2.mat");
load("data_q4.mat");
Mmunm=86;% 阵列天线数量
Kuaip=256;%快拍数
shef=78.8e9;%射频
c=299792458;
slope=78.986e12;
lambda=c/shef;
T=1.25e-7;
Fs=1/T;
d=0.0815/85;
X1=Z;
X2=reshape(Z_noisy,Mmunm*Kuaip,1);
R=X2*X2'/Kuaip;
Nfb=30;
L=226;
Mfb=Mmunm*Nfb;
IDi=Mmunm*Z+1;
Rfb=zeros(Mfb,Mfb);
for i=0:1:L
    Rii=R(Mmunm*i+1:Mfb+Mmunm*i,Mmunm*i+1:Mfb+Mmunm*i);
    Rfb=Rfb+Rii+Rii';
end
R=Rfb/2*L;
[vectorV,D]=eig(R);
eigE=diag(D);
eigE=sort(eigE,'descend'); 
vectorV=fliplr(vectorV);
numN=2;
Us=vectorV(:,1:numN*Mmunm);
Un=vectorV(:,numN*Mmunm+1:Mfb);
%%峰值搜索
taoMax=1/(T*slope);
rmax=c*taoMax/2;
rmax=10;
angle=[-10:0.5:10];
tao=[600*taoMax/1518:taoMax/1518:900*taoMax/1518];
lengthAng=length(angle);
lengthDis=length(tao);
ppPums=zeros(lengthAng,lengthDis);
a=zeros(Mmunm,1);
r=zeros(Nfb,1);
dis=c*tao/2;
count=1.0;
count=count-1;
for theta_1=1:lengthAng
    fprintf('百分之%.2f\n',count/(lengthDis*lengthAng))
    for dis_1=1:lengthDis
        count=count+1;
        %a=exp(1i*2*pi*sin(angle(theta_1))*(0:M-1)');
        a=exp(-1i*2*pi*d/lambda*sin(angle(theta_1)/180*pi)*(0:Mmunm-1)');
        r=exp(1i*2*pi*slope*tao(dis_1)*T*(0:Nfb-1)');
        kronResult=kron(r,a);
        ppPums(theta_1,dis_1)=1/(kronResult'*Un*Un'*kronResult);
    end
end
ppus=10*log10(abs(ppPums)/abs(max(max(ppPums))));
figure;
[X,Y]=meshgrid(dis,angle);
mesh(X,Y,ppus)
title('距离角度二维估计music空间谱')
xlabel('距离/（m）')
ylabel('角度/（°）')
zlabel('幅度/（db）')


