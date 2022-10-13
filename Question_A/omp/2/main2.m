clear all;clc;close all;

%参数设定
Ts=1.25e-7;
T=3.2E-5;

N=T/Ts;
Nf=32;
L=0.0815;
gamma=78.986e12;
Na=86;
f0=78.8e9;
c=3e8;
lamda=c/f0;

%载入数据
load('data_q2.mat');

%先对单个chirp信号内的中频信号进行FFT处理转换为距离估计
Z2=Z_noisy.*(ones(Na,1)*hamming(N)');
Zf1=fftshift((fft(Z2.')).',2).'/sqrt(N);

DisA=(-N/2:N/2-1)/N/(gamma*Ts/c);
[~,DelayInd]=max(sum(abs(Zf1),2));

Nfft=Na*64;
Zf2=fftshift((fft(Zf1.',Nfft)),1)/sqrt(Nfft);
AntDis=L/(Na-1);
PhiSeq=(-Nfft/2:Nfft/2-1)/Nfft/(AntDis/lamda);
PhiInd=find(abs(PhiSeq)<1);
AngSeq=asin(PhiSeq(PhiInd))/pi*180;

figure(1);
plot(AngSeq,abs(Zf2(PhiInd,DelayInd)));
ylabel('角度');xlabel('信号强度');

%构造稀疏矩阵字典 
AngSeq2=(-55:0.01:55);
Dic=DicGen(AntDis,lamda,Na,AngSeq2);

DataTmp=Zf1(DelayInd,:).';
%采用omp进行稀疏重构
[Sest,vr]=OMP(DataTmp,Dic,8);
figure(2);clf;hold on;
plot(AngSeq,abs(Zf2(PhiInd,DelayInd)));
plot(AngSeq2,abs(Sest).');
legend('FFT','OMP稀疏重构')
ylabel('角度');xlabel('信号强度');
