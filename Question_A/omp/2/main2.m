clear all;clc;close all;

%�����趨
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

%��������
load('data_q2.mat');

%�ȶԵ���chirp�ź��ڵ���Ƶ�źŽ���FFT����ת��Ϊ�������
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
ylabel('�Ƕ�');xlabel('�ź�ǿ��');

%����ϡ������ֵ� 
AngSeq2=(-55:0.01:55);
Dic=DicGen(AntDis,lamda,Na,AngSeq2);

DataTmp=Zf1(DelayInd,:).';
%����omp����ϡ���ع�
[Sest,vr]=OMP(DataTmp,Dic,8);
figure(2);clf;hold on;
plot(AngSeq,abs(Zf2(PhiInd,DelayInd)));
plot(AngSeq2,abs(Sest).');
legend('FFT','OMPϡ���ع�')
ylabel('�Ƕ�');xlabel('�ź�ǿ��');
