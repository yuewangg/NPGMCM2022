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
load('data_q1.mat');

%�ȶԵ���chirp�ź��ڵ���Ƶ�źŽ���FFT����ת��Ϊ�������
Z2=Z.*(ones(Na,1)*hamming(N)');
Zf1=fftshift((fft(Z2.')).',2).';

DisA=(-N/2:N/2-1)/N/(gamma*Ts/c);

figure(1);
contourf(1:Na,DisA,log10(abs(Zf1))*20,100,'linestyle','none');
axis([ 1 Na 0 DisA(end)]);
ylabel('����/m');xlabel('���߱��');
title('ÿ�����ߵ������в����');


[~,DelayInd]=max(sum(abs(Zf1),2));
%��86�����ߵ�chirp�ź�
Nfft=Na*16;
Zf2=fftshift((fft(Zf1.',Nfft)),1);

AntDis=L/(Na-1);
PhiSeq=(-Nfft/2:Nfft/2-1)/Nfft/(AntDis/lamda);
PhiInd=find(abs(PhiSeq)<1);
AngSeq=asin(PhiSeq(PhiInd));

figure(2);
plot(AngSeq/pi*180,abs(Zf2(PhiInd,DelayInd)));
ylabel('�Ƕ�');xlabel('�ź�ǿ��');

figure(3);
contourf(DisA,AngSeq,(abs(Zf2(PhiInd,:))),100,'linestyle','none');
axis([ 0 DisA(end)  AngSeq(1) AngSeq(end)]);
xlabel('����/m');ylabel('�Ƕ�');




