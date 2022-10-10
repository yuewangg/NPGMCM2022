clear all;clc;close all;

Ts=1.25e-7;
T=3.2E-5;

N=T/Ts;
Nf=32;
L=0.0815;
gamma=78.986e12;
Na=86;
f0=78.8e9;
c=299792458;
lamda=c/f0;

load('../data/data_q4.mat');

Z=Z_antnoisy;

Z2=Z.*(ones(Na,1)*hamming(N)');
Zf1=fftshift((fft(Z2.')).',2).'/sqrt(N);

D=(-N/2:N/2-1)/N/(gamma*Ts/c);
[~,DelayInd]=max(sum(abs(Zf1),2));

Nfft=Na*64;
Nf=fftshift((fft(Zf1.',Nfft)),1)/sqrt(Nfft);
AntDis=L/(Na-1);
Qes=(-Nfft/2:Nfft/2-1)/Nfft/(AntDis/lamda);
Pd=find(abs(Qes)<1);
Aq=asin(Qes(Pd))/pi*180;

figure(1);
plot(Aq,abs(Nf(Pd,DelayInd)));
xlabel('Angle');ylabel('RSSI');

Aq2=(-55:0.01:55);
Dic=DG(AntDis,lamda,Na,Aq2);

DataTmp=Zf1(DelayInd,:).';

[Sest,vr]=IRLS(DataTmp,Dic,8);
r_fft=abs(Nf(Pd,DelayInd));
r_omp=abs(Sest).';
[Ind]=find(abs(Sest)>0);
RSeq=[0:0.00136:15];
RTmp=RSeq(Ind);
[xpeak,locs] = findallpeaks(Aq,r_fft,12,1.5);
[xpeak2,locs2] = findallpeaks(Aq2,r_omp,12,1.5);
aaa=Aq(locs);
aab=Aq2(locs2);
figure(2);clf;hold on;
plot(Aq,r_fft);
plot(Aq2,r_omp);
plot(Aq(locs),r_fft(locs),'*r');
plot(Aq2(locs2),r_omp(locs2),'*g');
legend('MUSIC','OMP');
xlabel('Angle');ylabel('RSSI');
