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

MaxARange=55;
MinARange=8;
MaxRRange=15;
MinRRange=1;
Resl=0.001;
AntDis=L/(Na-1);

% %��������
load('../data/data_q4.mat');
Z=Z_antnoisy;

%% �Ľ�ǰ��FFT�㷨����
%�ȶԵ���chirp�ź��ڵ���Ƶ�źŽ���FFT����ת��Ϊ�������
Z2=Z.*(ones(Na,1)*hamming(N)');
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
contourf(1:Na,DisA,log10(abs(Zf1))*20,100,'linestyle','none');
axis([ 1 Na 0 DisA(end)]);
ylabel('����/m');xlabel('���߱��');
title('ÿ�����ߵ������в����');


figure(2);
plot(AngSeq/pi*180,abs(Zf2(PhiInd,DelayInd)));
ylabel('�Ƕ�');xlabel('�ź�ǿ��');
 AngSeq2=(-MaxARange:Resl:MaxARange);


figure(3);
contourf(DisA,AngSeq,(abs(Zf2(PhiInd,:))),100,'linestyle','none');
axis([ 0 DisA(end)  AngSeq(1) AngSeq(end)]);
xlabel('����/m');ylabel('�Ƕ�');

Dic=DicGen(AntDis,lamda,Na,AngSeq2);

DataTmp=Zf1(DelayInd,:).';
%����omp����ϡ���ع�
[Sest,vr]=OMP(DataTmp,Dic,8);
figure(4);clf;hold on;
plot(AngSeq,abs(Zf2(PhiInd,DelayInd)));
plot(AngSeq2,abs(Sest).');
legend('FFT','OMPϡ���ع�')
ylabel('�Ƕ�');xlabel('�ź�ǿ��');

%% �Ľ������ѹ����֪�Ŀ��Ŷ�����
%�ȶԵ���chirp�ź��ڵ���Ƶ�źŽ�ϡ���ع�
RSeq=[0:0.01:MaxRRange];
DicR=DicfGen(gamma,1/Ts,N,RSeq);
[Sest]=OMP(Z(1,:).',DicR,10);
[Ind]=find(abs(Sest)>0);
DataTmp=pinv(DicR(:,Ind))*Z.';
RTmp=RSeq(Ind);

AngSeq2=(-1:Resl:1);
AntLoc=[0:Na-1]*AntDis;
obj=@(x) fobj(x,AntLoc,lamda,Na,AngSeq2,DataTmp);

gbest=pso(obj,Na-1);

 AntLoc2=AntLoc+[0 gbest]*lamda;
for jj=1:length(RTmp)
    Dic=DicGenNew(AntLoc,lamda,Na,AngSeq2);
    %����omp����ϡ���ع�
    [Sest]=OMP(DataTmp(jj,:)',Dic,10);
    [Ind]=find(abs(Sest)>0);
    AngSelOld=AngSeq2(Ind);

    DisRes(jj)=RTmp(jj);
    AngRes(jj)=AngSelOld(1);
end

for jj=1:length(RTmp)
    Dic=DicGenNew(AntLoc2,lamda,Na,AngSeq2);
    %����omp����ϡ���ع�
    [Sest]=OMP(DataTmp(jj,:)',Dic,10);
    [Ind]=find(abs(Sest)>0);
    AngSelOld=AngSeq2(Ind);

    DisRes2(jj)=RTmp(jj);
    AngRes2(jj)=AngSelOld(1);
end
figure(5);clf;hold on;
plot(AntLoc,AntLoc*0,'bo')
plot(AntLoc2,AntLoc2*0,'rx');
legend('��������λ��','���Ƶ�ʵ������λ��');


figure(6);clf;hold on;
plot(DisRes,AngRes,'bo')
plot(DisRes2,AngRes2,'rx');
%legend('����ǰ����λ��','���������λ��');

