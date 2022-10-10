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

MaxARange=55;
MinARange=8;
MaxRRange=15;
MinRRange=1;
Resl=0.001;
AntDis=L/(Na-1);

% %载入数据
load('../data/data_q4.mat');
Z=Z_antnoisy;

%% 改进前的FFT算法处理
%先对单个chirp信号内的中频信号进行FFT处理转换为距离估计
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
ylabel('距离/m');xlabel('天线编号');
title('每个天线单独进行测距结果');


figure(2);
plot(AngSeq/pi*180,abs(Zf2(PhiInd,DelayInd)));
ylabel('角度');xlabel('信号强度');
 AngSeq2=(-MaxARange:Resl:MaxARange);


figure(3);
contourf(DisA,AngSeq,(abs(Zf2(PhiInd,:))),100,'linestyle','none');
axis([ 0 DisA(end)  AngSeq(1) AngSeq(end)]);
xlabel('距离/m');ylabel('角度');

Dic=DicGen(AntDis,lamda,Na,AngSeq2);

DataTmp=Zf1(DelayInd,:).';
%采用omp进行稀疏重构
[Sest,vr]=OMP(DataTmp,Dic,8);
figure(4);clf;hold on;
plot(AngSeq,abs(Zf2(PhiInd,DelayInd)));
plot(AngSeq2,abs(Sest).');
legend('FFT','OMP稀疏重构')
ylabel('角度');xlabel('信号强度');

%% 改进后基于压缩感知的抗扰动天线
%先对单个chirp信号内的中频信号进稀疏重构
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
    %采用omp进行稀疏重构
    [Sest]=OMP(DataTmp(jj,:)',Dic,10);
    [Ind]=find(abs(Sest)>0);
    AngSelOld=AngSeq2(Ind);

    DisRes(jj)=RTmp(jj);
    AngRes(jj)=AngSelOld(1);
end

for jj=1:length(RTmp)
    Dic=DicGenNew(AntLoc2,lamda,Na,AngSeq2);
    %采用omp进行稀疏重构
    [Sest]=OMP(DataTmp(jj,:)',Dic,10);
    [Ind]=find(abs(Sest)>0);
    AngSelOld=AngSeq2(Ind);

    DisRes2(jj)=RTmp(jj);
    AngRes2(jj)=AngSelOld(1);
end
figure(5);clf;hold on;
plot(AntLoc,AntLoc*0,'bo')
plot(AntLoc2,AntLoc2*0,'rx');
legend('理想天线位置','估计的实际天线位置');


figure(6);clf;hold on;
plot(DisRes,AngRes,'bo')
plot(DisRes2,AngRes2,'rx');
%legend('纠正前估计位置','纠正后估计位置');

