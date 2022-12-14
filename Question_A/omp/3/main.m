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
Resl=0.05;
AntDis=L/(Na-1);

% %载入数据
load('data_q3.mat');
% for ii=1:Nf
%     Z=squeeze(Z_time(ii,:,:));
%     %先对单个chirp信号内的中频信号进行FFT处理转换为距离估计
%     Z2=Z.*(ones(Na,1)*hamming(N)');
%     Zf1=fftshift((fft(Z2.')).',2).'/sqrt(N);
%    figure(1);
%    subplot(6,6,ii);
%     DisA=(-N/2:N/2-1)/N/(gamma*Ts/c);
%     contourf(1:Na,DisA,(abs(Zf1)),'linestyle','none');
%     axis([ 1 Na 0 DisA(end)]);
%     ylabel('距离/m');xlabel('天线编号');
%     title('每个天线单独进行测距结果');
% end

%进行逐帧数据处理
for ii=1:Nf
    Z=squeeze(Z_time(ii,:,:));
    %先对单个chirp信号内的中频信号进稀疏重构
    if (ii==1)
        RSeq=[0:0.01:MaxRRange];
    else
        Rseq=[];
        if (ii==2)
            targetROld2=(unique(Res(ii-1).Tdis));
        else
            targetROld2=union(unique(Res(ii-1).Tdis),unique(Res(ii-2).Tdis));
        end
        for jj=1:length(targetROld2)
            Rseq=[Rseq targetROld2(jj)+(-MinRRange:0.01:MinRRange)];
        end
    end
    DicR=DicfGen(gamma,1/Ts,N,RSeq);
    [Sest]=OMP(Z(1,:).',DicR,10);
    [Ind]=find(abs(Sest)>0);
    DataTmp=pinv(DicR(:,Ind))*Z.';
    RTmp=RSeq(Ind);
    
    if (ii==1)
        AngSeq2=(-MaxARange:Resl:MaxARange);
    else
        AngSeq2=[];
        if (ii==2)
            targetAOld=(unique(Res(ii-1).Tang));
        else
            targetAOld=union(unique(Res(ii-1).Tang),unique(Res(ii-2).Tang));
        end
        for jj=1:length(targetAOld)
            AngSeq2=[Rseq targetAOld(jj)+(-MinARange:Resl:MinARange)];
        end
    end
    Res(ii).TdisNum=length(RTmp);
    Res(ii).Tdis=[];
    Res(ii).Tang=[];
    for jj=1:length(RTmp)
        
        Dic=DicGen(AntDis,lamda,Na,AngSeq2);
        %采用omp进行稀疏重构
        [Sest]=OMP(DataTmp(jj,:)',Dic,10);
        [Ind]=find(abs(Sest)>0);
        AngSelOld=AngSeq2(Ind);
        
        Res(ii).AngNum(jj)=length(AngSelOld);
        Res(ii).Tdis=[Res(ii).Tdis RTmp(jj)*ones(1,length(AngSelOld))];
        Res(ii).Tang=[Res(ii).Tang AngSelOld];
    end
end

figure(2);clf;hold on;
for ii=1:Nf
    plot(Res(ii).Tdis,Res(ii).Tang,'rx')
end
ylabel('角度');xlabel('距离');
title('轨迹点');

%进行数据关链
TargetNum=length(Res(1).Tdis);
Target=[];
for ii=1:TargetNum
    Target(ii).Dis(1)=[Res(1).Tdis(ii)];
    Target(ii).Ang(1)=[Res(1).Tang(ii)];
end
DisThr=5;
for ii=2:Nf
    for ii3=1:TargetNum
        Target(ii3).Dis(ii)=Target(ii3).Dis(ii-1);
        Target(ii3).Ang(ii)=Target(ii3).Ang(ii-1);
    end
    for ii2=1:length(Res(ii).Tdis);
        LocNew=[Res(ii).Tdis(ii2) Res(ii).Tang(ii2)];
        Dis=[];
        for ii3=1:TargetNum
            Dis1= abs(Target(ii3).Dis(ii-1)-LocNew(1));
            Dis2= abs(Target(ii3).Ang(ii-1)-LocNew(2));
            
            Dis(ii3)=max(Dis1*2.5,Dis2);
        end
        [DisMin,Ind]=min(Dis);
        if (DisMin<DisThr)
            Target(Ind).Dis(ii)=LocNew(1);
            Target(Ind).Ang(ii)=LocNew(2);
        else
            TargetNum=TargetNum+1;
            for jj=1:ii
            Target(TargetNum).Dis(jj)=LocNew(1);
            Target(TargetNum).Ang(jj)=LocNew(2);
            end
        end
    end
end

figure(3);clf;hold on;
for ii=1:TargetNum
   plot(Target(ii).Dis,Target(ii).Ang);    
end
ylabel('角度');xlabel('距离');
title('轨迹点');
