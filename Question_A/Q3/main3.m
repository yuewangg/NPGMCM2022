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

MaxARange=55;
MinARange=8;
MaxRRange=15;
MinRRange=1;
Resl=0.05;
AntDis=L/(Na-1);

load('../data/data_q3.mat');

for jj=1:Nf
    Z=squeeze(Z_time(jj,:,:));

    if (jj==1)
        RSeq=[0:0.01:MaxRRange];
    else
        Rseq=[];
        if (jj==2)
            targetROld2=(unique(Res(jj-1).Tdis));
        else
            targetROld2=union(unique(Res(jj-1).Tdis),unique(Res(jj-2).Tdis));
        end
        for kk=1:length(targetROld2)
            Rseq=[Rseq targetROld2(kk)+(-MinRRange:0.01:MinRRange)];
        end
    end
    DicR=DfG(gamma,1/Ts,N,RSeq);
    [Sest]=OMP(Z(1,:).',DicR,10);
    [Ind]=find(abs(Sest)>0);
    DataTmp=pinv(DicR(:,Ind))*Z.';
    RTmp=RSeq(Ind);

    if (jj==1)
        AngSeq2=(-MaxARange:Resl:MaxARange);
    else
        AngSeq2=[];
        if (jj==2)
            targetAOld=(unique(Res(jj-1).Tang));
        else
            targetAOld=union(unique(Res(jj-1).Tang),unique(Res(jj-2).Tang));
        end
        for kk=1:length(targetAOld)
            AngSeq2=[Rseq targetAOld(kk)+(-MinARange:Resl:MinARange)];
        end
    end
    Res(jj).TdisNum=length(RTmp);
    Res(jj).Tdis=[];
    Res(jj).Tang=[];
    for kk=1:length(RTmp)

        Dic=DG(AntDis,lamda,Na,AngSeq2);

        [Sest]=OMP(DataTmp(kk,:)',Dic,10);
        [Ind]=find(abs(Sest)>0);
        AngSelOld=AngSeq2(Ind);

        Res(jj).AngNum(kk)=length(AngSelOld);
        Res(jj).Tdis=[Res(jj).Tdis RTmp(kk)*ones(1,length(AngSelOld))];
        Res(jj).Tang=[Res(jj).Tang AngSelOld];
    end
end

figure(2);clf;hold on;
for jj=1:Nf
    plot(Res(jj).Tdis,Res(jj).Tang,'rx')
end
ylabel('A');xlabel('D');
title('1');

TargetNum=length(Res(1).Tdis);
Target=[];
for jj=1:TargetNum
    Target(jj).Dis(1)=[Res(1).Tdis(jj)];
    Target(jj).Ang(1)=[Res(1).Tang(jj)];
end
DisThr=5;
for jj=2:Nf
    for jj3=1:TargetNum
        Target(jj3).Dis(jj)=Target(jj3).Dis(jj-1);
        Target(jj3).Ang(jj)=Target(jj3).Ang(jj-1);
    end
    for jj2=1:length(Res(jj).Tdis);
        LocNew=[Res(jj).Tdis(jj2) Res(jj).Tang(jj2)];
        Dis=[];
        for jj3=1:TargetNum
            Dis1= abs(Target(jj3).Dis(jj-1)-LocNew(1));
            Dis2= abs(Target(jj3).Ang(jj-1)-LocNew(2));

            Dis(jj3)=max(Dis1*2.5,Dis2);
        end
        [DisMin,Ind]=min(Dis);
        if (DisMin<DisThr)
            Target(Ind).Dis(jj)=LocNew(1);
            Target(Ind).Ang(jj)=LocNew(2);
        else
            TargetNum=TargetNum+1;
            for kk=1:jj
            Target(TargetNum).Dis(kk)=LocNew(1);
            Target(TargetNum).Ang(kk)=LocNew(2);
            end
        end
    end
end

figure(3);clf;hold on;
for jj=1:TargetNum
   plot(Target(jj).Dis,Target(jj).Ang);
end
ylabel('A');xlabel('D');
title('2');
