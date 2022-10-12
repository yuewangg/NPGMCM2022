function  calAngAndDis(X1,NUMk)

M=86;% 阵列天线数量
K=256;%快拍数
N=123;%信号源数目
f0=78.8e9;%射频
c=299792458;
slope=78.986e12;
lambda=c/f0;%波长
T=1.25e-7;
Fs=1/T;%采样频率
d=0.0815/85;%阵元间距
X2=reshape(X1,M*K,1);
%X1=Z;
R=X2*X2'/K;
%R=X1*X1'/K;
Nfb=30;
L=226;
Mfb=M*Nfb;

Rfb=zeros(Mfb,Mfb);
for i=0:1:L
    Rii=R(M*i+1:Mfb+M*i,M*i+1:Mfb+M*i);
    Rfb=Rfb+Rii+Rii';
end
R=Rfb/2*L;
[V,D]=eig(R);%eig 求特征值 D，特征向量 V
Eig=diag(D);%特征值
Eig=sort(Eig,'descend'); %从大到小排特征值
V=fliplr(V);%特征向量对应排
N=2;
Us=V(:,1:N*M);
Un=V(:,N*M+1:Mfb);
%%峰值搜索
taoMax=1/(T*slope);
rmax=c*taoMax/2;

angle=[-20:1:20];
tao=[500*taoMax/1518:taoMax/1518:670*taoMax/1518];
lengthAng=length(angle);
lengthDis=length(tao);
P_mus=zeros(lengthAng,lengthDis);
a=zeros(M,1);
r=zeros(Nfb,1);
dis=c*tao/2;
count=1.0;
count=count-1;
for theta_1=1:lengthAng%遍历每个角度
    fprintf('第%d次 进度百分之%.2f\n',NUMk,count/(lengthDis*lengthAng))
    for dis_1=1:lengthDis
        count=count+1;
        %a=exp(1i*2*pi*sin(angle(theta_1))*(0:M-1)');
        a=exp(-1i*2*pi*d/lambda*sin(angle(theta_1)/180*pi)*(0:M-1)');
        r=exp(1i*2*pi*slope*tao(dis_1)*T*(0:Nfb-1)');
        kronResult=kron(r,a);
        P_mus(theta_1,dis_1)=1/(kronResult'*Un*Un'*kronResult);
    end
end
pus=10*log10(abs(P_mus)/abs(max(max(P_mus))));
figure;
[X,Y]=meshgrid(dis,angle);
mesh(X,Y,pus)
title('距离角度二维估计music空间谱')
xlabel('距离/（m）')
ylabel('角度/（°）')
zlabel('幅度/（db）')

end

