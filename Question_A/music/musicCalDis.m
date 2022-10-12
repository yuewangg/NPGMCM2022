%% music 算法估计 DOA
function [outDis1,outDis2]=musicCalDis(X1)
%%阵列天线参数
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
%%接收信号
%X1=Z_noisy;
%X1=Z;
X2=X1.';
R=X2*X2'/K;
[V,D]=eig(R);%eig 求特征值 D，特征向量 V
Eig=diag(D);%特征值
Eig=sort(Eig,'descend'); %从大到小排特征值
V=fliplr(V);%特征向量对应排
%%信源数估计 SORTE 方法估计出相干信源组一个独立信源看作是一个相干信源组 in[2]
D=diag(D);
D=sort(D,'descend');

N=2;
Us=V(:,1:N);
Un=V(:,N+1:K);
%%峰值搜索
taoMax=1/(T*slope);
rmax=c*taoMax/2;
rmax=10;
P_mus=zeros(1,10000);
a=zeros(K,1);
%distance=[0:rmax/1500:rmax-rmax/1500];
tao=[0:taoMax/10000:taoMax-taoMax/10000];
%tao=[0:taoMax/1000:taoMax-1*taoMax/1000];
for dis_1=1:10000
a=exp(1i*2*pi*slope*tao(dis_1)*T*(0:K-1)');
P_mus(dis_1)=1/(a'*Un*Un'*a);
end
dis=c*tao/2;
%h=plot(dis,10*log10(abs(P_mus)/abs(max(P_mus))));

%[num,loc] = findpeaks(10*log10(abs(P_mus)/abs(max(P_mus)))) ;                        %寻找出全部峰值
%[a_num,a_loc] = max(num);                          %在全部峰值里面找出最大的一个a_num，包含其位置a_loc
%location_in_x_1 = loc(a_loc);  

[num,loc] = findpeakm( 10*log10(abs(P_mus)/abs(max(P_mus))));                       %寻找出全部峰值
[a_num,a_loc]=sort(loc,1,'descend');

maxPink=a_num(1);
scMaxPink=a_num(2);
[pag]=ind2sub(size(10*log10(abs(P_mus)/abs(max(P_mus)))),find(10*log10(abs(P_mus)/abs(max(P_mus)))==maxPink));
[pag1]=ind2sub(size(10*log10(abs(P_mus)/abs(max(P_mus)))),find(10*log10(abs(P_mus)/abs(max(P_mus)))==scMaxPink));

location_in_x_1 = pag; 
location_in_x_2 = pag1; 
dis1=c*tao(location_in_x_1)/2;
dis2=c*tao(location_in_x_2)/2;
fprintf('The dis is %8.5f and %8.5f\n',dis1,dis2);

%set(h,'Linewidth',2);
%xlabel('距离/(m)');
%ylabel('归一化空间谱/(dB)');
%title('music 算法估计距离');
%axis([7.5 9 -30 0]);
%grid on;
outDis1=dis1;
outDis2=dis2;
end