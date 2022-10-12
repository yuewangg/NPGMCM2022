%%ffT代码
clear; close all;
load("data_q1.mat");
load("data_q2.mat");
sampNum=256;
%data_radar=Z_noisy.';
DATArader=Z.';
rxNum=86;
c=299792458;
N=256;
angelNum=180*10;
T=1.25e-7;
fs=1/T;
d=0.0815/85;
f0=78.8e9;
lambda=c/f0;   
K=78.986e12;
%% 2维FFT处理
%距离FFT
range_win = hamming(sampNum);   %加海明窗
RangePro =[];
for k=1:rxNum
  temp=DATArader(:,k).*range_win;    %加窗函数
  temp_fft=fft(temp,N);    %对每个chirp做N点FFT
  RangePro(:,k)=temp_fft;
end
disArr=zeros(1,rxNum)
for n=1:rxNum
    [numN,mloc] = findpeaks( abs(RangePro(:,k))) ;                      
    [aNum,aLoc] = max(numN);                          
    dis = c*(((mloc(aLoc)-1)*fs)/N)/2/K; 
    disArr(1,n)=dis;
end
%角度FFT
AnglePro = [];
for n=1:N   %range
  temp=RangePro(n,:); 
  %temp=data_radar(n,:); 
  temp_fft=fftshift(fft(temp,angelNum));    %对2D FFT结果进行Q点FFT
  AnglePro(n,:)=temp_fft;  
end
%% 计算峰值位置
maxPink=0;
maxPinkID=0;
scMaxPink=0;
scMaxPinkID=0;
for n=1:N
    [numN,mloc] = findpeakm( abs(AnglePro(n,:))) ;                      %寻找出全部峰值
    [aNum,aLoc]=sort(mloc,1,'descend');
    if maxPink<aNum(1)
        maxPink=aNum(1);
        maxPinkID=aLoc(1);
        scMaxPink=aNum(2);
        scMaxPinkID=aLoc(2);
    end
end
fprintf('%5.4f   %5.4f\n',maxPink,maxPinkID); 
fprintf('%5.4f   %5.4f\n',scMaxPink,scMaxPinkID); 
AnglePro=abs(AnglePro);
pink=max(AnglePro(:))
[Row,Pag]=ind2sub(size(AnglePro),find(AnglePro==pink));
[Row1,Pag1]=ind2sub(size(AnglePro),find(AnglePro==scMaxPink));
figure(1)
f=linspace(0,fs-1,N);
dis=c*(f/N)/2/K;
plot(f,abs(RangePro))
title('距离fft转换图')
xlabel('距离/（m）')
ylabel('信号幅值')
axis([3.6e6 3.8e6 0 700]);
grid on
figure(2)
plot(0:1:rxNum-1,disArr,'linewidth',2)
title('各天线计算目标距离图')
xlabel('虚拟天线阵列编号')
ylabel('距离/（m）')
grid on
figure(3)
plot(-90:0.1:90-0.1,abs(AnglePro))
title('角度fft转换图')
xlabel('角度/（°）')
ylabel('信号幅值')
axis([-50 +50 0 3.5e4]);
grid on
%% 计算目标距离、速度、角度
FB = ((Row-1)*fs)/N ;        %差拍频率
FW = (Pag-angelNum/2-1)/angelNum;          %空间频率
FB1 = ((Row1-1)*fs)/N ;        %差拍频率
FW1 = (Pag1-angelNum/2-1)/angelNum;          %空间频率
R = c*(FB)/2/K;          %距离公式
theta = asin(FW*lambda/d);  %角度公式
angle = theta*180/pi;
R1 = c*(FB1)/2/K;          %距离公式
theta1 = asin(FW1*lambda/d);  %角度公式
angle1 = theta1*180/pi;
fprintf('目标距离： %f m %f m\n',R,R1);
fprintf('目标角度： %f° %f°\n',angle,angle1);
targetDis=[R,R1];
targetAng=[angle,angle1];
figure(4)
plot(targetDis,targetAng,'*')
title('二维极坐标图')
xlabel('距离/（m）')
ylabel('角度/（°）')
axis([0 10 -50 50]);
grid on
