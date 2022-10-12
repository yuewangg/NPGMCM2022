clear; close all;
load("data_q1.mat");
load("data_q2.mat");
n_samples=256;
data_radar=Z_noisy.';
%data_radar=Z.';
n_RX=86;
c=299792458;
N=256;
Q=180;
T=1.25e-7;
fs=1/T;
d=0.0815/85;
f0=78.8e9;
lambda=c/f0;   %雷达信号波长
K=78.986e12;
%% 2维FFT处理
%距离FFT
range_win = hamming(n_samples);   %加海明窗
range_profile =[];
for k=1:n_RX
  temp=data_radar(:,k).*range_win;    %加窗函数
  temp_fft=fft(temp,N);    %对每个chirp做N点FFT
  range_profile(:,k)=temp_fft;
end
disArr=zeros(1,n_RX)
for n=1:n_RX
    [num,loc] = findpeaks( abs(range_profile(:,k))) ;                      %寻找出全部峰值
    [a_num,a_loc] = max(num);                          %在全部峰值里面找出最大的一个a_num，包含其位置a_loc
    dis = c*(((loc(a_loc)-1)*fs)/N)/2/K; 
    disArr(1,n)=dis;
end
%角度FFT
angle_profile = [];
for n=1:N   %range
  temp=range_profile(n,:); 
  %temp=data_radar(n,:); 
  temp_fft=fftshift(fft(temp,Q));    %对2D FFT结果进行Q点FFT
  angle_profile(n,:)=temp_fft;  
end

%% 计算峰值位置
maxPink=0;
maxPinkID=0;
scMaxPink=0;
scMaxPinkID=0;
for n=1:N
    [num,loc] = findpeakm( abs(angle_profile(n,:))) ;                      %寻找出全部峰值
    [a_num,a_loc]=sort(loc,1,'descend');
    if maxPink<a_num(1)
        maxPink=a_num(1);
        maxPinkID=a_loc(1);
        scMaxPink=a_num(2);
        scMaxPinkID=a_loc(2);
    end
end
fprintf('%5.4f   %5.4f\n',maxPink,maxPinkID); 
fprintf('%5.4f   %5.4f\n',scMaxPink,scMaxPinkID); 

angle_profile=abs(angle_profile);
pink=max(angle_profile(:))
[row,pag]=ind2sub(size(angle_profile),find(angle_profile==pink));
[row1,pag1]=ind2sub(size(angle_profile),find(angle_profile==scMaxPink));
%%
figure(1)
f=linspace(0,fs-1,N);
plot(f,abs(range_profile))
title('距离fft转换图')
xlabel('频率/（hz）')
ylabel('信号幅值')
axis([4.2e6 4.5e6 0 1e3]);
grid on
figure(2)
plot(0:1:n_RX-1,disArr,'linewidth',2)
title('各天线计算目标距离图')
xlabel('虚拟天线阵列编号')
ylabel('距离/（m）')
grid on
figure(3)
plot(-90:1:90-1,abs(angle_profile))
title('角度fft转换图')
xlabel('空间频率/（hz）')
ylabel('信号幅值')
axis([-50 +50 0 3.5e4]);
grid on

%% 计算目标距离、速度、角度
fb = ((row-1)*fs)/N ;        %差拍频率
fw = (pag-Q/2-1)/Q;          %空间频率
fb1 = ((row1-1)*fs)/N ;        %差拍频率
fw1 = (pag1-Q/2-1)/Q;          %空间频率
R = c*(fb)/2/K;          %距离公式
theta = asin(fw*lambda/d);  %角度公式
angle = theta*180/pi;
R1 = c*(fb1)/2/K;          %距离公式
theta1 = asin(fw1*lambda/d);  %角度公式
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