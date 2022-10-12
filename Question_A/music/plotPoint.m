close all; clc

targetDis=[8.184077,8.184077];
targetAng=[-1.263110,1.263110];
MUSDis=[8.198,8.28637];
MUSAng=[-1.745,1.695];
figure(4)
plot(targetDis,targetAng,'*')
hold on;
plot(MUSDis,MUSAng,'+')
title('二维极坐标图')
xlabel('距离/（m）')
ylabel('角度/（°）')
legend('fft方法','music方法')
grid on

