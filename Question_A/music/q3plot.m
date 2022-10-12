load('q3DisCal.mat')
load('q3AngelCal.mat')

Ang=angelResult;
DIS=disResult;
for i=1:32
    if angelResult(i,1)>angelResult(i,2)
        Ang(i,1)=angelResult(i,2);
        Ang(i,2)=angelResult(i,1);
    else
        Ang(i,1)=angelResult(i,1);
        Ang(i,2)=angelResult(i,2);
    end
    if disResult(i,1)>disResult(i,2)
        DIS(i,1)=disResult(i,2);
        DIS(i,2)=disResult(i,1);
    else
        DIS(i,1)=disResult(i,1);
        DIS(i,2)=disResult(i,2);
    end
end

figure(1);
plot(DIS,Ang,'*')
title('二维极坐标散点图')
xlabel('距离/（m）')
ylabel('角度/（°）')
axis([0 10 -50 50]);
grid on
legend('物体1','物体2')        

figure(2);
plot(DIS,Ang)
title('二维极坐标轨迹图')
xlabel('距离/（m）')
ylabel('角度/（°）')
axis([0 10 -50 50]);
grid on
legend('物体1','物体2')     
