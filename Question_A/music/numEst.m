clear all;clc;close all ;
load("data_q1.mat");
load("data_q2.mat");
%%阵列天线参数
M=86;% 阵列天线数量
N=256;%快拍数

f0=78.8e9;%射频
c=299792458;

lambda=c/f0;%波长
T=1.25e-7;
Fs=1/T;%采样频率
d=0.0815/255;%阵元间距
%%接收信号
X2=Z_noisy;
X1=Z;


R=X1*X1'/N;
R2=X2*X2'/N;

[row_length col_length]=size(R2);
R_new=R2(1:row_length-1,1:col_length-1);
[V D]=eig(R_new);
D =diag(D).';
[D,I0] = sort(D);
D=fliplr(D);
V=fliplr(V(:,I0));
U=[V zeros(row_length-1,1);zeros(1,col_length-1) 1];
S=U'*R*U;
[row_lengthS col_lengthS]=size(S);
P=diag(R);

k=1;
while 1
    GDE=abs(S(k,col_length))-1/(2*N*(col_length-1))*sum(abs(S(1:row_length-1,col_length)));%调整因子取为1/N
    if GDE<0 ||k>col_length-2
        break;
    end
    k=k+1;
end
num=k-1
    
