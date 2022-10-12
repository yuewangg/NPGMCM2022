
function [num]=AIC(R,N,M,D)
%AIC准则估计信源数目
%input R 阵元接收信号的协方差矩阵
%input N 采样数
%input M 阵元数
%input D 1：一维线阵 2：二维方阵
%output num 信源数目
if D==1
    [u,s,v] = svd(R);
    sd = diag(s);
    sd=sort(sd,'descend');
    a = zeros(1,M);
    for m = 0:M-1
        negv = sd(m+1:M);
        Tsph = mean(negv)/((prod(negv))^(1/(M-m)));
        a(m+1) = 2*N*(M-m)*log(Tsph) + m*(2*M-m)*2;
    end
    [y,b] = min(a);
    num = b - 1;
else
    [u,s,v] = svd(R);
    sd = diag(s);
    sd=sort(sd,'descend');
    a = zeros(1,prod(M));
    for m = 0:prod(M)-1
        negv = sd(m+1:prod(M));
        Tsph = mean(negv)/((prod(negv))^(1/(prod(M)-m)));
        a(m+1) = N*(prod(M)-m)*log(Tsph) + m*(2*prod(M)-m);
    end
    [y,b] = min(a);
    num = b - 1;
end
end
