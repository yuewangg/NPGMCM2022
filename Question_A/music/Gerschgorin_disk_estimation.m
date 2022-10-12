
function [num]=Gerschgorin_disk_estimation (R,N)
%盖尔圆准则估计信源数目
%input R 协方差矩阵
%input N  快拍数
%output num 信源数
    [row_length col_length]=size(R);
    R_new=R(1:row_length-1,1:col_length-1);
    [V D]=eig(R_new);
    D =diag(D).';
    [D,I0] = sort(D);
    D=fliplr(D);
    V=fliplr(V(:,I0));
    U=[V zeros(row_length-1,1);zeros(1,col_length-1) 1];
    S=U'*R*U;
    [row_lengthS col_lengthS]=size(S);
    P=diag(R)

    k=1;
    while 1
        GDE=abs(S(k,col_length))-1/(2*N*(col_length-1))*sum(abs(S(1:row_length-1,col_length)));%调整因子取为1/N
        if GDE<0 ||k>col_length-2
            break;
        end
        k=k+1;
    end
    num=k-1;
end
