function  calAngAndDis(X1,NUMk)

M=86;% ������������
K=256;%������
N=123;%�ź�Դ��Ŀ
f0=78.8e9;%��Ƶ
c=299792458;
slope=78.986e12;
lambda=c/f0;%����
T=1.25e-7;
Fs=1/T;%����Ƶ��
d=0.0815/85;%��Ԫ���
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
[V,D]=eig(R);%eig ������ֵ D���������� V
Eig=diag(D);%����ֵ
Eig=sort(Eig,'descend'); %�Ӵ�С������ֵ
V=fliplr(V);%����������Ӧ��
N=2;
Us=V(:,1:N*M);
Un=V(:,N*M+1:Mfb);
%%��ֵ����
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
for theta_1=1:lengthAng%����ÿ���Ƕ�
    fprintf('��%d�� ���Ȱٷ�֮%.2f\n',NUMk,count/(lengthDis*lengthAng))
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
title('����Ƕȶ�ά����music�ռ���')
xlabel('����/��m��')
ylabel('�Ƕ�/���㣩')
zlabel('����/��db��')

end

