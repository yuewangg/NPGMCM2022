function Dic=DG(D,Lamda,Number,A)

Dic=zeros(Number,length(A));

for jj=1:length(A)
    Dic0=exp(1i*2*pi*(0:Number-1)*D/Lamda*sind(A(jj)));
    
    Dic(:,jj)=Dic0;
end