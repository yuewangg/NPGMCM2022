function DicA=DG(Dis,Lamda,Number,A)

DicA=zeros(Number,length(A));

for jj=1:length(A)
    DicA0=exp(1i*2*pi*(0:Number-1)*Dis/Lamda*sind(A(jj)));
    
    DicA(:,jj)=DicA0;
end