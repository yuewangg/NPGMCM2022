function DicR=DfG(gamma,Fs,Number,D)

DicR=zeros(Number,length(D));
c=299792458;

for jj=1:length(D)
    DicR0=exp(1i*2*pi*(0:Number-1)*gamma/Fs*D(jj)/c);
    
    DicR(:,jj)=DicR0;
end