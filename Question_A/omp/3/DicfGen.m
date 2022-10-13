function Dic=DicfGen(gamma,Fs,Num,R)

Dic=zeros(Num,length(R));

for ii=1:length(R)
    Dic0=exp(1i*2*pi*(0:Num-1)*gamma/Fs*R(ii)/3e8);
    
    Dic(:,ii)=Dic0;
end