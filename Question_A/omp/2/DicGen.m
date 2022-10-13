function Dic=DicGen(Dis,Lamda,Num,Ang)

Dic=zeros(Num,length(Ang));

for ii=1:length(Ang)
    Dic0=exp(1i*2*pi*(0:Num-1)*Dis/Lamda*sind(Ang(ii)));
    
    Dic(:,ii)=Dic0;
end