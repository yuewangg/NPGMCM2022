function Dic=DicGenNew(AntLoc,Lamda,Num,Ang)

Dic=zeros(Num,length(Ang));

for ii=1:length(Ang)
    Dic0=exp(1i*2*pi*AntLoc/Lamda*sind(Ang(ii)));
    
    Dic(:,ii)=Dic0;
end