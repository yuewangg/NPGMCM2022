function f=fobj(x,AntLoc,lamda,Na,AngSeq2,DataTmp)


AntLoc=AntLoc+[0 x]*lamda;
Dic=DicGenNew(AntLoc,lamda,Na,AngSeq2);

for jj=1:size(DataTmp,1)
    %采用omp进行稀疏重构
    [Sest,vr]=OMP(DataTmp(jj,:)',Dic,10);
    res(jj)=norm(vr);
end
f=max(res);