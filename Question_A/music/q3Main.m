load("data_q3.mat");

X1=zeros(86,256);
angelResult=zeros(32,2);
disResult=zeros(32,2);
for k=1:32
    for i=1:86
       for j=1:256
            X1(i,j)=Z_time(k,i,j);
        end
    end
    calAngAndDis(X1,k);
%    [angel1,angel2]=musicCalAng(X1);
%    [dis1,dis2]=musicCalDis(X1);
%    angelResult(k,1)=angel1;
%    angelResult(k,2)=angel2;
%    disResult(k,1)=dis1;
%    disResult(k,2)=dis2;
end
