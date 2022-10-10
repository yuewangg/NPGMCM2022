function [xhat,vr]=OMP(y,A,Mt)

[M,N]=size(A);
xhat=zeros(N,1);
At=zeros(M,Mt);
pxt=zeros(1,Mt);
r=y;  
vr=norm(r);
for i=1:Mt
    product=A'*r;  
    [val,pos]=max(abs(product)); 
    pos=pos(end);
    At(:,i)=A(:,pos);
    pxt(i)=pos; 

    A(:,pos)=zeros(M,1);
    xs=(At(:,1:i)'*At(:,1:i))^(-1)*At(:,1:i)'*y;  
    r=y-At(:,1:i)*xs;  

    if (i>2)
        rold= y-At(:,1:i-1)*(At(:,1:i-1)'*At(:,1:i-1))^(-1)*At(:,1:i-1)'*y;
        if (norm(rold-r)<norm(r)*2)
            i=i-1;
            break;
        end
    end

    vr(i)=norm(r);
end
xhat(pxt(1:i))=xs(1:i);
xhat(isnan(xhat))=0;
