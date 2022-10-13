function [xhat,vr]=OMP(y,A,MaxIt)
%y:measurement
%A:sensing matrix
%t:number of iteration
%k:sparse level

[M,N]=size(A);
xhat=zeros(N,1);
At=zeros(M,MaxIt);
pos_xhat=zeros(1,MaxIt);
r=y;  %Initialization of residual
vr=norm(r);
for i=1:MaxIt
    product=A'*r;  %Inner product of each column in sensing matrix and residual
    [val,pos]=max(abs(product)); %Find the column most relevant to residual
    pos=pos(end);
    At(:,i)=A(:,pos);
    pos_xhat(i)=pos; %save the column index
    
    A(:,pos)=zeros(M,1);
    xhat_ls=(At(:,1:i)'*At(:,1:i))^(-1)*At(:,1:i)'*y;  %Least Square
    r=y-At(:,1:i)*xhat_ls;   %update residual
    
    if (i>2)
        rold= y-At(:,1:i-1)*(At(:,1:i-1)'*At(:,1:i-1))^(-1)*At(:,1:i-1)'*y;
        if (norm(rold-r)<norm(r)*1.5)
            i=i-1;
            break;
        end
    end
    
    vr(i)=norm(r);
end
xhat(pos_xhat(1:i))=xhat_ls(1:i);
xhat(isnan(xhat))=0;