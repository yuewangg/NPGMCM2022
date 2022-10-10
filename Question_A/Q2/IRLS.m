functionu=IRLS(Phi,b,p)

epsilon= 1;
u= pinv(Phi)*b;
m= norm(u);
Em=[];
 while epsilon>1e-8
    w=(u.^2+epsilon).^(p/2-1);
    Q  =diag(1./w);
    t = u;
    u = (Q * Phi') *((Phi*Q*Phi') \ b);
    if m<=(sqrt(epsilon)/100)  
         epsilon = epsilon/10;
    end
    m=norm(t-u);
    Em = [Em;m];
 end
epsilon
