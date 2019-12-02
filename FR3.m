function [x] = FR3(F,B,F_next,W_old,lambda)
x0 = W_old;
maxk=5;
rho=0.6;
sigma=0.1;
k=0; 
epsilon=1e-4; 
mid = max(x0*F+B,0);
midg = mid;
midg(midg>0) = 1;
g = 2*((mid-F_next).*midg)*F' + 2*lambda*(x0 - W_old);
d = -g;
while(k<maxk)
    m=0; 
    mk=0;
    B0_1 = norm(max(x0*F+B,0)-F_next,'fro')^2 + lambda * norm(x0-W_old,'fro')^2;
    gsearch = sigma*sum(sum(g.*d)); 
    m_on = 0;
    while(m<50)   %Armijo
        x1 = x0+rho^m*d;
        A0 = norm(max(x1*F+B,0)-F_next,'fro')^2 + lambda * norm(x1-W_old,'fro')^2;
        %B0 = norm(max(F*x0+B,0)-F_next,'fro')^2 + lambda * norm(x0-W_old,'fro')^2 + gsearch * rho^m; %sigma*rho^m*g'*d
        B0 = B0_1 + gsearch * rho^m;
        if m == 0 && m_on==0
            m = max(floor(log(B0_1*0.01/abs(gsearch))/log(rho)),0);
            m_on = 1;
        end
        if A0 < B0
            mk=m; break;
        end
        m=m+1;
    end
    if mk > 0
        x0=x0+rho^mk*d;
    end
    mid = max(x0*F+B,0);
    midg = mid;
    midg(midg>0) = 1;
    gnew = 2*((mid-F_next).*midg)*F' + 2*lambda*(x0 - W_old);
    beta = norm(gnew,'fro')^2 / norm(g,'fro')^2;
    d = -gnew + beta * d;
    g = gnew;
    k = k + 1;
end
x=x0;
end
