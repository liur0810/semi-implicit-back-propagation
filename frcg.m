function [x]=frcg(W_old,F,B,F_next,lambda)
W_old = W_old';
F = F';
B = B';
F_next = F_next';
x0 = W_old;
maxk=3;   %����������
rho=0.6;
sigma=0.4;
k=0;  
epsilon=1e-4; 
n=length(x0);
while(k<maxk)
    %g=feval(funs,x0);  %�����ݶ�
    mid = max(F*x0+B,0);
    midg = mid;
    midg(midg>0) = 1;
    g = 2*F'*((mid-F_next).*midg) + 2*lambda*(x0 - W_old);
    itern=k-(n+1)*floor(k/(n+1));
    itern=itern+1;
    %������������
    if(itern==1)  
        d=-g;  
    else
        beta=(g'*g)/(g0'*g0);
        d=-g+beta*d0;
        gd=g'*d;
        if(gd>=0.0)
            d=-g;  
        end
    end
    if(norm(g)<epsilon), break; end   %������ֹ����
    m=0; mk=0;
    B0_1 = max(F*x0+B,0)-F_next;
    B0_2 = x0-W_old;
    B0_3 = norm(B0_1,'fro')^2 + lambda * norm(B0_2,'fro')^2;
    while(m<50)   %Armijo����
        x1 = x0+rho^m*d;
        gsearch = sigma*g'*d;
        A0_1 = max(F*x1+B,0)-F_next;
        A0_2 = x1-W_old;
        A0 = norm(A0_1,'fro')^2 + lambda * norm(A0_2,'fro')^2;
        %A0 = norm(max(F*x1+B,0)-F_next,'fro')^2 + lambda * norm(x1-W_old,'fro')^2;
        %B0 = norm(max(F*x0+B,0)-F_next,'fro')^2 + lambda * norm(x0-W_old,'fro')^2 + gsearch * rho^m; %sigma*rho^m*g'*d
        B0 = B0_3 + gsearch * rho^m;
        %if(feval(fun,x0+rho^m*d)<feval(fun,x0)+sigma*rho^m*g'*d)
        if A0 < B0
            mk=m; break;
        end
        m=m+1;
    end
    if mk > 0
        x0=x0+rho^mk*d;
    end
    g0=g;  
    d0=d; 
    k=k+1;
end
x=x0;
x=x';
end