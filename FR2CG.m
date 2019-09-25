function [WW] = FR2CG(F,B,F_next,W_old,lambda)
[n,m] = size(W_old);
WW = zeros([n,m]);
for i = 1:n
    W{i} = zeros(1,m);
end
A0 = F * F';
A0 = eye(size(A0)) + A0;
A0 = A0';
B0 = (F_next-B)*F' + lambda * W_old; % W*A0 = B0
B0 = B0';
for i = 1:n
    W{i} = cg(A0',B0(:,i),W_old(i,:)');
    WW(i,:) = W{i};
end
end

