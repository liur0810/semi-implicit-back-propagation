function [WW] = FR(F,B,F_next,W_old,lambda)
[n,m] = size(W_old);
WW = zeros([n,m]);
for i = 1:n
    W{i} = zeros(1,m);
end
for i = 1:n
    W{i} = frcg(W_old(i,:),F,B(i,:),F_next(i,:),lambda);
    WW(i,:) = W{i};
end
end

