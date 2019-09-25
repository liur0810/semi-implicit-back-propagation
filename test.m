function [acc,loss] = test(Weight,b,X,Y)
ls = X;
N = length(b) + 1;
[~,Xnumber] = size(X);
for i=1:(N-2)
    ls = lReLu(Weight{i}*ls + b{i}*ones(1,Xnumber));
end
ls = softmax(Weight{N-1}*ls + b{N-1}*ones(1,Xnumber));
[~,m1] = max(ls);
[~,m2] = max(Y);
m1 = m1 -m2;
m1(m1~=0) = 1;
calloss = sum(sum(Y.*log(ls)));
acc = (Xnumber-sum(m1))/Xnumber;
loss = -calloss/Xnumber;