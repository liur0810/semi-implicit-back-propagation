function [W] = lReLu(W)
W1 = W;
W2 = W;
W1(W1<0) = 0;
%W2(W2>0) = 0;
W = W1;
end

