function [Weight,b] = sibp_inner(X,Y,Weight,b,lambda,decay,ita) 
[~,Xnumber] = size(X);
N = length(b)+1;
G{1} = 0;
F{1} = X;
for i=2:N-1
    G{i} = Weight{i-1}*F{i-1}+b{i-1}*ones(1,Xnumber);
    F{i} = lReLu(G{i});
end
F{N} = Weight{N-1}*F{N-1}+b{N-1}*ones(1,Xnumber);
iter = 1;
%ita = 1*0.99^decay; %decay
for j=1:iter
    gd = softmax(F{N}) - Y;
    F{N} =  F{N} - gd * ita;
    Weight{N-1} = FR2(F{N-1},b{N-1}*ones(1,Xnumber),F{N},Weight{N-1},lambda{N-1}(1));
    b{N-1} = FR2(ones(1,Xnumber),Weight{N-1}*F{N-1},F{N},b{N-1},lambda{N-1}(2));
    gd = Weight{N-1}' * gd;
    mid = F{N-1};
    F{N-1} = F{N-1} - ita * gd;
end
for i=(N-2):-1:1
    for j=1:iter
        Weight{i} = FR(F{i},b{i}*ones(1,Xnumber),F{i+1},Weight{i},lambda{i}(1));
        b{i} = FR(ones(1,Xnumber),Weight{i}*F{i},F{i+1},b{i},lambda{i}(2));
        if i > 1
            mid(mid>0) = 1;
            gd = Weight{i}' * (gd .* mid);
            mid = F{i};
            F{i} = F{i} - ita * gd;
        end
    end
end
end