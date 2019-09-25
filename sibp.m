function [Weight,b,valaccbatch,vallossbatch,trainaccbatch,trainlossbatch] = sibp(X,Y,maxiter,batchsize,Weight,b,lambda,layersize,ita)
N = length(layersize);
load mnist_val.mat;
disp('A3S Training...')
InitWeight = Weight;
Initb = b;
[~,Xnumber] = size(X);
decay = 0;
Weight_best = Weight;
b_best = b;
%valaccbatch = 0;
for l=1:maxiter
    %perm = randperm(Xnumber/batchsize);
    for set = 1:(Xnumber/batchsize)
        tic;
        decay = decay + 1;
        for set2 = 1:batchsize
            rdset = randi(Xnumber);
            Xin(:,set2) = X(:,rdset);
            Yin(:,set2) = Y(:,rdset);
        end
        %Xin = X;
        %Yin = Y;  %for  full-batch
        %for i=2:N
        %    lambda{i-1} = lambda_0{i-1}/decay;
        %end
        [Weight,b] = sibp_inner(Xin,Yin,Weight,b,lambda,decay,ita);
        [valaccbatch(decay),vallossbatch(decay)] = test(Weight,b,X_val,Y_val);
        [trainaccbatch(decay),trainlossbatch(decay)] = test(Weight,b,X,Y);
        if valaccbatch(decay) == max(valaccbatch)
            Weight_best = Weight;
            b_best = b;
        end
        timetoc = toc;
        fprintf('epoch %d batch %d training accuracy %f validation accuracy %f time cost %f\n',l,set,trainaccbatch(decay),valaccbatch(decay),timetoc);
        %fprintf('epoch %d batch %d training accuracy %f time cost %f\n',l,set,trainaccbatch(decay),timetoc);
    end
    trainacc(l) = trainaccbatch(decay);
    trainloss(l) = trainlossbatch(decay);
end
end