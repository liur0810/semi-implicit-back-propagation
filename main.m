clearvars;
load mnist_train.mat; %cifar10_train.mat
[Xsize,Xnumber] = size(X);
[Ysize,~] = size(Y);
layersize = [Xsize,500,Ysize]; 
N = length(layersize);
for i=2:N
    Weight{i-1} =  normrnd(0,0.01,[layersize(i) layersize(i-1)]);
    b{i-1} = normrnd(0,0.01,[layersize(i) 1]);
    lambda{i-1} = [1,1,1,1];
end
InitWeight = Weight;
Initb = b;
batchsize = 100;
maxiter = 2;
ita = 0.1;
[Weight,b,valaccbatch,valloss,trainaccbatch,trainloss] = sibp(X,Y,maxiter,batchsize,InitWeight,Initb,lambda,layersize,ita);
save('data.mat','Weight','b','valaccbatch','valloss','trainaccbatch','trainloss','InitWeight','Initb');
