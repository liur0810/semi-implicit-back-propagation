clearvars;
load mnist_train.mat; %cifar10_train.mat
X = gpuArray(X);
Y = gpuArray(Y);
[Xsize,Xnumber] = size(X);
[Ysize,~] = size(Y);
layersize = [Xsize,500,Ysize]; 
N = length(layersize);
for i=2:N
    Weight{i-1} =  gpuArray(normrnd(0,0.05,[layersize(i) layersize(i-1)]));
    b{i-1} = gpuArray(normrnd(0,0.05,[layersize(i) 1]));
    lambda{i-1} = [1,1,1,1];
end
InitWeight = Weight;
Initb = b;
batchsize = 100;
maxiter = 5;
ita = 0.1;
[Weight,b,valaccbatch,valloss,trainaccbatch,trainloss,batchtime] = sibp(X,Y,maxiter,batchsize,InitWeight,Initb,lambda,layersize,ita);
save('data','Weight','b','valaccbatch','valloss','trainaccbatch','trainloss','InitWeight','Initb','batchtime');
