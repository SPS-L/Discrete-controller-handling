%% This code is just for practice code of ANN
%%% An example with 3 inputs and 1 output
%%% y=2a+3b+5c

% We need input and outputs for training and then we will be able to obtain
% output for any inputs without knowing the equation

%Assumbtion the inputs range from 0 to 1 [0-1]
% clear all
% clc
%Input data
i1=sqrt(xhist(11,:).^2+xhist(12,:).^2); %Input for 3Bus System
%i1=sqrt(xhist(14,:).^2+xhist(15,:).^2); %Input for SMIB System
% avrovhisti=[avrov,avrovhist];
% avrovhisti(:,end)=[];
% i2=avrovhisti;
%Output data
y1=hhist(5,:);
% y2=avrovhist;

% CI=[i1;i2];
% CO=[y1;y2];
CI=i1;
CO=y1;

% net=newff(minmax(CI),[9,10,9],{'logsig','tansig','purelin'},'trainlm'); % 3 5 1 InputLayer HiddenLayer OutputLayer
net=newff(minmax(CI),[1,3,1],{'logsig','tansig','purelin'},'trainlm'); % 3 5 1 InputLayer HiddenLayer OutputLayer
% net=newff(minmax(CI),[1,3,1],{'logsig','elliot2sig','purelin'},'trainlm'); % 3 5 1 InputLayer HiddenLayer OutputLayer
% net=newff(minmax(CI),[1,3,1],{'logsig','elliot2sig','purelin'},'trainlm'); % 3 5 1 InputLayer HiddenLayer OutputLayer
% net=feedforwardnet(17,'trainlm');
net=init(net);                          %Initialize the Network (weight and biases)
net.trainParam.show=1;                  %The result of error at each iteration
net.trainParam.epochs=10000;              % Maximum limit of the network training itiration
net.trainParam.goal=1e-12;              % Stopping criterion based on error goal
net=train(net,CI,CO);                   % Train the network
tic
y=net(1)
toc

