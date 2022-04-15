function Y = neuralNetwork(trainVec, trainY, testVec)
trainData = table(trainVec,trainY);
layers = [
    featureInputLayer(1,'Name','position')
    fullyConnectedLayer(32,'Name','fc1')
    reluLayer('Name','rl1')
    fullyConnectedLayer(32,'Name','fc2')
    reluLayer('Name','rl2')
    fullyConnectedLayer(16,'Name','fc3')
    reluLayer('Name','rl3')
    fullyConnectedLayer(16,'Name','fc4')
    reluLayer('Name','rl4')
    fullyConnectedLayer(8,'Name','fc5')
    reluLayer('Name','rl5')
    fullyConnectedLayer(8,'Name','fc6')
    reluLayer('Name','rl6')
    fullyConnectedLayer(4,'Name','fc7')
    reluLayer('Name','rl7')
    fullyConnectedLayer(4,'Name','fc8')
    reluLayer('Name','rl8')
    fullyConnectedLayer(1,'Name','fc9')
    regressionLayer('Name','output')];
lgraph = layerGraph(layers);
% Define trainingOptions and also set 'Shuffle' to 'never' for this workaround to work
options = trainingOptions('adam', ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise',...
    'MaxEpochs',300, ...
    'MiniBatchSize',1024, ...
    'Verbose',1, ...
    'Plots','training-progress',...
    'Shuffle','never');
net = trainNetwork(trainData,lgraph,options);
Y = predict(net,testVec);
end