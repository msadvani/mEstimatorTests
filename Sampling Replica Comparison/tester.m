
numKappa = 10;
kappaSet = linspace(.1,2,numKappa)


numSamples = 5;
kappaMat = repmat(kappaSet, numSamples, 1);

kappaSetVec = reshape(kappaMat,1,numKappa*numSamples)


hold on;

scatter(kappaSetVec,kappaSetVec,'r.')

plotProp = plot(kappaSet, kappaSet,'k');

set(plotProp,'Linewidth',3)