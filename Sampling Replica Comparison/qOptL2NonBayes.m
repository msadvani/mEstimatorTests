fNoise = @(x) (1/2)*exp(-abs(x));
kappa = .6;

%% Optimal error

L0 = 10; %initial interval length
N0=10^3; %initial num point in integral approx
dblErr=.001; %permissible error in I
IQ = @(x) IQfunc(fNoise,x,L0,N0, dblErr);

fToMin = @(z)(IQ(z) - kappa).^2;

qOpt = gridMinSearch(fToMin,0,5,100,.001);

qOpt

% 
% xSet = linspace(0,5,100);
% 
% [val,arg]=min(fToMin(xSet));

%based on Val...once you've achived a certain threshold, you
%stop...otherwise you iterate zooming in on the area...
%xMin = xSet(arg);
%qOpt = xMin;
% plot(xSet, fToMin(xSet))






%Now find where IQ intersects with kappa...
%funcOptQ = @(z) ((IQ(z)-kappa).^2);
%disp('qOpt approx');
%qOpt = fminunc(funcOptQ,rand)
%qOpt = fmincon(funcOptQ,rand,-1,0);

% qOpt
% % 
% % 
qL2 = @(k) 2*k./(1-k);
disp('qL2')
% qL2(kappa)
 disp('qOpt/qL2')
 
qOpt./qL2(kappa)