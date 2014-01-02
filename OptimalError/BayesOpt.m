clear all;
close all;

kappa = .5;


fNoise = @(x) (2*pi)^(-1/2)*exp(-x.^2/2);
g = @(x) (1/2)*exp(-abs(x));


%% Theory for Optimal error

L0 = 10; %initial interval length
N0=10^3; %initial num point in integral approx
dblErr = .0001; %permissible error in I


%Ooops...maybe need to make the vector veions more general...
I = @(x) calculateInfoMat(fNoise,x,L0,N0, dblErr);
J = @(x) calculateInfoMat(g,x, L0, N0, dblErr);



%plotting the function
%a = linspace(0,50,200);
%plot(a,a - a.^2.*J(a))

%two constraints
F1 = @(q,a) (a - a.^2.*J(a) - q).^2;
F2 = @(q,a) (I(q).*a - kappa).^2;


% numA = 50;
% numQ = 50;
% aSet =repmat(linspace(0,400,numA),numQ,1);
% qSet = repmat(linspace(0,10,numQ)',1,numA);


aMin1 = @(q)gridMinSearch(@(x)F1(q,x),0,5,10,.01);
aMin2 = @(q)gridMinSearch(@(x)F2(q,x),0,5,10,.01);


d = @(q)(aMin1(q)-aMin2(q)).^2;

qOpt = gridMinSearchNonVec(d,.05,1,10,.005)

% qSet = linspace(.1,1,8);
%  
% dSet = zeros(size(qSet));
% for cnt = 1:length(qSet)
%     dSet(cnt) = d(qSet(cnt));
% end
% 
% plot(qSet, dSet)






%qSet = linspace(.1,1,5);
% aMin1Set = zeros(size(qSet));
% aMin2Set = zeros(size(qSet));
% 
% for cnt = 1:length(qSet)
%     aMin1Set(cnt) = aMin1(qSet(cnt));
%     aMin2Set(cnt) = aMin2(qSet(cnt));
% end
% 
% hold on;
% 
% plot(qSet, aMin1Set)
% plot(qSet, aMin2Set,'r')




%qSet = linspace(.1,1,5);
% aMin1Set = zeros(size(qSet));
% aMin2Set = zeros(size(qSet));
% 
% for cnt = 1:length(qSet)
%     aMin1Set(cnt) = aMin1(qSet(cnt));
%     aMin2Set(cnt) = aMin2(qSet(cnt));
% end
% 
% hold on;
% 
% plot(qSet, aMin1Set)
% plot(qSet, aMin2Set,'r')



%plot()
% qSet = linspace(.1,1,10);
% 
% 
% 
% 
% 
% 
% aMin1Set = zeros(size(qSet));
% aMin2Set = zeros(size(qSet));
% 
% for cnt = 1:length(qSet)
%     aMin1Set(cnt) = aMin1(qSet(cnt));
%     aMin2Set(cnt) = aMin2(qSet(cnt));
% end
% 
% hold on;
% 
% plot(qSet, aMin1Set)
% plot(qSet, aMin2Set,'r')
% 
% legend('aMin1','aMin2');

%% Compute the integral estimating the error for lasso










%% Theory for a suboptimal case (L2 norm) [This comes from replica equations]

% f=@(lam)qThyL2(kappa,lam)
% 
% lamMin = gridMinSearch(f,0,10,100,.01);
% 
% qThyMin = f(lamMin)





