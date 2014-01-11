clear all;
close all;


%fNoise = @(x) exp(x)./((1+exp(x)).^2);

%fNoise = @(x) (1/2)*exp(-abs(x));


fNoise = @(x) (2*pi).^(-1/2)*exp(-x.^2./(2));


g = @(x) (1/2)*exp(-abs(x));

gauss = @(x,sigma) (2*pi*sigma.^2).^(-1/2)*exp(-x.^2./(2*sigma.^2));

xSet = linspace(-5,5,500);
%plot(xSet,gauss(xSet,2))

zeta = @(x,sigma) sumIntIndef(@(y)(fNoise(x-y).*gauss(y,sigma)), 5, 20, .0001);

zetaSet = zeros(size(xSet));
for cnt = 1:length(xSet)
   xVal = xSet(cnt);
   zetaSet(cnt) = zeta(xVal,2); 
end

plot(xSet, zetaSet)

z2 = convolNumeric(xSet,@(y)fNoise(y),@(x)gauss(x,.1));


hold on;
plot(xSet,z2,'r')
