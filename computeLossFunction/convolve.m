function [ out ] = convolve(xSet,fun1,fun2,L)
%Convolve 
out = zeros(size(xSet));




% for cnt = 1:numel(xSet)
%     xVal = xSet(cnt);
%     toConv  = @(x,y)fun1(x-y).*fun2(y);
%     out(cnt) = quad(@(y)toConv(xVal,y),-L,L);
% end

N=10^4;
for cnt = 1:numel(xSet)
    xVal = xSet(cnt);
    toConv  = @(x,y)fun1(x-y).*fun2(y);
    ySet = linspace(-L,L,N);
    out(cnt) = sum(toConv(xVal,ySet))*(2*L/N);
end





end

