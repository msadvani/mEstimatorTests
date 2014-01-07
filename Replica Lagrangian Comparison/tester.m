%Testing out different things


%sigma = .000000001
%gauss = @(x)((2*pi.*sigma^2)^(-.5)*exp(-x.^2./(2.*sigma^2)));



clear all;
func = @(x) (2*x-1);

[output]=findZeroBB(func,.8,2,.001)


