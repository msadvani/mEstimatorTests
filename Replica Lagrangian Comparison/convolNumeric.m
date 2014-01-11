function [ zSet] = convolNumeric(xSet,f1,f2)
%convolNumeric. Numerically convolves two functions and evaluates at x
%Note this code is pretty slow, that's ok...for now it does not need to be
%fast, if you want to speed it up you can use the convol function and
%compute the appropriate mapping

    zSet = zeros(size(xSet));

    zeta = @(x) sumIntIndef(@(y)(f1(x-y).*f2(y)),5, 20, .0001);

    for cnt =1:numel(xSet)
        x = xSet(cnt);
        zSet(cnt) = zeta(x);
    end
end

