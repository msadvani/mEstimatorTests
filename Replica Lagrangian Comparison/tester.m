fNoise = @(x) exp(x(x>0))+x(x<0)

fNoise(-1)
fNoise(1)