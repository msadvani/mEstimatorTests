hold on;


f = @(x) exp(-abs(x))/2;

gauss = @(x,s) (2*pi*s^2)^(-1/2)*exp(-x.^2/(2*s^2));


L=15;

 numX=1000;
 xSet  = linspace(-5,5,numX);
 plot(xSet,-log(f(xSet)),'k');


numX = 200;
xSet  = linspace(-5,5,numX);
q=2.5
%qSet = [.5,1,2,4];

    
    
    toMin = @(x,y) (log(convolve(y,f,@(w)gauss(w,q),L)) + (x-y).^2/(2*q));

    argmin = @(x)fminunc(@(y)toMin(x,y),randn);

    rho = @(x)-toMin(x,argmin(x));


   
    rhoSet = zeros(1,numX);


    for cnt = 1:numX
       x = xSet(cnt);   
       rhoSet(cnt) = rho(x);
    end


%    plot(xSet, rhoSet,'c')
    plot(xSet, rhoSet,'c')