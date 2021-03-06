function [x,y,val] = gridMin2d(fun,ax,ay,bx,by,tol)
%Attempts to find the smalles value of a function between a and b using a
%grid search. Reports the x and y location as well as the value of the
%function at this point.

N=10;



funSet = zeros(N);
while (abs(ax-bx)>tol || abs(ay-by)>tol)
    [ax,bx]
    [ay,by]
    xSet = repmat(linspace(ax,bx,N),N,1);
    ySet = repmat(linspace(ay,by,N)',1,N);
    
    for cnt = 1:numel(xSet);
        %[cnt,numel(xSet)]
        x = xSet(cnt);
        y=ySet(cnt);
        funSet(cnt) = fun(x,y);
    end

    
    [val,col]=min(min(funSet));
    [~,r]=min(funSet);
    row=r(col);
    
    
    if (col==1)
        ax = xSet(1,1);
        bx = xSet(1,2);
    elseif (col ==N)
        ax = xSet(1,N-1);
        ax = xSet(1,N);
    else
        ax = xSet(1,col-1);
        bx = xSet(1,col+1);
    end
    
    if (row==1)
        ay = ySet(1,1);
        by = ySet(2,1);
    elseif(row ==N)
        ay = ySet(1,N-1);
        by = ySet(1,N);
    else
        ay = ySet(row-1,1);
        by = ySet(row+1,1);
    end
    
        
    
%     
%     ay = ySet(row-1,1);
%     by = ySet(row+1,1);
end



x=mean([ax,bx]);
y = mean([ay,by]);

if(val>tol)
   error('I thought 0 error was expected minima') 
end


end

