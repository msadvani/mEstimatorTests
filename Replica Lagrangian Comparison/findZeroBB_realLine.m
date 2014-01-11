function [out] = findZeroBB_realLine(func,a,b,tol)
%FindZeroBB - uses a branch and bound method to find zeros of a monotonic
%increasing function passing zero.
    maxIter = 10^6;  %This is a global variable
    
    if (a<=0)
        error('a must be greater than 0')
    end
    
    if (b<a)
        error('b must be greater than a')
    end


%% find a large lower bound a    
    cnt=1;
    while(func(a)>0)
        a=a-(b-a);
        cnt=cnt+1;
        if(cnt>maxIter)
           error('maxIteration exceeded in finding largeNegPt - the function may not be negative between 0 and a') 
        end
    end

%% find a small upper bound on b
    cnt = 1;
    b0=b;
    while(func(b)<0)
        b=b+(b0-a);
        cnt=cnt+1;
        if(cnt>maxIter)
           error('maxIteration exceeded in finding smallPosPt - the function may not be positive betweennegative between 0 and b') 
        end
    end
    
    
    if(a>b ||func(a)>0 || func(b)<0)
       error('The function may not be monotonic increasing'); 
    end
    
    
    
    
    d = b-a;
     
%     disp('a')
%     disp(a)
%     
%     disp('b')
%     disp(b)
    
    while(d>tol)
        
        m = (b+a)/2;
        func(m);
        if (func(m)>0)
           b=m; 
        else
           a=m;
        end
        
        d=b-a;
    end
    
    out = m;
end

