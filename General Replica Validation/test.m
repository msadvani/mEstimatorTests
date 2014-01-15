function [ out ] = test(x)
%to test breaking out of a function
if(x>0)
   out = x;
   return
end

out =0;

disp('hello')

end

