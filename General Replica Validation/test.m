
%func = @(x) (x-2)^2;
%soln=gridMin(func,-5,5,.001)


% 
%  B=[3,4,5;2,6,8;2,1,3]
%  
%  [val,col]=min(min(B))
%  [~,r]=min(B);
%  row=r(col)


func2d = @(x,y) (x-1).^2 +(y-2).^2;

[arg,err]=gridMin2d(func2d,0,0,3,3,.001)



x<g