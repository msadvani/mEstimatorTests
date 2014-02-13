load qOptaOptdata.mat

hold on;
h1 = plot(kappaSet, qOpt,'-')
h2 = plot(kappaSet, mean(aOpt),'r-')


set(h1,'LineWidth= 2.5')