lb = [-10];
ub = [20];
%[x,fval,exitflag] = ga(@lincontest6,... 
%    2,A,b,[],[],lb)
[x fval,exitflag] = ga(@rast,1,[],[],[],[],lb,ub)




