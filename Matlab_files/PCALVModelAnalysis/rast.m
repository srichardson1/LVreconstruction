function y = rast(x)
% the default value of n = 2.
n = 1;
s = -1;
for j = 1:n
%s = s+(x(j)^2-10*cos(2*pi*x(j)));
s=1+x(j)^2;
end
y = 10*n+s;
y=s
