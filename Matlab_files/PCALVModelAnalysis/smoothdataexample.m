close all;
x = 1:100;
A = cos(2*pi*0.05*x+2*pi*rand) + 0.5*randn(1,100);
%B = smooth(A(1:100),'lowess');
%B = smooth(A(1:100),'lowess',1);
B = smooth(A(1:100),'sgolay',0);
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
%B = filter(b,a,A(1:100));
plot(x,A,'-or',x,B,'-x')
legend('Original Data','Smoothed Data')