function val=Gprior(x,mu,sigma) %GAUSSIAN PRIOR
val=mvnpdf(x,mu,sigma);
end