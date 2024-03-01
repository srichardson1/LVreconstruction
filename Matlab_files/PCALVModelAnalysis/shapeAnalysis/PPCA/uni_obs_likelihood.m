function val=uni_obs_likelihood(obs,lat,load,mean,sig) %Univariate observation log-likelihood p(y|i) where i is mixture component
C=sig*diag(ones(numel(obs),1))+load*load';
C_inv=(eye(numel(obs))-eye(numel(obs))*load*inv(sig^2*eye(size(load,2))+load'*eye(numel(obs))*load)*load'*eye(numel(obs)))/sig^2;
E=(obs-mean)*inv(C)*(obs-mean)';
val=(2*pi)^(-numel(obs)/2)*det(C)^(-1/2);
end 