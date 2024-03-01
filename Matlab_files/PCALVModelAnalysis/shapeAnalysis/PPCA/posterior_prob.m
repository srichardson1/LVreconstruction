function val=posterior_prob(obs,lat,load,mean,sig,prop)
lik=uni_obs_likelihood(obs,lat,load,mean,sig);
marg=marginal_likelihood(obs,lat,load,mean,sig,prop);
val=li*prop/marg;
end