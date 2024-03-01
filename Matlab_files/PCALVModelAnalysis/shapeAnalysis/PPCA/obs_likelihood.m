function val=obs_likelihood(obs,lat,load,mean,sig) %Log-likelihood for set of observations where obs contains observations in rows
val=(2*pi*sig^2)^-(size(obs,1)/2)*exp(-1/(2*sig^2)*norm(obs-(lat*load+mean)));
end