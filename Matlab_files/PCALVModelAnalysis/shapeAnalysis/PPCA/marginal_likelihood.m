function val=marginal_likelihood(obs,lat,load,mean,sig,prop) %likelihood marginalised over mixtures.  
likelihoods=zeros(numel(obs),1)
for i=1:numel(mean)
    likelihoods(i)=uni_obs_likelihood(obs,lat{i},load{i},mean(i),sig(i));
end
val=sum(likelihoods*prop);
end