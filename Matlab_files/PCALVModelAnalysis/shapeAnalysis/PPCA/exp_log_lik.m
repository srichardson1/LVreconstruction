function val=exp_log_lik(observations,latents,loadings,means,sigmas,props,posterior_prob)
likelihoods=cell(numel(latents),1); %We have one likelihood array for each mixture component
for i=1:numel(latents)
 likelihood=zeros(numel(observations),1);
 for j=1:numel(likelihood)
     likelihood(j)=uni_obs_likelihood(observations(j),latents{i},loadings{i},means(i),sigmas(i));
 end
 likelihoods{i}=likelihood;
end
val=0;
for i=1:numel(observations)
    for j=1:numel(latents)
        val=val+posterior_prob(i,j)*ln(props(j)*uni_obs_likelihood(observations(i),latents{j},loadings{j},means(j),sigmas(j)));
    end
end
val;
end