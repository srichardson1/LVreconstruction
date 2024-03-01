function val=post_prob(observations,latents,loadings,means,sigmas,props)
val=zeros(numel(observations),numel(latents));
for k=1:numel(observations)
    for l=1:numel(latents)
        val(k,l)=posterior_prob(observations(k,:),latents{l},loadings{l},means(l,:),sigmas(l),props(l));
    end
end
val;
end