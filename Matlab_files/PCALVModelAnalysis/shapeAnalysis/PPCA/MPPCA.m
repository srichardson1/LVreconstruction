%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% MPPCA EM Algorithm %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% function val=Gprior(x,mu,sigma) %GAUSSIAN PRIOR
% val=mvnpdf(x,mu,sigma);
% end
% 
% 
% function val=uni_obs_likelihood(obs,lat,load,mean,sig) . %Univariate observation log-likelihood p(y|i) where i is mixture component
% C=sig*diag(ones(numel(obs),1))+load*load'
% C_inv=(eye(numel(obs))-load*inv(sig^2*eye(numel(obs))+load'*load)*load')/sig^2
% E=(obs-mean)'*inv(C)*(obs-mean)
% val=(2*pi)^(-numel(obs)/2)*det(C)^(-1/2);
% end 
% 
% function val=obs_likelihood(obs,lat,load,mean,sig) %Log-likelihood for set of observations where obs contains observations in rows
% val=(2*pi*sig^2)^-(size(obs,1)/2)*exp(-1/(2*sig^2)*norm(obs-(lat*load+mean)));
% end
% 
% function val=marginal_likelihood(obs,lat,load,mean,sig,prop) %likelihood marginalised over mixtures.  
% likelihoods=zeros(numel(obs),1)
% for i=1:numel(mean)
%     likelihoods(i)=uni_obs_likelihood(obs,lat{i},load{i},mean(i),sig(i));
% end
% val=sum(likelihoods*prop);
% end
%     
% function val=posterior_prob(obs,lat,load,mean,sig,prop);
% lik=uni_obs_likelihood(obs,lat,load,mean,sig);
% marg=marginal_likelihood(obs,lat,load,mean,sig,prop);
% val=li*prop/marg;
% end
% 
% function val=post_prob(observations,latents,loadings,means,sigmas,props)
% val=zeros(numel(observations),numel(latents));
% for k=1:numel(observations)
%     for l=1:numel(latents)
%         val(k,l)=posterior_prob(observations(k),latents{l},loadings{l},means(l),sigmas(l),props(l));
%     end
% end
% val;
% end
% 
% function val=exp_log_lik(observations,latents,loadings,means,sigmas,props,posterior_prob)
% likelihoods=cell(numel(latents),1); %We have one likelihood array for each mixture component
% for i=1:numel(latents)
%  likelihood=zeros(numel(observations),1);
%  for j=1:numel(likelihood)
%      likelihood(j)=uni_obs_likelihood(observations(j),latents{i},loadings{i},means(i),sigmas(i));
%  end
%  likelihoods{i}=likelihood;
% end
% val=0;
% for i=1:numel(observations)
%     for j=1:numel(latents)
%         val=val+posterior_prob(i,j)*ln(props(j)*uni_obs_likelihood(observations(i),latents{j},loadings{j},means(j),sigmas(j)));
%     end
% end
% val;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Read in Data  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

workingDir = pwd();

EndoSurfaceNodes = [];
SubjecIndex = 0;

subjectindex=[];
counter=1;
for i = 1:90;
    if i<10;
    if exist(strcat('/Users/alan/ownCloud/Results/HV0',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
    end;
    else;
        if exist(strcat('/Users/alan/ownCloud/Results/HV',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
        end;
    end;
end;

for i = 1:(numel(subjectindex)-1); %Extracting Endo surface nodes
    if subjectindex(i)<10
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('/Users/alan/ownCloud/Results/HV0',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,1:numel(sendo)) = sendo;
    else
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,1:numel(sendo)) = sendo;
    end
end

[sendo, sepi, abaqusInputData] = load_LVGeoData_1(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(numel(subjectindex))),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,:) = [sendo];

for i=1:24
    if exist(strcat('/Users/alan/ownCloud/MI/Patient',num2str(i),'_prolatePara.mat'))==2;
    [sendo, sepi, ~]=load_LVGeoData_1(strcat('/Users/alan/ownCloud/MI/Patient',num2str(i),'/earlyDiastole'), ...
                     workingDir);
             EndoSurfaceNodes=[EndoSurfaceNodes; sendo];
    end
end

    load('/Users/alan/ownCloud/PCA_Data/case4-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint1=reshape(sendo',1,2896*3);

    load('/Users/alan/ownCloud/PCA_Data/case4-followup-9months/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint2=reshape(sendo',1,2896*3);
AMYpoints=[AMYpoint1; AMYpoint2];

points=[EndoSurfaceNodes; AMYpoints];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% MPPCA with K means Clustering %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = kmeans(points,2);

set1=points(idx==1,:);
set2=points(idx==2,:);

propns=[numel(set1)/numel(points) numel(set2)/numel(points)];

mean1=mean(points(idx==1,:));
mean2=mean(points(idx==2,:));
means=[mean1; mean2];
sigmas=rand(2,1);

for row = 1 : size(set1,1)
    for col = 1 : size(set1,2)
        set1_standardize(row, :) = set1(row, :) -  mean1;
    end
end

for row = 1 : size(set2,1)
    for col = 1 : size(set2,2)
        set2_standardize(row, :) = set2(row, :) -  mean2;
    end
end

C1=set1_standardize'*set1_standardize;
C2=set2_standardize'*set2_standardize;
[u, d, v] = svd(C1, 0);
Loadings1=u(:,1:5);
[u, d, v] = svd(C2, 0);
Loadings2=u(:,1:5);


Latent_var1=normrnd(0,1,5,1);
Latent_var2=normrnd(0,1,5,1);
Latents=cell(2,1);
Latents{1}=Latent_var1;
Latents{2}=Latent_var2
%Loadings1=normrnd(0,1,8688,5);
%Loadings2=normrnd(0,1,8688,5);
Loadings=cell(2,1);
Loadings{1}=Loadings1;
Loadings{2}=Loadings2;
%save('/Users/alan/ownCloud/PhD/Papers/BHF/PCALVModelAnalysis/shapeAnalysis/PPCA/svd.mat',Loadings);
%%
load('/Users/alan/ownCloud/PhD/Papers/BHF/PCALVModelAnalysis/shapeAnalysis/PPCA/svd.mat');
idx = kmeans(points,2);

set1=points(idx==1,:);
set2=points(idx==2,:);
sigmas=[100 100];
propns=[numel(set1)/numel(points) numel(set2)/numel(points)];

mean1=mean(points(idx==1,:));
mean2=mean(points(idx==2,:));
means=[mean1; mean2];
Latent_var1=normrnd(0,1,3,1);
Latent_var2=normrnd(0,1,3,1);
Latents=cell(2,1);
Latents{1}=Latent_var1;
Latents{2}=Latent_var2
Loadings{1}=Loadings{1}(:,1:3);
Loadings{2}=Loadings{2}(:,1:3);
diff=1
for iter=1:100
    n_obs=size(points,1);
    dim_obs=size(points,2);
    dim_lat=size(Loadings{1},2);
    M_invs=cell(2,1);
M_invs{1}=inv(sigmas(1)*eye(dim_lat)+Loadings{1}'*Loadings{1});
M_invs{2}=inv(sigmas(2)*eye(dim_lat)+Loadings{2}'*Loadings{2});
inv_covs=cell(2,1); %Arrays of inverse covariance matrices
inv_covs{1}=(eye(8688)-Loadings{1}*M_invs{1}*Loadings{1}')/sigmas(1)^2
inv_covs{2}=(eye(8688)-Loadings{2}*M_invs{2}*Loadings{2}')/sigmas(2)^2;
mixture_likelihoods=cell(2,1); %Arrays of likelihood of each observation in the different mixtures
for mixture=1:2 %change 2 for higher number of mixtures
    det1=det(sigmas(mixture)*eye(dim_obs)+Loadings{mixture}*Loadings{mixture}');
    mixture_likelihood=zeros(size(points,1),1); 
    for obsind=1:n_obs
        %mixture_likelihood(obsind)=(2*pi)^(-size(points,2)/2)*1/sqrt(det1)*exp((points(obsind,:)-means(mixture,:))*inv_covs{mixture}*(points(obsind,:)-means(mixture,:))');
        mixture_likelihood(obsind)=(points(obsind,:)-means(mixture,:))*inv_covs{mixture}*(points(obsind,:)-means(mixture,:))';
    end
    mixture_likelihoods{mixture}=mixture_likelihood;
end

%marginal_likelihoods=zeros(size(points,1),1);
%for i=1:size(points,1)
    marginal_likelihoods=propns(1)*mixture_likelihoods{1}+propns(2)*mixture_likelihoods{2};
%end

posterior_probs=zeros(n_obs,numel(sigmas)); %Calculating matrix of R_ni values
 for k=1:n_obs
     for l=1:numel(sigmas)
         posterior_probs(k,l)=mixture_likelihoods{l}(k)*propns(l)/marginal_likelihoods(k);
     end
 end
 
propns=sum(posterior_probs)/n_obs; %Updated mixture proportions
mean1=sum(bsxfun(@times,points,posterior_probs(:,1)))/sum(posterior_probs(:,1));
mean2=sum(bsxfun(@times,points,posterior_probs(:,2)))/sum(posterior_probs(:,2));
means=[mean1; mean2]; %Updated mixture means

%Posterior sufficient statistics

x=cell(numel(sigmas),1);
for i=1:numel(sigmas)
    x_ind=zeros(dim_lat,n_obs);
    for j=1:n_obs
        x_ind(:,j)=M_invs{i}*Loadings{i}'*(points(j,:)-means(i,:))';
    end
    x{i}=x_ind;
end

x_2=cell(numel(sigmas),n_obs);
for i=1:numel(sigmas)
    for j=1:n_obs
        x_2{i,j}=sigmas(i)^2*M_invs{i}+x{i}(:,j)*x{i}(:,j)';
    end
end
        

Sample_covariance=cell(2,1);
obs_standardiseds=cell(2,1); %Standardised observations for calculating sample covariance
for mixture=1:numel(sigmas)
    obs_standardised=zeros(size(points,1),size(points,2));
    for j=1:n_obs
    obs_standardised(j,:)=points(j,:)-means(mixture,:);
    end
    obs_standardiseds{mixture}=obs_standardised;
end
Sample_covariance{1}=zeros(size(points,2),size(points,2));
Sample_covariance{2}=zeros(size(points,2),size(points,2));
for j=1:numel(sigmas) %Sample covariance for each of the mixtures
for i=1:n_obs
    Sample_covariance{j}=Sample_covariance{j}+1/(propns(j)*n_obs)*posterior_probs(i,j)*obs_standardiseds{j}'*obs_standardiseds{j};
end
end
Loadings_update=cell(2,1);
for mixture=1:numel(sigmas)
    Loadings_update{mixture}=Sample_covariance{mixture}*Loadings{mixture}*inv(sigmas(mixture)^2*eye(size(M_invs{mixture},1))+M_invs{mixture}*Loadings{mixture}'*Sample_covariance{mixture}*Loadings{mixture});
end

for mixture=1:numel(sigmas)
sigmas(mixture)=1/size(points,2)*trace(Sample_covariance{mixture}-Sample_covariance{mixture}*Loadings{mixture}*M_invs{mixture}*Loadings_update{mixture}');
end
Loadings=Loadings_update;
end


%Using the Loadings to reconstruct LVs and project to principal component
%space

reconstruct=Loadings{1}'*points';
