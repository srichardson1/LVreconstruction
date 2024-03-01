clear all; close all; clc;
workingDir = pwd();

SurfaceNodes = [];
SubjecIndex = 0;

subjectindex=[];
counter=1;
for i = 1:90
    if i<10
    if exist(strcat('/Users/alan/ownCloud/Results/HV0',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
    end
    else
        if exist(strcat('/Users/alan/ownCloud/Results/HV',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
        end
    end
end

for i = 1:(numel(subjectindex)-1); %Extracting Endo surface nodes
    if subjectindex(i)<10
 [sendo, sepi, ~] = load_LVGeoData(strcat('/Users/alan/ownCloud/Results/HV0',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
SurfaceNodes(SubjecIndex,1:(2*numel(sepi))) = [sendo' sepi'];
    else
 [sendo, sepi, ~] = load_LVGeoData(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
SurfaceNodes(SubjecIndex,1:(2*numel(sepi))) = [sendo' sepi'];
    end
end
       
[sendo, sepi, abaqusInputData] = load_LVGeoData(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(numel(subjectindex))),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
SurfaceNodes(SubjecIndex,1:(2*numel(sepi))) = [sendo' sepi'];
[sendo, sepi, abaqusInputData] = load_LVGeoData(strcat('/Users/alan/ownCloud/Results/HV02/earlyDiastole'), ...
                     workingDir);
%%calculate average shape
SurfaceNodes_ave = mean(SurfaceNodes);

%%
%%minus the average shape from EndoSurfaceNodes
for row = 1 : size(SurfaceNodes,1)
    for col = 1 : size(SurfaceNodes,2)
        SurfaceNodes_standardize(row, col) = SurfaceNodes(row, col) -  SurfaceNodes_ave(col);
    end
end

%%calculate the covariance matrix 
N = size(SurfaceNodes,1);
C = SurfaceNodes_standardize'*SurfaceNodes_standardize;

%%that is coming from http://mghassem.mit.edu/pcasvd/
[u, d] = svd(C); %
% Pull out eigen values and vectors
eigVals = diag(d);
eigVecs = u;

%%the cumulative energy content for the m'th eigenvector is the sum of the
%%energy content across eigenvalues 1: m
for i = 1 : N
    energy(i) = sum(eigVals(1:i));
end
propEnergy = energy./energy(end);

%%determin the number of principal components required to model 90% of data
%%variance
percentMark = min(find(propEnergy > 0.9));
%percentMark = 1;
eigVecsSelected = u(:,1:percentMark);


%%we can calculate all the weights      

Nweights = (SurfaceNodes_standardize)*eigVecsSelected;
cd('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Radial')
bar(propEnergy(1:20))
title('Percentage Variation Explained')
ylabel('Variation Explained')
xlabel('Number of Principal Components')
%print('/Users/alan/ownCloud/PhD/Reports/PCA_BHF_2018_01_17/Radial_PC_variation_bar','-dpng')
%%next one we will need to output the mesh using tecplot
load('/Users/alan/ownCloud/Results/alpha0_global.mat')
alpha0 = alpha0_global;
tecFid = fopen('mean_mesh_endo.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,SurfaceNodes_ave, alpha0,'endo', tecFid);
fclose(tecFid);
tecFid = fopen('mean_mesh_epi.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,SurfaceNodes_ave(2897:5792), alpha0,'epi', tecFid);
fclose(tecFid);

tecFid = fopen('ori_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,[], alpha0,'endo', tecFid);
fclose(tecFid);

%%

%%let us reconstruct the endo face *of one subject* with percentMark eigen vectors
for mode=1:percentMark
eigVecsSelected = u(:,1:percentMark);
endoFaceSample = EndoSurfaceNodes(1,:);
endoFaceSample = endoFaceSample - EndoSurfaceNodes_ave;
faceWeight = eigVecsSelected'*(endoFaceSample');
identity=zeros(percentMark,1);
identity(mode)=1
    mapped_endo = -2*eigVecsSelected*std(Nweights(:,mode))*identity;  
    % mapped_endo = mapped_endo+faceWeight(i).*eigVecsSelected(:,i);
tecFid = fopen(strcat('reconstructed_mesh_',num2str(mode),'_minus2.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave'), alpha0,'endo', tecFid,3);
fclose(tecFid);

eigVecsSelected = u(:,1:percentMark);
endoFaceSample = EndoSurfaceNodes(1,:);
endoFaceSample = endoFaceSample - EndoSurfaceNodes_ave;
faceWeight = eigVecsSelected'*(endoFaceSample');
identity=zeros(percentMark,1);
identity(mode)=1;
    mapped_endo = 2*eigVecsSelected*std(Nweights(:,mode))*identity;  
    % mapped_endo = mapped_endo+faceWeight(i).*eigVecsSelected(:,i);
tecFid = fopen(strcat('reconstructed_mesh_',num2str(mode),'_plus2.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave'), alpha0,'endo', tecFid,3);
fclose(tecFid);
end
for i =percentMark %[1 5 10 15 20 25 30 35 40 45 50]
    eigVecsSelected=u(:,1:i);
endoFaceSample = EndoSurfaceNodes(1,:);
endoFaceSample = endoFaceSample - EndoSurfaceNodes_ave;
faceWeight = eigVecsSelected'*(endoFaceSample');
mapped_endo = zeros(size(eigVecsSelected(:,1)));
mapped_endo = endoFaceSample*eigVecsSelected*eigVecsSelected';


tecFid = fopen(strcat('reconstructed_mesh.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave)', alpha0,'endo', tecFid,3);
fclose(tecFid);
end


%%

figure;
subplot(231);
plot(Nweights(:,1),'.');
subplot(232);
plot(Nweights(:,2),'.');
subplot(233);
plot(Nweights(:,3),'.');
subplot(234);
plot(Nweights(:,4),'.');


k = find(baseline_scan==1);
k2 = find(baseline_scan==0);


figure;
subplot(231);
plot(Nweights(k,1),'.'); hold on;plot(Nweights(k2,1),'r+'); hold on;
subplot(232);
plot(Nweights(k,2),'.'); hold on;plot(Nweights(k2,2),'r+'); hold on;
subplot(233);
plot(Nweights(k,3),'.'); hold on;plot(Nweights(k2,3),'r+'); hold on;
subplot(234);
plot(Nweights(k,4),'.'); hold on;plot(Nweights(k2,4),'r+'); hold on;





% % Do something with them; for example, project each of the neutral and smiling faces onto the corresponding eigenfaces
% neurtalFaces = faces(:, neutral); smileFaces = faces(:, smile);
% neutralWeights = eigenVecs' * neutralFaces;
% smileWeights = eigenVecs' * smileFaces;
% 
% % Use the coefficients of these projections to classify each smiling face
% for i = 1:numFaces
% weightDiff = repmat(smileWeights(:, i), 1, numFaces) - neutralWeights;
% [val, ind] = min(sum(abs(weightDiff), 1));
% bestMatch(i) = ind;
% end






