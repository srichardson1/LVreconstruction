clear all; close all; clc;
workingDir = pwd();
%%%Saved data load
%load('/Users/alan/Desktop/PhD/Papers/BHF/PCA_reconstructions/Prolate/SVD.mat')

%%%generate data

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
EndoSurfaceNodes = zeros(numel(subjectindex),2896*3);
%load the data
%[sendo, sepi, ~] = load_LVGeoData('/Users/alan/Desktop/PhD/Papers/BHF/Results/HV02/earlyDiastole', ...
%                      workingDir);
%SubjecIndex = SubjecIndex + 1;
%EndoSurfaceNodes(SubjecIndex,:) = sendo;


%[sendo, sepi, ~] = load_LVGeoData_3('/Users/alan/Desktop/PhD/Papers/BHF/Results/HV03/earlyDiastole', ...
%                      workingDir);
%SubjecIndex = SubjecIndex + 1;
%EndoSurfaceNodes(SubjecIndex,:) = sendo;

%[sendo, sepi, abaqusInputData] = load_LVGeoData('/Users/alan/Desktop/PhD/Papers/BHF/Results/HV04/earlyDiastole', ...
%                      workingDir);
%SubjecIndex = SubjecIndex + 1;
%EndoSurfaceNodes(SubjecIndex,:) = sendo;

%[sendo, sepi, abaqusInputData] = load_LVGeoData('/Users/alan/Desktop/PhD/Papers/BHF/Results/HV05/earlyDiastole', ...
 %                     workingDir);
%SubjecIndex = SubjecIndex + 1;
%EndoSurfaceNodes(SubjecIndex,:) = sendo;

for i = 1:(numel(subjectindex)-1); %Extracting Endo surface nodes
    if subjectindex(i)<10
 [sendo, sepi, ~] = load_LVGeoData_3(strcat('/Users/alan/ownCloud/Results/HV0',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,:) = sendo';
    else
 [sendo, sepi, ~] = load_LVGeoData_3(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(i)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,:) = sendo';
    end
end
        
[sendo, sepi, abaqusInputData] = load_LVGeoData_3(strcat('/Users/alan/ownCloud/Results/HV',num2str(subjectindex(end)),'/earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,:) = sendo';

[sendo, sepi, abaqusInputData] = load_LVGeoData_3(strcat('/Users/alan/ownCloud/Results/HV02/earlyDiastole'), ...
                     workingDir);

                 EndoSurfaceNodes(:,2:3:size(EndoSurfaceNodes,2))=EndoSurfaceNodes(:,2:3:size(EndoSurfaceNodes,2))*pi/180;
                 EndoSurfaceNodes(:,1:3:size(EndoSurfaceNodes,2))=EndoSurfaceNodes(:,2:3:size(EndoSurfaceNodes,2))*pi/180;
                 abaqusInputData.node(:,4:5)=abaqusInputData.node(:,4:5)*pi/180;
      
%%calculate average shape
EndoSurfaceNodes_ave = mean(EndoSurfaceNodes);
%%minus the average shape from EndoSurfaceNodes

for row = 1 : size(EndoSurfaceNodes,1)
        EndoSurfaceNodes_standardize(row, :) = EndoSurfaceNodes(row, :) -  EndoSurfaceNodes_ave;
end

%%calculate the covariance matrix 
N = size(EndoSurfaceNodes,1);
C = EndoSurfaceNodes_standardize'*EndoSurfaceNodes_standardize;

%%that is coming from http://mghassem.mit.edu/pcasvd/
[u, d, v] = svd(C, 0); %
% Pull out eigen values and vectors
eigVals = diag(d);
eigVecs = u;
eigVals
%save('/Users/alan/Desktop/PhD/Papers/BHF/PCA_reconstructions/Prolate/SVD.mat')
%%
%%the cumulative energy content for the m'th eigenvector is the sum of the
%%energy content across eigenvalues 1: m

for i = 1 : N
    energy(i) = sum(eigVals(1:i));
end
propEnergy = energy./energy(end);

%%determine the number of principal components required to model 90% of data
%%variance

percentMark = min(find(propEnergy > 0.95));

%percentMark = 1;
eigVecsSelected = u(:,1:percentMark);

%%we can calculate all the weights 
Nweights = (EndoSurfaceNodes_standardize)*eigVecsSelected;

%%next one we will need to output the mesh using tecplot
load('/Users/alan/ownCloud/Results/alpha0_global.mat')
cd('/Users/alan/Desktop/PhD/Papers/BHF/PCA_reconstructions/Prolate')
bar(propEnergy(1:20))
title('Percentage Variation Explained')
ylabel('Variation Explained')
xlabel('Number of Principal Components')
print('/Users/alan/Desktop/PhD/Reports/PCA_BHF_2018_01_17/Prolate_PCA/Prolate_PC_variation_bar','-dpng')
alpha0 = alpha0_global;
%%
tecFid = fopen('mean_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,EndoSurfaceNodes_ave, alpha0,'endo', tecFid,2);
fclose(tecFid);

tecFid = fopen('ori_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,[], alpha0,'endo', tecFid,2);
fclose(tecFid);
%%
%%let us reconstruct the endo face with percentMark eigen vectors
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
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave'), alpha0,'endo', tecFid,2);
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
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave'), alpha0,'endo', tecFid,2);
fclose(tecFid);
end

endoFaceSample = EndoSurfaceNodes(1,:);
endoFaceSample = endoFaceSample - EndoSurfaceNodes_ave;
faceWeight = eigVecsSelected'*(endoFaceSample');
mapped_endo = zeros(size(eigVecsSelected(:,1)));
for i = 1 : 5
  %    mapped_endo = mapped_endo-2*std(Nweights(:,1))*eigVecsSelected(:,i);
     mapped_endo = mapped_endo+faceWeight(i).*u(:,i);
end

tecFid = fopen('reconstructed_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo+EndoSurfaceNodes_ave', alpha0,'endo', tecFid,2);
fclose(tecFid);

eigVecsSelected=u(:,1:size(u,1));
endoFaceSample = EndoSurfaceNodes(1,:);
endoFaceSample = endoFaceSample - EndoSurfaceNodes_ave;
faceWeight = eigVecsSelected'*(endoFaceSample');
mapped_endo = zeros(size(eigVecsSelected(:,1)));
mapped_endo = endoFaceSample*eigVecsSelected*eigVecsSelected';

tecFid = fopen('reconstructed_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,(mapped_endo+EndoSurfaceNodes_ave)', alpha0,'endo', tecFid,2);
fclose(tecFid);



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