clear all; close all; clc;

workingDir = pwd();

EndoSurfaceNodes = [];
SubjecIndex = 0;

subjectindex=[];
counter=1;
for i = 1:90;
    if i<10;
    if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV0',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
    end;
    else;
        if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
        end;
    end;
end;

for i = 1:(numel(subjectindex)-1); %Extracting Endo surface nodes
    if subjectindex(i)<10
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV0',num2str(subjectindex(i)),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,1:numel(sendo)) = sendo;
    else
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(subjectindex(i)),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,1:numel(sendo)) = sendo;
    end
end

[sendo, sepi, abaqusInputData] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(subjectindex(numel(subjectindex))),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EndoSurfaceNodes(SubjecIndex,:) = [sendo];

for i=1:24
    if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\MI\Patient',num2str(i),'_prolatePara.mat'))==2;
    [sendo, sepi, ~]=load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\MI\Patient',num2str(i),'\earlyDiastole'), ...
                     workingDir);
             EndoSurfaceNodes=[EndoSurfaceNodes; sendo];
    end
end
%{
    load('/Users/alan/ownCloud/PCA_Data/case1-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint1=reshape(sendo',1,2896*3);
EndoSurfaceNodes=[EndoSurfaceNodes; AMYpoint1];
    load('/Users/alan/ownCloud/PCA_Data/case2-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint2=reshape(sendo',1,2896*3);
 EndoSurfaceNodes=[EndoSurfaceNodes; AMYpoint2];
     load('/Users/alan/ownCloud/PCA_Data/case3-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint3=reshape(sendo',1,2896*3);
EndoSurfaceNodes=[EndoSurfaceNodes; AMYpoint3];
    load('/Users/alan/ownCloud/PCA_Data/case4-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sendo=abaqusInputData.node(abaqusInputData.endoNodes,1:3);
AMYpoint4=reshape(sendo',1,2896*3);
 EndoSurfaceNodes=[EndoSurfaceNodes; AMYpoint4];
%}

EpiSurfaceNodes = [];
SubjecIndex = 0;

subjectindex=[];
counter=1;
for i = 1:90;
    if i<10;
    if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV0',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
    end;
    else;
        if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(i),'_prolatePara.mat')) == 2;
        subjectindex(counter)=i;
        counter=counter+1;
        end;
    end;
end;

for i = 1:(numel(subjectindex)-1); %Extracting Endo surface nodes
    if subjectindex(i)<10
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV0',num2str(subjectindex(i)),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EpiSurfaceNodes(SubjecIndex,1:numel(sendo)) = sepi;
    else
 [sendo, sepi, ~] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(subjectindex(i)),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EpiSurfaceNodes(SubjecIndex,1:numel(sepi)) = sepi;
    end
end

[sendo, sepi, abaqusInputData] = load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\Results\HV',num2str(subjectindex(numel(subjectindex))),'\earlyDiastole'), ...
                     workingDir);
SubjecIndex = SubjecIndex + 1;
EpiSurfaceNodes(SubjecIndex,:) = [sepi];

for i=1:24
    if exist(strcat('C:\Users\2026068l\ownCloud\PCA_Data\MI\Patient',num2str(i),'_prolatePara.mat'))==2;
    [sendo, sepi, ~]=load_LVGeoData_1(strcat('C:\Users\2026068l\ownCloud\PCA_Data\MI\Patient',num2str(i),'\earlyDiastole'), ...
                     workingDir);
             EpiSurfaceNodes=[EpiSurfaceNodes; sepi];
    end
end
%{
    load('/Users/alan/ownCloud/PCA_Data/case1-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sepi=abaqusInputData.node(abaqusInputData.epiNodes,1:3);
AMYpoint1=reshape(sepi',1,2896*3);
EpiSurfaceNodes=[EpiSurfaceNodes; AMYpoint1];
    load('/Users/alan/ownCloud/PCA_Data/case2-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sepi=abaqusInputData.node(abaqusInputData.epiNodes,1:3);
AMYpoint2=reshape(sepi',1,2896*3);
 EpiSurfaceNodes=[EpiSurfaceNodes; AMYpoint2];
    load('/Users/alan/ownCloud/PCA_Data/case3-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sepi=abaqusInputData.node(abaqusInputData.epiNodes,1:3);
AMYpoint3=reshape(sepi',1,2896*3);
EpiSurfaceNodes=[EpiSurfaceNodes; AMYpoint3];
    load('/Users/alan/ownCloud/PCA_Data/case4-baseline/early-diastole/abaqusInputData_OneLayerMesh.mat');
    sepi=abaqusInputData.node(abaqusInputData.epiNodes,1:3);
AMYpoint4=reshape(sepi',1,2896*3);
 EpiSurfaceNodes=[EpiSurfaceNodes; AMYpoint4];
%}
 
SurfaceNodes=[EndoSurfaceNodes'; EpiSurfaceNodes']';

%SurfaceNodes=[EndoSurfaceNodes];
%SurfaceNodes=SurfaceNodes(1:94,:);

%%

%%calculate average shape
SurfaceNodes_ave = mean(SurfaceNodes);

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
[u, d, v] = svd(C, 0); %
% Pull out eigen values and vectors
eigVals = diag(d);
eigVecs = u(:,1:5);

%save('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian/Eigvecs_HVMI.mat','eigVecs')
%save('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian/SVD_withAMY.mat')

%load('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian/SVD.mat')
%%
load('C:\Users\2026068l\ownCloud\PhD\Papers\BHF\PCA_reconstructions\Cartesian\Eigvecs_HVMI.mat')
%%the cumulative energy content for the m'th eigenvector is the sum of the
%%energy content across eigenvalues 1: m
for i = 1 : N
    energy(i) = sum(eigVals(1:i));
end
propEnergy = energy./energy(end);

%%determin the number of principal components required to model 90% of data
%%variance
%percentMark = min(find(propEnergy > 0.9));
%percentMark = 1;
%eigVecsSelected = u(:,1:percentMark);
eigVecsSelected = eigVecs(:,1:5);

%%we can calculate all the weights 
Nweights = (SurfaceNodes_standardize)*eigVecsSelected;
%{
cd('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian')
bar(propEnergy(1:20))
title('Percentage Variation Explained')
ylabel('Variation Explained')
xlabel('Number of Principal Components')
print('/Users/alan/ownCloud/PhD/Reports/Report_2018_02_28_PCA_LDA/Cartesian_PC_variation_bar','-dpng')
%%next one we will need to output the mesh using tecplot
load('/Users/alan/ownCloud/Results/alpha0_global.mat')
alpha0 = alpha0_global;
tecFid = fopen('mean_mesh_endo.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,SurfaceNodes_ave, alpha0,'endo', tecFid);
fclose(tecFid);
tecFid = fopen('mean_mesh_epi.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,SurfaceNodes_ave(8689:17376), alpha0,'epi', tecFid);
fclose(tecFid);
tecFid = fopen('ori_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,[], alpha0,'endo', tecFid);
fclose(tecFid);
%}

%%
load('C:\Users\2026068l\ownCloud\PCA_Data\Results\alpha0_global.mat')
alpha0 = alpha0_global;
%%let us reconstruct the endo face *of one subject* with percentMark eigen vectors

cd('C:\Users\2026068l\ownCloud\PhD\Papers\BHF\PCA_reconstructions\Cartesian')
mapped_endo = zeros(size(eigVecsSelected(:,1)));
for i = 1:5
    mapped_endo = zeros(size(eigVecsSelected(:,1)));
    mapped_endo = mapped_endo-2*std(Nweights(:,i))*eigVecsSelected(:,i);  
    % mapped_endo = mapped_endo+faceWeight(i).*eigVecsSelected(:,i);
    tecFid = fopen(strcat('reconstructed_endo_minus2',num2str(i),'.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo(1:8688)+SurfaceNodes_ave(1:8688)', alpha0,'endo', tecFid);
fclose(tecFid);
    tecFid = fopen(strcat('reconstructed_epi_minus2',num2str(i),'.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo(8689:17376)+SurfaceNodes_ave(8689:17376)', alpha0,'epi', tecFid);
fclose(tecFid);
end


mapped_endo = zeros(size(eigVecsSelected(:,1)));
for i = 1:5
    mapped_endo = zeros(size(eigVecsSelected(:,1)));
    mapped_endo = mapped_endo+2*std(Nweights(:,i))*eigVecsSelected(:,i);  
    % mapped_endo = mapped_endo+faceWeight(i).*eigVecsSelected(:,i);
    tecFid = fopen(strcat('reconstructed_endo_plus2',num2str(i),'.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo(1:8688)+SurfaceNodes_ave(1:8688)', alpha0,'endo', tecFid);
fclose(tecFid);
    tecFid = fopen(strcat('reconstructed_epi_plus2',num2str(i),'.dat'),'w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo(8689:17376)+SurfaceNodes_ave(8689:17376)', alpha0,'epi', tecFid);
fclose(tecFid);
end

%{
size(u)

projection=EndoSurfaceNodes_standardize(1,:)*u(:,1:5)

mapped_endo = projection(1)*u(:,1)+projection(2)*u(:,2)+projection(3)*u(:,3)+projection(4)*u(:,4)+projection(5)*u(:,5);


tecFid = fopen('reconstructed_mesh.dat','w');
tectPlotFileGenerationForOneSurface(abaqusInputData,mapped_endo+EndoSurfaceNodes_ave', alpha0,'endo', tecFid,1);
fclose(tecFid);
%}
%%
%noamy
%load('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian/SVD_endo_epi.mat')

newpoints=SurfaceNodes*u(:,1:3)
newpointsAMY=[SurfaceNodes(95,:); SurfaceNodes(96,:); SurfaceNodes(97,:); SurfaceNodes(98,:)]*u(:,1:3)
scatter3(newpoints(1:72,1),newpoints(1:72,2),newpoints(1:72,3),'blue')
hold on 
scatter3(newpoints(73:94,1),newpoints(73:94,2),newpoints(73:94,3),'red')
hold off
hold on 
scatter3(newpointsAMY(:,1),newpointsAMY(:,2),newpointsAMY(:,3),'green')
hold off
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
view(0,0)

%saveas(gcf,'/Users/alan/ownCloud/PhD/Reports/Report_2018_02_21_PCA_LDA/Amyloidosis_PCA/3_PC_Projection_noAMY_PC13.png');
%%
%withamy
load('/Users/alan/ownCloud/PhD/Papers/BHF/PCA_reconstructions/Cartesian/SVD_withAMY.mat')

newpoints=SurfaceNodes*u(:,1:3)

scatter3(newpoints(1:72,1),newpoints(1:72,2),newpoints(1:72,3),'blue')
hold on 
scatter3(newpoints(73:94,1),newpoints(73:94,2),newpoints(73:94,3),'red')
hold off
hold on 
scatter3(newpoints(95,1),newpoints(95,2),newpoints(95,3),'green')
scatter3(newpoints(96,1),newpoints(96,2),newpoints(96,3),'black')
hold off
view(0,0)
%slice 21
    [sendo, sepi, ~]=load_LVGeoData_1('/Users/alan/ownCloud/PCA_Data/case4-baseline/early-diastole', ...
                     workingDir);
             AMY1_slice21=[sendo];
         [sendo, sepi, ~]=load_LVGeoData_1('/Users/alan/ownCloud/PCA_Data/case4-followup-9months/early-diastole', ...
                     workingDir);
             AMY2_slice21=[sendo];
             
AMY1_slice21points=AMY1_slice21*u(:,1:3);
AMY2_slice21points=AMY2_slice21*u(:,1:3);

scatter3(newpoints(1:72,1),newpoints(1:72,2),newpoints(1:72,3),'blue')
hold on 
scatter3(newpoints(73:94,1),newpoints(73:94,2),newpoints(73:94,3),'red')
hold off
hold on 
scatter3(newpoints(95,1),newpoints(95,2),newpoints(95,3),'green')
scatter3(newpoints(96,1),newpoints(96,2),newpoints(96,3),'black')
hold off
hold on 
scatter3(AMY1_slice21points(1),AMY1_slice21points(2),AMY1_slice21points(3),'black')
scatter3(AMY2_slice21points(1),AMY2_slice21points(2),AMY2_slice21points(3),'black')
hold off

saveas(gcf,'/Users/alan/ownCloud/PhD/Reports/MPPCA_2018_02_12/Cartesian_PCA/3_PC_Projection_AMY_slice21.png');
