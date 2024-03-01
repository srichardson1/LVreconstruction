
%clear all; close all; clc;


%%%find the proper directory 
% LVWM_config;%%hao yu jue


%LVWM_config;

%%%% load parameter for one study

%cd(resultDir);

%SetParameters;

%load abaqusInputData;
%cd(workingDir);

nodeMat = abaqusInputData.node;
endoNodes = (abaqusInputData.endoNodes)';
epiNodes = (abaqusInputData.epiNodes)';

EndoNodesXYZ = nodeMat(endoNodes,:);
EpiNodesXYZ  = nodeMat(epiNodes,:);

DT_endo =  DelaunayTri(EndoNodesXYZ(:,1),EndoNodesXYZ(:,2),EndoNodesXYZ(:,3));
[~, vol_endo] = convexHull(DT_endo);

DT_epi =  DelaunayTri(EpiNodesXYZ(:,1),EpiNodesXYZ(:,2),EpiNodesXYZ(:,3));
[~, vol_epi] = convexHull(DT_epi);

vol_LVwall = vol_epi - vol_endo;



