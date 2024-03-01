%%%this is used to divide LV into different regions
clear all; 
close all; 
clc;

% LVWM_config;
LVWM_config;
cd(resultDir);
load AHADefinition;
load abaqusInputData;
load XYZEndoFitted;
cd(workingDir);

%%here is the special treatment for HV25, because the slice 5 and slice 4
%%are same, need to interpolate 5 from 4 and 6
% if strcmp(patientConfigs(patientIndex,1).name, 'HV25')
%     u_slices_4 = XYZEndoFitted(1,4).u;
%     u_slices_6 = XYZEndoFitted(1,6).u;
%     
%     u_slices_5 = (u_slices_4 + u_slices_6)./2.0;
%     XYZEndoFitted(1,5).u = u_slices_5;
%     
% end


%% that is related to a specific patient, but considering it is general, 
%% basal slice: [1 2]; middle slice [3 4 5]; apical slices [6 7]
totalSASlices = length(basalSlices) + length(middlSlices) + length(apicaSlices);
slices_with_6regions = [basalSlices middlSlices];
slices_with_4regions = apicaSlices;
SASlicePositionApex = length(slices_with_6regions) + 1; %%the start index for apical slices

for i = 1 : totalSASlices
   u_slices(i) = mean(XYZEndoFitted(1,i).u); %%this is in radian, from 0 to -pi/2
end

%%now setup the degree boundary for each slices 
for i = 1 : totalSASlices
   if i == 1
       %u_up = u_slices(i);
       u_up = max([XYZEndoFitted(1,i).u, pi/2]);%make sure the top plane is always included
       u_bo = (u_slices(i+1)+u_slices(i))/2;
       %u_bo = min(XYZEndoFitted(1,i).u);
       u_slices_bc(i,:)=[u_up u_bo]*180/pi;
   elseif i==totalSASlices
       u_up = (u_slices(i)+u_slices(i-1))/2;
       u_bo = u_slices(i) + (u_slices(i)-u_slices(i-1))/2; 
       u_slices_bc(i,:)=[u_up u_bo]*180/pi;
   else
       u_up = (u_slices(i)+u_slices(i-1))/2;
       u_bo = (u_slices(i+1)+u_slices(i))/2;
       u_slices_bc(i,:)=[u_up u_bo]*180/pi;
   end 
end

% %%this is a specific setting up for subject 18
%  u_slices_bc(1,:) = [ 90.0000   -6.0];
%  u_slices_bc(2,:) = [ -6.0  -11.7194];

%%this part is for generating sliceRegions which will be used for LV
%%optimizaiton procedure
nodeMat=abaqusInputData.node;
%%this is for slice above apex region 
for sliceIndex = 1 : SASlicePositionApex - 1
   u_up = u_slices_bc(sliceIndex,1);
   u_bo = u_slices_bc(sliceIndex,2);
   segRegions = zeros([1,size(nodeMat, 1)]);
    k6 = find(slices_with_6regions==sliceIndex);
    AHAMidConfig = BaseMidConfig(1,k6);   
   for nodeIndex = 1 : size(nodeMat, 1)
      u_node = nodeMat(nodeIndex,4);
      if u_node > u_bo && u_node <= u_up
        theta = nodeMat(nodeIndex,5);%%in degree the circumferential degree
        regionValue = assignSegAccordingToThetaForMiddleRegion(theta,AHAMidConfig);
        segRegions(nodeIndex)=regionValue;
      end
     
   end
   sliceRegions(sliceIndex).segRegions = segRegions;
end
%%this is for slice in apex region 
if ~isempty(apicaSlices)
    for sliceIndex = SASlicePositionApex : apicaSlices(end)
        u_up = u_slices_bc(sliceIndex,1);
        u_bo = u_slices_bc(sliceIndex,2);
        segRegions = zeros(1,size(nodeMat, 1));
         for nodeIndex = 1 : size(nodeMat, 1)
            u_node = nodeMat(nodeIndex,4);
           if u_node > u_bo && u_node <= u_up
              theta = nodeMat(nodeIndex,5);%%in degree the circumferential degree
              k4 = find(slices_with_4regions==sliceIndex);
              AHAApexConfig = ApexConfig(1,k4); 
              regionValue = assignSegAccordingToThetaForApexRegion(theta,AHAApexConfig);
              segRegions(nodeIndex)=regionValue;
           end
         end
         sliceRegions(sliceIndex).segRegions = segRegions;
    end
end

%%for apex point
segRegions = zeros([1 size(nodeMat,1)]);
for nodeIndex = 1 : size(nodeMat,1)
    if nodeMat(nodeIndex,4)<u_bo
        segRegions(nodeIndex)=7;
    end
end
sliceRegions(sliceIndex+1).segRegions = segRegions;

%%%now combine segRegions
segRegions = zeros([1 size(nodeMat,1)]);
for sliceIndex = 1 : size(sliceRegions,2)
    segRegionsT = sliceRegions(sliceIndex).segRegions;
    for nodeIndex = 1 : length(segRegionsT)
        if segRegionsT(nodeIndex)>0
            segRegions(nodeIndex) = segRegionsT(nodeIndex);
        end        
    end    
end


%%%output the regions
nodeMat  = abaqusInputData.node;
elemMat = abaqusInputData.elem;
NodeMat(:,2) = nodeMat(:,1)*abaqusInputData.scaleTomm;
NodeMat(:,3) = nodeMat(:,2)*abaqusInputData.scaleTomm;
NodeMat(:,4) = nodeMat(:,3)*abaqusInputData.scaleTomm;
uvw = nodeMat(:,4:6);
nodeMatMRI = rotationBackToMRICoordinateSystemt(NodeMat,resultDir);

uvw(:,1) = 0;
uvw(:,2) = segRegions;
uvw(:,3) = 0;

cd(resultDir)
fid_regions = fopen('sliceSegments.dat','w');
cd(workingDir);
TecplotHexMeshVec(nodeMatMRI, elemMat,uvw,fid_regions);
fclose(fid_regions);
%%%%%%%%%**************************stop here, bellow code is for generating
%%%%%%%%%AHA 17 segments which will be used for summarization and others

cd(resultDir);
save LVMeshSegDivisions sliceRegions;
cd(workingDir);


%%%first we will need to divide each node into different slices
node=abaqusInputData.node;
for ni = 1 : size(node,1)
    slice_index_node = 0;
    u_node = node(ni,4);%%in degree
    for sliceIndex = 1 : totalSASlices
        u_up = u_slices_bc(sliceIndex,1);
        u_bo = u_slices_bc(sliceIndex,2);
        if u_node > u_bo && u_node <= u_up
            slice_index_node = sliceIndex;
            break;
        end
    end
    if u_node < u_slices_bc(totalSASlices,2) %%adding the value for the region bellow the last slice, the apical region
        slice_index_node = totalSASlices+1;
    end
    slice_No_assigned_Node(ni,1) = slice_index_node;
end

%%now we know the node associating with which plane, then we can assign the
%%regions
segRegions=[];
for ni = 1 : size(node,1)
   sliceIndex = slice_No_assigned_Node(ni,1);
   k6 = find(slices_with_6regions==sliceIndex);
   k4 = find(slices_with_4regions==sliceIndex);
   theta = node(ni,5);%%in degree the circumferential degree
   if ~isempty(k6) %%6 region definition
       AHAMidConfig = BaseMidConfig(1,k6);
       regionValue = assignSegAccordingToThetaForMiddleRegion(theta,AHAMidConfig);
       segRegions(ni,1)=regionValue;
   end
   if~isempty(k4) %%4 region definition
      AHAApexConfig = ApexConfig(1,k4); 
      regionValue = assignSegAccordingToThetaForApexRegion(theta,AHAApexConfig);
      segRegions(ni,1)=regionValue;
   end
   if isempty(k6) && isempty(k4)
       regionValue = 7;
       segRegions(ni,1)=regionValue;
   end   
end

%%output for check
cd(resultDir);
load rotationConfig;
cd(workingDir);

nodeMat  = abaqusInputData.node;
elemMat = abaqusInputData.elem;
NodeMat(:,2) = nodeMat(:,1)*abaqusInputData.scaleTomm;
NodeMat(:,3) = nodeMat(:,2)*abaqusInputData.scaleTomm;
NodeMat(:,4) = nodeMat(:,3)*abaqusInputData.scaleTomm;
uvw = nodeMat(:,4:6);
nodeMatMRI = rotationBackToMRICoordinateSystemt(NodeMat,resultDir);

% cd(resultDir);
% fid = fopen('LVMesh_AHADefinition_PCAReconstructed.dat','w');
% cd(workingDir);
% uvw(:,1) = slice_No_assigned_Node;
% uvw(:,2) = segRegions;
% uvw(:,3) = 0;
% TecplotHexMeshVec(nodeMatMRI, elemMat,uvw,fid);%ouput the v degree from the fitting
% fclose(fid);


%%%now we figure out the node id sequence for each region using AHA17
%%%definition
nodeRegions = [];
basa_InfSept=[];basa_AntSept=[];basa_Ant=[];base_AntLat=[];base_InfLat=[];base_Inf=[];
midd_InfSept=[];midd_AntSept=[];midd_Ant=[];midd_AntLat=[];midd_InfLat=[];midd_Inf=[];
apex_Sept=[];apex_Ant=[];apex_Lat=[];apex_Inf=[];
apicalRegion = [];
for ni = 1 : size(node,1)
   sliceIndex = slice_No_assigned_Node(ni,1);
   k_base = find(basalSlices==sliceIndex);
   k_mid = find(middlSlices==sliceIndex);
   k_apex = find(apicaSlices==sliceIndex);
   regionValue = segRegions(ni,1);
   
   if ~isempty(k_base)
       if regionValue == 1
           basa_InfSept = [basa_InfSept; ni];
           nodeRegions(ni,1) = 1;
       elseif regionValue == 2
           basa_AntSept = [basa_AntSept; ni];
           nodeRegions(ni,1) = 2;
       elseif regionValue == 3
           basa_Ant = [basa_Ant; ni];
           nodeRegions(ni,1) = 3;
       elseif regionValue == 4
           base_AntLat = [base_AntLat; ni];
           nodeRegions(ni,1) = 4;
       elseif regionValue == 5
           base_InfLat = [base_InfLat; ni];
           nodeRegions(ni,1) = 5;
       elseif regionValue == 6
            base_Inf = [base_Inf; ni];
            nodeRegions(ni,1) = 6;
       end   
   end %%k_base
   
   if ~isempty(k_mid)
       if regionValue == 1
           midd_InfSept = [midd_InfSept; ni];
           nodeRegions(ni,1) = 1+6;
       elseif regionValue == 2
           midd_AntSept = [midd_AntSept; ni];
           nodeRegions(ni,1) = 2+6;
       elseif regionValue == 3
           midd_Ant = [midd_Ant; ni];
           nodeRegions(ni,1) = 3+6;
       elseif regionValue == 4
           midd_AntLat = [midd_AntLat; ni];
           nodeRegions(ni,1) = 4+6;
       elseif regionValue == 5
           midd_InfLat = [midd_InfLat; ni];
           nodeRegions(ni,1) = 5+6;
       elseif regionValue == 6
            midd_Inf = [midd_Inf; ni];
            nodeRegions(ni,1) = 6+6;
       end   
   end %%k_base
   
   if ~isempty(k_apex)
       if regionValue == 1
           apex_Sept = [apex_Sept; ni];
           nodeRegions(ni,1) = 1+12;
       elseif regionValue == 2
           apex_Ant = [apex_Ant; ni];
           nodeRegions(ni,1) = 2+12;
       elseif regionValue == 3
           apex_Lat = [apex_Lat; ni];
           nodeRegions(ni,1) = 3+12;
       elseif regionValue == 4
           apex_Inf = [apex_Inf; ni];
           nodeRegions(ni,1) = 4+12;
       end
   end%%k_apex
   
   if regionValue == 7
       apicalRegion = [apicalRegion; ni];
       nodeRegions(ni,1) = 17;
   end
   
end

AHALVMeshDivision.basa_InfSept = basa_InfSept;
AHALVMeshDivision.basa_AntSept = basa_AntSept;
AHALVMeshDivision.basa_Ant     = basa_Ant;
AHALVMeshDivision.base_AntLat  = base_AntLat;
AHALVMeshDivision.base_InfLat  = base_InfLat;
AHALVMeshDivision.base_Inf     = base_Inf;

AHALVMeshDivision.midd_InfSept = midd_InfSept;
AHALVMeshDivision.midd_AntSept = midd_AntSept;
AHALVMeshDivision.midd_Ant     = midd_Ant;
AHALVMeshDivision.midd_AntLat  = midd_AntLat;
AHALVMeshDivision.midd_InfLat  = midd_InfLat;
AHALVMeshDivision.midd_Inf     = midd_Inf;

AHALVMeshDivision.apex_Sept = apex_Sept;
AHALVMeshDivision.apex_Ant = apex_Ant;
AHALVMeshDivision.apex_Lat = apex_Lat;
AHALVMeshDivision.apex_Inf = apex_Inf;

AHALVMeshDivision.apicalRegion = apicalRegion;
AHALVMeshDivision.nodeRegions = nodeRegions;


%% output to check whether it is right
cd(resultDir);
fid = fopen('LVMesh_AHADefinition.dat','w');
cd(workingDir);
uvw(:,1) = slice_No_assigned_Node;
uvw(:,2) = segRegions;
uvw(:,3) = nodeRegions;
TecplotHexMeshVec(nodeMatMRI, elemMat,uvw,fid);%ouput the v degree from the fitting
fclose(fid);



%%for element
%%%now we figure out the elem id sequence for each region using AHA17
%%%definition
elem = abaqusInputData.elem;
elem_basa_InfSept=[];elem_basa_AntSept=[];elem_basa_Ant=[];elem_base_AntLat=[];elem_base_InfLat=[];elem_base_Inf=[];
elem_midd_InfSept=[];elem_midd_AntSept=[];elem_midd_Ant=[];elem_midd_AntLat=[];elem_midd_InfLat=[];elem_midd_Inf=[];
elem_apex_Sept=[];elem_apex_Ant=[];elem_apex_Lat=[];elem_apex_Inf=[];
elem_apicalRegion = [];
elRegions = [];
for eli = 1 : size(elem,1)
   el_nodes_ID = elem(eli,:);
   sliceIndex_nodes = slice_No_assigned_Node(el_nodes_ID,1);
   sliceIndex = round(mean(sliceIndex_nodes));
   sliceIndex_max = max(sliceIndex_nodes);
   sliceIndex_min = min(sliceIndex_nodes);
   if sliceIndex_min > 0 %% try to use the miniumum sliceindex
           sliceIndex = sliceIndex_min;
   end
   
   node_coords_xyz = node(el_nodes_ID,1:3);
   node_coords_uvw = node(el_nodes_ID, 4:6);
   el_shapeCentre_xyz = [mean(node_coords_xyz(:,1)) mean(node_coords_xyz(:,2)) mean(node_coords_xyz(:,3))];
       
   el_shapeCentre_uvw = [mean(node_coords_uvw(:,1)) mean(node_coords_uvw(:,2)) mean(node_coords_uvw(:,3))];
   %% there will be difficulty of dealing with element across the zero
   el_shape_v_max = max(node_coords_uvw(:,2));
   el_shape_v_min = min(node_coords_uvw(:,2));
   el_cross_zero = 0;
   if el_shape_v_max - el_shape_v_min > 180 %%decide whether cross zero point
           el_cross_zero = 1;
   end
   
   k_base = find(basalSlices==sliceIndex);
   k_mid = find(middlSlices==sliceIndex);
   k_apex = find(apicaSlices==sliceIndex);
   
   %%initialize to be 17
   elRegions(eli, 1) = 17;
   
   
           
   if ~isempty(k_base)
       theta = el_shapeCentre_uvw(2);
       regionValue = assignSegAccordingToThetaForMiddleRegion(theta,AHAMidConfig);
       if el_cross_zero 
               regionValue_min = assignSegAccordingToThetaForMiddleRegion(el_shape_v_min,AHAMidConfig);
               regionValue_max = assignSegAccordingToThetaForMiddleRegion(el_shape_v_max,AHAMidConfig);
               regionValue = min([regionValue_min, regionValue_max]);
       end
       if regionValue>= 0.5 && regionValue < 1.5
           elem_basa_InfSept = [elem_basa_InfSept; eli];
           elRegions(eli,1) = 1;
       elseif regionValue >=1.5 && regionValue<2.5 
           elem_basa_AntSept = [elem_basa_AntSept; eli];
           elRegions(eli,1) = 2;
       elseif regionValue >=2.5 && regionValue < 3.5
           elem_basa_Ant = [elem_basa_Ant; eli];
           elRegions(eli,1) = 3;
       elseif regionValue >=3.5 && regionValue < 4.5
           elem_base_AntLat = [elem_base_AntLat; eli];
           elRegions(eli,1) = 4;
       elseif regionValue >=4.5 && regionValue <5.5
           elem_base_InfLat = [elem_base_InfLat; eli];
           elRegions(eli,1) = 5;
       elseif regionValue >=5.5 && regionValue <=6.5
            elem_base_Inf = [elem_base_Inf; eli];
            elRegions(eli,1) = 6;
       end   
   end %%k_base
   
   if ~isempty(k_mid)
       theta = el_shapeCentre_uvw(2);
       regionValue = assignSegAccordingToThetaForMiddleRegion(theta,AHAMidConfig);
       if el_cross_zero 
               regionValue_min = assignSegAccordingToThetaForMiddleRegion(el_shape_v_min,AHAMidConfig);
               regionValue_max = assignSegAccordingToThetaForMiddleRegion(el_shape_v_max,AHAMidConfig);
               regionValue = min([regionValue_min, regionValue_max]);
       end
       if regionValue  < 1.5
           elem_midd_InfSept = [elem_midd_InfSept; eli];
           elRegions(eli,1) = 1+6;
       elseif regionValue >=1.5 && regionValue<2.5 
           elem_midd_AntSept = [elem_midd_AntSept; eli];
           elRegions(eli,1) = 2+6;
       elseif regionValue >=2.5 && regionValue < 3.5
           elem_midd_Ant = [elem_midd_Ant; eli];
           elRegions(eli,1) = 3+6;
       elseif regionValue >=3.5 && regionValue < 4.5
           elem_midd_AntLat = [elem_midd_AntLat; eli];
           elRegions(eli,1) = 4+6;
       elseif regionValue >=4.5 && regionValue <5.5
           elem_midd_InfLat = [elem_midd_InfLat; eli];
           elRegions(eli,1) = 5+6;
       elseif regionValue >=5.5 && regionValue <=6.5
            elem_midd_Inf = [elem_midd_Inf; eli];
            elRegions(eli,1) = 6+6;
       end   
   end %%k_base
   
   if ~isempty(k_apex)
       theta = el_shapeCentre_uvw(2);
       regionValue = assignSegAccordingToThetaForApexRegion(theta,AHAApexConfig);
       if el_cross_zero 
               regionValue_min = assignSegAccordingToThetaForApexRegion(el_shape_v_min,AHAApexConfig);
               regionValue_max = assignSegAccordingToThetaForApexRegion(el_shape_v_max,AHAApexConfig);
               regionValue = min([regionValue_min, regionValue_max]);
       end
       if regionValue < 1.5 
           elem_apex_Sept = [elem_apex_Sept; eli];
           elRegions(eli,1) = 1+12;
       elseif regionValue >=1.5 && regionValue <2.5
           elem_apex_Ant = [elem_apex_Ant; eli];
           elRegions(eli,1) = 2+12;
       elseif regionValue >=2.5 && regionValue < 3.5
           elem_apex_Lat = [elem_apex_Lat; eli];
           elRegions(eli,1) = 3+12;
       elseif regionValue >= 3.5 && regionValue <= 4.5
           elem_apex_Inf = [elem_apex_Inf; eli];
           elRegions(eli,1) = 4+12;
       end
   end%%k_apex
   
   regionValue_nodes = segRegions(el_nodes_ID,1);
   regionValue = (mean(regionValue_nodes));
   if regionValue > 6.5 || ~isempty(find(regionValue_nodes>6,1))
       elem_apicalRegion = [elem_apicalRegion; eli];
       elRegions(eli,1) = 17;
   end
   
end


%% now we tried to generate cell centre data to show the regions, whether it
%% is right or not
%% a rough check 
nodeElRegion = [];
for i = 1 : size(elRegions,1)
    nl = elemMat(i,:);
    nodeElRegion(nl,1) = elRegions(i,1);
end



cd(resultDir);
fid = fopen('LVMesh_AHADefinition.dat','w');
cd(workingDir);
uvw(:,1) = slice_No_assigned_Node;
uvw(:,2) = nodeElRegion;
uvw(:,3) = nodeRegions;
TecplotHexMeshVec(nodeMatMRI, elemMat,uvw,fid);%ouput the v degree from the fitting

fclose(fid);

% AHALVMeshDivision.basa_InfSept = basa_InfSept;
% AHALVMeshDivision.basa_AntSept = basa_AntSept;
% AHALVMeshDivision.basa_Ant     = basa_Ant;
% AHALVMeshDivision.base_AntLat  = base_AntLat;
% AHALVMeshDivision.base_InfLat  = base_InfLat;
% AHALVMeshDivision.base_Inf     = base_Inf;
% 
% AHALVMeshDivision.midd_InfSept = midd_InfSept;
% AHALVMeshDivision.midd_AntSept = midd_AntSept;
% AHALVMeshDivision.midd_Ant     = midd_Ant;
% AHALVMeshDivision.midd_AntLat  = midd_AntLat;
% AHALVMeshDivision.midd_InfLat  = midd_InfLat;
% AHALVMeshDivision.midd_Inf     = midd_Inf;
% 
% AHALVMeshDivision.apex_Sept = apex_Sept;
% AHALVMeshDivision.apex_Ant = apex_Ant;
% AHALVMeshDivision.apex_Lat = apex_Lat;
% AHALVMeshDivision.apex_Inf = apex_Inf;
% 
% AHALVMeshDivision.apicalRegion = apicalRegion;

%%for element based 
AHALVMeshDivision.elem_basa_InfSept = elem_basa_InfSept;
AHALVMeshDivision.elem_basa_AntSept = elem_basa_AntSept;
AHALVMeshDivision.elem_basa_Ant     = elem_basa_Ant;
AHALVMeshDivision.elem_base_AntLat  = elem_base_AntLat;
AHALVMeshDivision.elem_base_InfLat  = elem_base_InfLat;
AHALVMeshDivision.elem_base_Inf     = elem_base_Inf;

AHALVMeshDivision.elem_midd_InfSept = elem_midd_InfSept;
AHALVMeshDivision.elem_midd_AntSept = elem_midd_AntSept;
AHALVMeshDivision.elem_midd_Ant     = elem_midd_Ant;
AHALVMeshDivision.elem_midd_AntLat  = elem_midd_AntLat;
AHALVMeshDivision.elem_midd_InfLat  = elem_midd_InfLat;
AHALVMeshDivision.elem_midd_Inf     = elem_midd_Inf;

AHALVMeshDivision.elem_apex_Sept = elem_apex_Sept;
AHALVMeshDivision.elem_apex_Ant = elem_apex_Ant;
AHALVMeshDivision.elem_apex_Lat = elem_apex_Lat;
AHALVMeshDivision.elem_apex_Inf = elem_apex_Inf;

AHALVMeshDivision.elem_apicalRegion = elem_apicalRegion;
AHALVMeshDivision.elRegions = elRegions;

AHALVMeshDivision.segRegions = segRegions;
AHALVMeshDivision.slice_No_assigned_Node = slice_No_assigned_Node;


%% here we define region for each slice, not following AHA17 segments
%totalSASlices
%slices_with_6regions = [basalSlices middlSlices];
%slices_with_4regions = apicaSlices;

%%here we will need to combine the first and second slices together due to
%%(here is a question, do we need to?)
%%the manipulation, the first slice only have very small region. 
combine_1_2 = false;
if length(slices_with_6regions)<= 4
    disp('slices_with_6regions is less than 4 slices, need to check');
    pause;
else 
    combine_1_2 = true;
end

elRegionsFull = [];
nodeElRegionFull = [];
% ElListRegionFull

%if combine_1_2 
if 1  %not combine_1_2
    totalElRegion = 6*length(slices_with_6regions) + 4*length(slices_with_4regions)+ 1;
    for eli = 1 : size(elem,1)
       el_nodes_ID = elem(eli,:);
       sliceIndex_nodes = slice_No_assigned_Node(el_nodes_ID,1);
       sliceIndex = round(mean(sliceIndex_nodes));
       sliceIndex_max = max(sliceIndex_nodes);
       sliceIndex_min = min(sliceIndex_nodes);
       if sliceIndex_min > 0 %% try to use the miniumum sliceindex
           sliceIndex = sliceIndex_min;
       end
       
       node_coords_xyz = node(el_nodes_ID,1:3);
       node_coords_uvw = node(el_nodes_ID, 4:6);
       el_shapeCentre_xyz = [mean(node_coords_xyz(:,1)) mean(node_coords_xyz(:,2)) mean(node_coords_xyz(:,3))];
       
       el_shapeCentre_uvw = [mean(node_coords_uvw(:,1)) mean(node_coords_uvw(:,2)) mean(node_coords_uvw(:,3))];
       %% there will be difficulty of dealing with element across the zero
       el_shape_v_max = max(node_coords_uvw(:,2));
       el_shape_v_min = min(node_coords_uvw(:,2));
       el_cross_zero = 0;
       if el_shape_v_max - el_shape_v_min > 180 %%decide whether cross zero point
           el_cross_zero = 1;
       end
       
       if sliceIndex == 0
           %disp('in apex');
       elseif sliceIndex > totalSASlices
           disp('in apex');
       end

       k_base = find(basalSlices==sliceIndex);
       k_mid =  find(middlSlices==sliceIndex);
       k_apex = find(apicaSlices==sliceIndex);

       %%initialize to be 17
       elRegionsFull(eli, 1) = -1;

       %regionValue_nodes = segRegions(el_nodes_ID,1);
       %regionValue = (mean(regionValue_nodes));
       
       sliceIndexT = 0;
       if ~isempty(k_base) || ~isempty(k_mid) 
           if sliceIndex == 1 || sliceIndex == 2
               sliceIndexT = sliceIndex;
           else
               sliceIndexT = sliceIndex;
           end
           theta = el_shapeCentre_uvw(2);
           regionValue = assignSegAccordingToThetaForMiddleRegion(theta,AHAMidConfig);
           if el_cross_zero 
               regionValue_min = assignSegAccordingToThetaForMiddleRegion(el_shape_v_min,AHAMidConfig);
               regionValue_max = assignSegAccordingToThetaForMiddleRegion(el_shape_v_max,AHAMidConfig);
               regionValue = min([regionValue_min, regionValue_max]);
           end
        
               if regionValue>= 0.0 && regionValue < 1.5
        %            elem_basa_InfSept = [elem_basa_InfSept; eli];
                   elRegionsFull(eli,1) = 1 + (sliceIndexT-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               elseif regionValue >=1.5 && regionValue<2.5 
        %            elem_basa_AntSept = [elem_basa_AntSept; eli];
                   elRegionsFull(eli,1) = 2 + (sliceIndexT-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               elseif regionValue >=2.5 && regionValue < 3.5
        %            elem_basa_Ant = [elem_basa_Ant; eli];
                   elRegionsFull(eli,1) = 3 + (sliceIndexT-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               elseif regionValue >=3.5 && regionValue < 4.5
        %            elem_base_AntLat = [elem_base_AntLat; eli];
                   elRegionsFull(eli,1) = 4+(sliceIndexT-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               elseif regionValue >=4.5 && regionValue <5.5
        %            elem_base_InfLat = [elem_base_InfLat; eli];
                   elRegionsFull(eli,1) = 5+(sliceIndexT-1)*6;
                   nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               elseif regionValue >=5.5 && regionValue <=6.5
        %             elem_base_Inf = [elem_base_Inf; eli];
                    elRegionsFull(eli,1) = 6+(sliceIndexT-1)*6;
                    nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
               end  
               
       end %%k_base and k_mid

       if ~isempty(k_apex)
           theta = el_shapeCentre_uvw(2);
           regionValue = assignSegAccordingToThetaForApexRegion(theta,AHAApexConfig);
           if el_cross_zero 
               regionValue_min = assignSegAccordingToThetaForApexRegion(el_shape_v_min,AHAApexConfig);
               regionValue_max = assignSegAccordingToThetaForApexRegion(el_shape_v_max,AHAApexConfig);
               regionValue = min([regionValue_min, regionValue_max]);
           end
           if regionValue < 1.5 
    %            elem_apex_Sept = [elem_apex_Sept; eli];
               elRegionsFull(eli,1) = 1+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
           elseif regionValue >=1.5 && regionValue <2.5
    %            elem_apex_Ant = [elem_apex_Ant; eli];
               elRegionsFull(eli,1) = 2+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
           elseif regionValue >=2.5 && regionValue < 3.5
    %            elem_apex_Lat = [elem_apex_Lat; eli];
               elRegionsFull(eli,1) = 3+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
           elseif regionValue >= 3.5 && regionValue <= 4.5
    %            elem_apex_Inf = [elem_apex_Inf; eli];
               elRegionsFull(eli,1) = 4+(length(slices_with_6regions))*6+ (sliceIndex- length(slices_with_6regions)-1)*4;
               nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
           end
       end%%k_apex

       if sliceIndex > totalSASlices 
    %        elem_apicalRegion = [elem_apicalRegion; eli];
           elRegionsFull(eli,1) = 6*(length(slices_with_6regions)) + 4*length(slices_with_4regions)+ 1;
           nodeElRegionFull(el_nodes_ID,1) = elRegionsFull(eli,1);
       end

    end
end %combine_1_2



%% now we tried to generate cell centre data to show the regions, whether it
%% is right or not
%% a rough check 
% nodeElRegionFull = [];
% for i = 1 : size(elRegionsFull,1)
%     nl = elemMat(i,:);
%     nodeElRegionFull(nl,1) = elRegionsFull(i,1);
% end



cd(resultDir);
fid = fopen('LVMesh_AHASegments_full.dat','w');
cd(workingDir);
uvw(:,1) = slice_No_assigned_Node;
uvw(:,2) = nodeElRegionFull;
uvw(:,3) = nodeRegions;
TecplotHexMeshVec(nodeMatMRI, elemMat,uvw,fid);%ouput the v degree from the fitting
fclose(fid);



AHALVMeshDivision.elRegionsFull = elRegionsFull;


cd(resultDir);
save AHALVMeshDivision AHALVMeshDivision;
cd(workingDir);





