%%%this is used to divide LV into different regions
clear all; 
close all; 
clc;

% LVWM_config;%%%Hao yu jue
LVWM_config;

timeInstanceSelected = 1; %%%using diastole to define the degree 
SASlicePositonMiddle = floor(usuableSXSlice/2); %%%the middle slice position
SASlicePositonApex = SASlicePositionApex;
ImgDivisionManualB = 1;


path(path, './PostProcessingMatlabIB/IBpostWithMRI');

if ImgDivisionManualB == 1
	% cd(IBResultDir)
	% load HeartMesh;
	% cd(workingDir);
	%%%extrat the mesh with nodes and elems
	% elemMat = HeartMesh.elemMat;
	% nodeMat = HeartMesh.nodeMat;  %%%be aware that node in mm

	%%%now load the first sa image
	%%%now need to get the 3D dicom image
	cd(resultDir);
	load imDesired;
	load rotationConfig;
	cd(workingDir);

	for positionIndex = 1 : 2
		if positionIndex ==1 
			imFileName = SXSliceSorted(SASlicePositonMiddle).Time(timeInstanceSelectedDiastile).name;
		else
			imFileName = SXSliceSorted(SASlicePositonApex).Time(timeInstanceSelectedDiastile).name;
		end
		imFileName = sprintf('%s/%s',dicomDir,imFileName);
		imInfo1 = dicominfo(imFileName);
		imInfo = infoExtract(imInfo1);
		imData = dicomread(imFileName);


		h1=figure(); imshow(imData,[]);hold on;pause;
		%%define 6 points or 4 points 
		%%%starting from inferior insertion go through anterior insertion
		if positionIndex ==1 
			endo_c=[];
			but = 1;
			n=1;
			while but ==1 && n<=6
					[x y but]=ginput(1);
					if n==1 
						plot(x,y,'b.');hold on;
					elseif n==2
						plot(x,y,'r+'); hold on;
					elseif n==3
						plot(x,y,'y*');hold on;
					elseif n==4
						plot(x,y,'b<'); hold on;
					else
						plot(x,y,'r.');hold on;
					end
					endo_c(1,n)=x;
					endo_c(2,n)=y;
					n=n+1;
			end

			%%%degree calculation
			theta = [];
			centerPoint = [mean(endo_c(1,:)) mean(endo_c(2,:))];
			for i = 1 : size(endo_c,2)
				p = [endo_c(1,i),endo_c(2,i)];
				theta(i) = degreeCalculationPointBased(p,centerPoint)*180/pi; %%in the range of 0-360
			end

			%%%segment region
			MidConfig.InfSeptTheta = degreeReOrder(theta(2),theta(1));
			MidConfig.AntSeptTheta = degreeReOrder(theta(3),theta(2));
			MidConfig.AntTheta = degreeReOrder(theta(4),theta(3));
			MidConfig.AntLatTheta = degreeReOrder(theta(5),theta(4));
			MidConfig.InfLatTheta = degreeReOrder(theta(6),theta(5));
			MidConfig.InfTheta = degreeReOrder(theta(6),theta(1));
			MidConfig.endo_c = endo_c;
			MidConfig.theta = theta;
		else %%this is for the apex slices
			endo_c=[];
			but = 1;
			n=1;
			while but ==1 && n<=4
					[x y but]=ginput(1);
					if n==1 
						plot(x,y,'b.');hold on;
					elseif n==2
						plot(x,y,'r+'); hold on;
					elseif n==3
						plot(x,y,'y*');hold on;
					elseif n==4
						plot(x,y,'b<'); hold on;
					else
						plot(x,y,'r.');hold on;
					end
					endo_c(1,n)=x;
					endo_c(2,n)=y;
					n=n+1;
			end

			%%%degree calculation
			centerPoint = [mean(endo_c(1,:)) mean(endo_c(2,:))];
			theta = [];
			for i = 1 : size(endo_c,2)
				p = [endo_c(1,i),endo_c(2,i)];
				theta(i) = degreeCalculationPointBased(p,centerPoint)*180/pi; %%in the range of 0-360
			end

			%%%segment region
			ApexConfig.SeptTheta = degreeReOrder(theta(2),theta(1));
			ApexConfig.AntTheta = degreeReOrder(theta(3),theta(2));
			ApexConfig.Lat = degreeReOrder(theta(4),theta(3));
			ApexConfig.Inf = degreeReOrder(theta(1),theta(4)); 
			ApexConfig.endo_c = endo_c;
			ApexConfig.theta = theta;
		end
	end   

	cd(resultDir)
	save DivisionConfig MidConfig ApexConfig;
	cd(workingDir);
else
	cd(resultDir)
	load DivisionConfig;
	cd(workingDir);
end



cd(resultDir);
if exist('abaqusInputData.mat','file')
   cd(workingDir);
   LVMeshDivisionAHA(resultDir,MidConfig, ApexConfig, SASlicePositionApex, SASliceDistance,usuableSXSlice); 
end
cd(workingDir);




% h3D = figure(); hold on;set(h3D, 'renderer', 'opengl');
% tetramesh(elemMat(:,2:5), nodeMat(:, 2: 4),'FaceColor','r');

%  segPoints_cReal = TransformCurvesFromImToRealSpace(endo_c,imInfo);
% 
%  endo_cT(1,:) = segPoints_cReal(1,:)-LVUpperCenter(1);
%  endo_cT(2,:) = segPoints_cReal(2,:)-LVUpperCenter(2);
%  endo_cT(3,:) = segPoints_cReal(3,:)-LVUpperCenter(3);
%  segPoints_cRealRotated = RotationMatrix*endo_cT;
% 
% h3D = figure(); hold on;
% set(h3D, 'renderer', 'opengl');
% trisurf(HeartMesh.epiFace,nodeMat(:, 2), nodeMat(:, 3),nodeMat(:, 4), 'FaceColor', 'r', 'FaceAlpha',0.5, 'EdgeColor','none');
% view([0 90]);
% for i = 1 : size(segPoints_cRealRotated,2)
%     x = segPoints_cRealRotated(1,i);
%     y = segPoints_cRealRotated(2,i);
%     z = segPoints_cRealRotated(3,i);
%     if i==1 
%             plot3(x,y,z,'Color','b','Marker','+');hold on;
%     elseif i==2
%             plot3(x,y,z,'Color','k','Marker','+'); hold on;
%     elseif i==3
%             plot3(x,y,z,'Color','y','Marker','*');hold on;
%     elseif i==4
%             plot3(x,y,z,'Color','b','Marker','<'); hold on;
%     else
%             plot3(x,y,z,'Color','b','Marker','>');hold on;
%     end
%     
% end









