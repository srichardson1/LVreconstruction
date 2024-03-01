%%%SA segmentation
clear all;
close all;
clc;

% LVWM_config;
LVWM_config;
segB = 1;

cd(resultDir);

load imDesired;
cd(workingDir); 

cd(resultDir);
load DataSegSAOri;
cd(workingDir); 

%figure
%for imIndex = 1:1
%    endo_c = DataSegSA(imIndex).endo_cReal;
%    epi_c = DataSegSA(imIndex).epi_cReal;
%   plot 3D curves
%    plot3(endo_c(1,:),endo_c(2,:), endo_c(3,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
%    hold on
%    plot3(epi_c(1,:),epi_c(2,:), epi_c(3,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
%    hold off
%end
%
%   estimeate the centre
%
    imIndex=1;
    endo_c = DataSegSA(imIndex).endo_c;
    epi_c = DataSegSA(imIndex).epi_c;
    centrex=0.5*(mean(endo_c(1,:)+epi_c(1,:)));
    centrey=0.5*(mean(endo_c(2,:)+epi_c(2,:)));


totalSXSliceLocation = size(SXSliceSorted,2);
timeInstanceSelected = patientConfigs.timeInstanceSelected;
sampleN = patientConfigs.sampleN;


    for imIndex = 1 : 1
        imData =  SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imData;
        imInfo1 = SXSliceSorted(1,imIndex).SXSlice(timeInstanceSelected).imInfo;
        imInfo = infoExtract(imInfo1);
        sliceLocationStr = sprintf('%s',imInfo.SliceLocation);
        [imCropData, rect] = imcrop(imData,[]);

        h1=figure();
        imshow(imCropData,[]);hold on;
   end
        [xa1, ya1] = ginput(1);
        plot(xa1, ya1,'r+');
        [xa2, ya2] = ginput(1);
        plot(xa2, ya2,'r+');
        xa1=xa1+rect(1);
        ya1=ya1+rect(2);
        xa2=xa2+rect(1);
        ya2=ya2+rect(2);        
        
        dx=xa1-centrex;
        dy=ya1-centrey;
      
      if dx>=0 & dy==0
          theta1=0;
      elseif dx>0 & dy<0
          theta1 = -2*p+atan(abs(dy/dx));
      elseif dx==0 & dy>0
          theta1 = -pi/2;
      elseif dx<0 & dy<0
          theta1 = -pi-atan(abs(dy/dx));
      elseif dx<0 & dy==0
          theta1=pi;
      elseif dx<0 & dy>0
          theta1=-pi+atan(abs(dy/dx));
      elseif dx==0 & dy>0
          theta1=-1.5*pi;
      elseif dx>0 & dy>0
          theta1=-atan(abs(dy/dx));
      else
          theta1=0;
      end
      
        dx=xa2-centrex;
        dy=ya2-centrey;
      
      if dx>=0 & dy==0
          theta2=0;
      elseif dx>0 & dy<0
          theta2 = -2*pi+atan(abs(dy/dx));
      elseif dx==0 & dy>0
          theta2 = -pi/2;
      elseif dx<0 & dy<0
          theta2 = -pi-atan(abs(dy/dx));
      elseif dx<0 & dy==0
          theta2=pi;
      elseif dx<0 & dy>0
          theta2=-pi+atan(abs(dy/dx));
      elseif dx==0 & dy>0
          theta2=-1.5*pi;
      elseif dx>0 & dy>0
          theta2=-atan(abs(dy/dx));
      else
          theta2=0;
      end
        
%      theta1=theta1*180/pi;
%      theta2=theta2*180/pi;
      
      x(1)=xa1;
      y(1)=ya1;
      x(2)=centrex;
      y(2)=centrey;
      x(3)=xa2;
      y(3)=ya2;
      
%      figure
%      plot(endo_c(1,:),endo_c(2,:),'LineStyle', '-', 'Color', 'b', 'LineWidth',2);
%      hold on
%      plot(epi_c(1,:),epi_c(2,:),'LineStyle', '-', 'Color', 'r', 'LineWidth',2);
%      plot(centrex,centrey,'mx','MarkerSize', 10)
%      plot(xa1,ya1,'bs','MarkerSize', 10)
%      plot(xa2,ya2,'kd','MarkerSize', 10)
%      plot(x,y,'LineStyle', '--', 'Color', 'r', 'LineWidth',2)
%      legend('Endo-contour','Epi-contour','Centre','Start-RV','End-RV')
%      hold off

      h1=figure();
      imshow(imData,[]); hold on;
      title('Definition of common segment between RV and LV');
      plot(endo_c(1,:),endo_c(2,:),'b');
      plot(epi_c(1,:),epi_c(2,:),'r');
      plot(centrex,centrey,'mo','MarkerSize',10)
      plot(xa1,ya1,'bd','MarkerSize',8,'Color','g')
      plot(xa2,ya2,'ks','MarkerSize',8,'Color','g')
      plot(x,y,'LineStyle','--','Color','r','LineWidth',1.5)
      legend('Endo-contour','Epi-contour','Centre','Start-RV','End-RV')

      RVstartendangleWL(1).start = theta1;
      RVstartendangleWL(1).end = theta2;
%
%     make a data file for RV angles
%
      cd(resultDir);
      save RVstartendangleWL;
      cd(workingDir); 

    










