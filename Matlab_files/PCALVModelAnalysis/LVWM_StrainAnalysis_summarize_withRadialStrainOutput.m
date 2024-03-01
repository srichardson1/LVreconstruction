clear all; close all; clc;
%%this will post-processing the strain data in order to decide which we
%%should use and which we should not 
% LVWM_config;
LVWM_config;
path(path, './demonCode');

cd(resultDir);
fid = fopen('deformRes.dat', 'r');
fidS = fopen('deformResAllSlices.dat', 'w');
fidES = fopen('systolicStrain.dat', 'w');
%%for radial strain
fidS_rad = fopen('deformResAllSlices_rad.dat', 'w');
fidES_rad = fopen('systolicStrain_rad.dat', 'w');
cd(workingDir);

%%%in the sequence of infsetpal antseptal, ant, antlat, inflat, inf
%%%currently only consider circumferential
tline = fgetl(fid);
while ~feof(fid)
    cd(resultDir);
    cd(tline(1:end-1));
    tline = fgetl(fid);
    sliceNo = sscanf(tline, '%d');
    load BsplineResult_slice;
    cd(workingDir);
    
    index_cirInfSeptTotal = find_next_index_near_min_strain(cirInfSeptTotal);
    index_cirAntSeptTotal = find_next_index_near_min_strain(cirAntSeptTotal);

    index_cirAntTotal = find_next_index_near_min_strain(cirAntTotal);
    index_cirAntLatTotal = find_next_index_near_min_strain(cirAntLatTotal);
    
    index_cirInfLatTotal = find_next_index_near_min_strain(cirInfLatTotal);
    index_cirInfTotal = find_next_index_near_min_strain(cirInfTotal);
    
    %%for radial strain
    index_radInfSeptTotal = find_next_index_near_max_rad_strain(radInfSeptTotal);
    index_radAntSeptTotal = find_next_index_near_max_rad_strain(radAntSeptTotal);

    index_radAntTotal = find_next_index_near_max_rad_strain(radAntTotal);
    index_radAntLatTotal = find_next_index_near_max_rad_strain(radAntLatTotal);
    
    index_radInfLatTotal = find_next_index_near_max_rad_strain(radInfLatTotal);
    index_radInfTotal = find_next_index_near_max_rad_strain(radInfTotal);
    
    
%     fprintf(fidS, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, min(cirInfSeptTotal), min(cirAntSeptTotal), ...
%                     min(cirAntTotal), min(cirAntLatTotal), min(cirInfLatTotal), min(cirInfTotal));
    fprintf(fidS, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, mean(cirInfSeptTotal(index_cirInfSeptTotal)), mean(cirAntSeptTotal(index_cirAntSeptTotal)), ...
                    mean(cirAntTotal(index_cirAntTotal)), mean(cirAntLatTotal(index_cirAntLatTotal)), mean(cirInfLatTotal(index_cirInfLatTotal)), mean(cirInfTotal(index_cirInfTotal)));
    fprintf(fidS, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, 1,1,1,1,1,1);
    fprintf(fidS, '\n');
    fprintf(fidS, '\n');
    
    %%for end-systolic strain estimation 
    fprintf(fidES, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, min(cirInfSeptTotal), min(cirAntSeptTotal), ...
                    min(cirAntTotal), min(cirAntLatTotal), min(cirInfLatTotal), min(cirInfTotal));
    fprintf(fidES, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, 1,1,1,1,1,1);
    fprintf(fidES, '\n');
    fprintf(fidES, '\n');
    
    
    %%%output for rad strain 
    %     fprintf(fidS, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, min(cirInfSeptTotal), min(cirAntSeptTotal), ...
%                     min(cirAntTotal), min(cirAntLatTotal), min(cirInfLatTotal), min(cirInfTotal));
    fprintf(fidS_rad, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, mean(radInfSeptTotal(index_radInfSeptTotal)),...
                    mean(radAntSeptTotal(index_cirAntSeptTotal)), ...
                    mean(radAntTotal(index_cirAntTotal)), mean(radAntLatTotal(index_cirAntLatTotal)), ...
                    mean(radInfLatTotal(index_cirInfLatTotal)), mean(radInfTotal(index_cirInfTotal)));
    fprintf(fidS_rad, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, 1,1,1,1,1,1);
    fprintf(fidS_rad, '\n');
    fprintf(fidS_rad, '\n');
    
    %%for end-systolic strain estimation 
    fprintf(fidES_rad, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, max(radInfSeptTotal), max(radAntSeptTotal), ...
                    max(radAntTotal), max(radAntLatTotal), max(radInfLatTotal), max(radInfTotal));
    fprintf(fidES_rad, '%d,\t %f,\t %f,\t %f,\t %f,\t %f,\t %f,\t\n', sliceNo, 1,1,1,1,1,1);
    fprintf(fidES_rad, '\n');
    fprintf(fidES_rad, '\n');
    
    
    tline = fgetl(fid);
end

fclose(fidS);
fclose(fidES);
fclose(fid);



%%%summarize the average peak ES strain
cd(resultDir);
data_strain_rad = load('systolicStrain_rad.dat');
data_strain = load('systolicStrain.dat');
cd(workingDir);

N = size(data_strain,1)/2;
SegN = size(data_strain,2);

strainES = [];
strainES_rad = []; 
for i = 1 : N
    for j =  2 : SegN
        if data_strain(2*i,j) > 0.9 %&& data_strain(2*i-1,j) < -0.1
            strainES = [strainES; data_strain(2*i-1,j)];
        end
        if data_strain_rad(2*i,j) > 0
            strainES_rad = [strainES_rad; data_strain_rad(2*i-1,j)];
        end
    end
end

aveStrainEs = [mean(strainES)  mean(strainES_rad)]

