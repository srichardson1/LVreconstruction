if ~exist(folder, 'dir')

mkdir(folder);

xlswrite(fullFileName,{'Patient id'},'Sheet1','A1')
xlswrite(fullFileName,{'Mean(Fcc)'},'Sheet1','B1')
xlswrite(fullFileName,{'Mean(Fcr)'},'Sheet1','C1')
xlswrite(fullFileName,{'Mean(Fcl)'},'Sheet1','D1')
xlswrite(fullFileName,{'Mean(Frc)'},'Sheet1','E1')
xlswrite(fullFileName,{'Mean(Frr)'},'Sheet1','F1')
xlswrite(fullFileName,{'Mean(Frl)'},'Sheet1','G1')
xlswrite(fullFileName,{'Mean(Fcl)'},'Sheet1','H1')
xlswrite(fullFileName,{'Mean(Frl)'},'Sheet1','I1')
xlswrite(fullFileName,{'Mean(Fll)'},'Sheet1','J1')
xlswrite(fullFileName,{'I1'},'Sheet1','K1')
xlswrite(fullFileName,{'I2'},'Sheet1','L1')
xlswrite(fullFileName,{'I3'},'Sheet1','M1')

end








