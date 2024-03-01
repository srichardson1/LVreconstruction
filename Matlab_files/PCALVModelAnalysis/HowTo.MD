
## LV geometry reconstruction 

### LVWM_config    : config file, the one needs to be consistent for each running
  
adding config file in subdirectory projectFiles



### LVWM_DicomSampleSelection:  choose images

 if the images have been loaded previously, then no need to do it again, you can skip this and go to LVWM_SASegManualUsingImpoint



### LVWM_SASegManualUsingImpoint:           short-axis segmentation
* DataSegSA(imIndex).endo_cReal 
* DataSegSA(imIndex).epi_cReal 
* save the boundary points starting from x-axis and toward y+, x- and y- to itself

### LVWM_LASegManual: 

### LVWM_Apex_LASegManual

### SAAdjustmentBYLASeries: along 3 long axis, 2 itereations

### All_fitting_reconstruction

* LVWM_SALASeg_LongAxisAlignment
* LVWM_SA_GuidPointsGeneration_LongAxisAlignment
* LV_EndoFitting
* LV_EpiFitting
* LV_WholeMesh



## abaqusFilePreparation and myofiber construction

 run ABAQUS_FilePreparation.m

## AHA definition 

### LVDivision_AHADefinition:      
AHA division, current only using one mid slice and one apex slice, later needs to update, the whole LV also will be divied according to the definition. 
* LVWM_Division_AHADefinition 
* LVWM_AHAApply

## Segment MI regions and interpolated it into the mesh for simulation
see the folder MIMapping for details


## mainAPPStrainAnalysis_splineDeform: BSplineStrain recovery, needs to change slice manually since it may not work for all slices
* copy deformRes.dat (the guidance folder) to the earlyDiastole folder, update it according to the patient data 
* run LVWM_StrainAnalysis_splineDeform.m for all basal and middle slices 
* LVWM_StrainAnalysis_summarize_withRadialStrainOutput.m


