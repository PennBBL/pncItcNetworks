---
title: "PNC-ITC Replication"
authors: Kahini
output: html_document
---
#### Based off pncitc.md


 
# 1. CWASMDMR

The computation of  CWASMDMR was  done with  `cwasmdr` singularity image (`/cbica/projects/pncitc/cwasmdmr.simg`).

Distance matrix was first computed with the following script: 

```
#!/bin/bash
#$ -l h_vmem=320G
#$ -l tmpfree=200G

#SING=/share/apps/singularity/2.5.1/bin/singularity # singularity directory on chead 
singimage=/cbica/projects/pncitc/cwasmdmr.simg # cwasmdmr image on chead
scriptdir=/usr/local/bin # script to run inside the image 


mdmrouput=/cbica/projects/pncitc/mehtareplicate/cwas307 #output directory. note that i added /home/username in front 
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc

imagelist=/cbica/projects/pncitc/subjectData/imageinput_rest2.csv # list of image in nifti # check image list, /home/username  
#modelcsv=/cbica/home//adebimpa/cbica/projects/GURLAB/projects/pncitc/subjectData/demographics/model.csv # demographic and factors, 


rm  -rf $mdmrouput # remove previous run if exist 
metric=pearson # pearson correlation 

# compute distance matrix

singularity exec -e -B /cbica/projects/pncitc $singimage $scriptdir/Rscript $scriptdir/connectir_subdist.R \
        $mdmrouput \
	    --infuncs1=$imagelist \
	    --ztransform \
        --automask1  -c 3 -t 4 \
	    --brainmask1=$brainmask \
        --method="$metric"  \
        --bg=$bgim   --overwrite \
	    --memlimit=200

```

The output of distance matrix: `/cbica/projects/pncitc/mehtareplicate/cwas307`
   
The  distance matrix  was used for mdmr computation with `logk` as the main factor.
other covariates used are `sex`, `age`, and `relative rms`:

 ```math  
 distancematrix = f(logk)+relMeanRMSmotion+sex+age 
 ```
   
The script used for mdmr computation is as below: 
```
#!/bin/bash
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/mehtareplicate/cwas307 -f 'logk+relMeanRMSmotion+sex+age' -m /cbica/projects/pncitc/demographics/n307_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```
 
# The following will be updated as the steps are replicated
### 2. Significant clusters from mdmr
The cluster analysis was computed  with the script `scripts/grf_fslcluster.sh`, written based on  [FSL cluster analysis](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster) with  Gaussian Random Field (GRF) theory. 

The script `scripts/clusterz3.0.9.sh` was computed on the z-value obtained from mdmr  with the threshold `z=3.09`
Two clusters was obtained: one at the at the frontal region the order at the TPJ. 



The output of cluster masks: `/cbica/projects/GURLAB/projects/pncitc/output/cluster_Z3.09`




### 3. Seed-based correlation 
The two  masks were upsample from 4mm to 2mm and were used as seeds for seed-based correlation.

the two seeds: 
`/cbica/projects/GURLAB/projects/pncitc/output/cluster_Z3.09/mask1_2mm.nii.gz `
`/cbica/projects/GURLAB/projects/pncitc/output/cluster_Z3.09/mask2_2mm.nii.gz` 

The seed-based correlation was computed with the `scripts/seedcorrelations.sh`

### 4. Linear regression with FSL `flameo` 

Flameo regression computation requires `design`,`contrast` and `group`. The script `scripts/makeflameodesig.R`

The flameo linear regression was computed with the script: `scripts/flameo`
The outputs of flameo regression: 
`/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/mask1`

`/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/mask2`

In each directory, there are zvalues: 

  `zstat1 : average `

  `zstat2 : logk `

  `zstat3 : sex`

  `zstat4 : motion`

  `zstat5 : age`


The zstats  were FDR corrected with this script `scripts/flameoutputfdrcorrection`. The outputs are located here:

`/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/mask1/logk`
`/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/mask2/logk`

FDR corrected z-values. 

  `zfdr1 : average `

  `zfdr2 : logk `

  `zfdr3 : sex`

  `zfdr4 : motion`

  `zfdr5 : age`

### 5. Vizualisation of Results

All computations were done in PNC template. For vizualisation, all the nifti files  were tranformed to MNI before. 

  a. for clusters and mean of seed-based correlation : `notebook/seed-basedcorrelation.ipynb`

  b. for mask1 : `notebook/flameomask1.ipynb`

  c. for mask2 : `notebook/flameomask2.ipynb`




### 6. Regional plot of signitiifcant regions of logk 
 
 The positive and negative zvalues of seed-based correlation regression with logk was extracted with the script `scripts/extractsignificantcluster.R` for both seed masks 

 The results were vizualised with `notebook/meanseedcorrelationplot.Rmd`





 
