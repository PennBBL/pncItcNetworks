---
title: "PNC-ITC Replication"
authors: Kahini
output: html_document
---

Most recent version of manuscript available [here](https://docs.google.com/document/d/1zxtCFqp-8m46C0Al7LhRQeBWoTX5PBvJ/edit#)
#### Copied and modified from pncitc.md

Sample replication: n307 did not exclude those with the `health_exclude` criteria. Analyses were re-run on n293, using Pehlivanova et al's n427 sample and then running restQA exclusions on them (all this information was from Pehlivanova .csvs). The code and .csvs for this are available at `cbica/projects/pncitc/samplerecreation`. All final analyses were run in `cbica/projects/pncitc/ignore` (yes, pretty poorly named, but this was initially just a safety check), using the same steps as below  - additionally, any new .csvs should be pointed to in the scripts in that directory. I also moved the bblid_scanid .csv to demographics, and created a folder within subjectData called rest293 for the n = 293 replication.
**Results were similar in N293 for the second cluster, things changed for the first cluster. Visualizations are available in the .html format within ignore**

Additionally: "~/" means that the path must be supplied by the user. 

The code for sample replication is: 

```
setwd("/Users/kahinim/Desktop")

# read the subject demographics
restdatapnc=read.csv('ITC/demographics/n2416_RestQAData_20170714.csv') # pnc QA for resting-state data
nmel=read.csv('ITC/demographics/n452_pnc_itc_whole_sample_20160825.csv') # Marieta final subject QA  
z = read.csv('/Users/kahinim/Desktop/n427_fsSubcortVol.csv')
pncitc=merge(nmel,restdatapnc, by=c('bblid','scanid')) # merge by Ids  
pncitc=merge(pncitc,z, by=c('bblid','scanid')) # merge by Ids  
# select the neccessary variable for screening and further processing
# age, logk, sex, rest exclusion  variables: voxelwise and motion
pncit1 <- data.frame(
  pncitc$bblid,
  pncitc$healthExclude,
  pncitc$scanid,pncitc$logk,pncitc$ageAtScan,pncitc$logAlpha,pncitc$sex,pncitc$race,pncitc$race2,pncitc$restExclude,pncitc$restExcludeVoxelwise,
  pncitc$restNoDataExclude,pncitc$relMeanRMSmotion,pncitc$restNSpikesMotion,pncitc$restNSpikesMotionExclude,pncitc$restRelMeanRMSMotionExclude, pncitc$restNoDataExclude, pncitc$restVoxelwiseCoverageExclude
)
colnames(pncit1)=c('bblid','healthExclude',
                   'scanid','logk','ageAtScan','logAlpha','sex','race','race2','restExclude','restExcludeVoxelwise',
                   'restNoDataExclude','relMeanRMSmotion','restNSpikesMotion','restNSpikesMotionExclude','restRelMeanRMSMotionExclude','restNoDataExclude', 'restVoxelwiseCoverageExclude')

pncit1=pncit1[which(pncit1$restExcludeVoxelwise==0),]
pncit1=pncit1[which(pncit1$restNoDataExclude==0),]
pncit1=pncit1[which(pncit1$restRelMeanRMSMotionExclude==0),]
pncit1=pncit1[which(pncit1$restVoxelwiseCoverageExclude==0),]
pncit1=pncit1[which(pncit1$healthExclude==0),]
pncit1=pncit1[which(pncit1$restNSpikesMotionExclude==0),]
pncit1=pncit1[-which(is.na(pncit1$relMeanRMSmotion)),]


# during manual checking, one subject (id:96832) has data points 90 less than 120 expected 
pncit1=pncit1[-which(pncit1$bblid==96832),]

#get the ids of final subjects
ids=data.frame(pncit1$bblid,pncit1$scanid) # get bblid and scanid for futher analyses 

# write out demographics and bblid and scanid
write.csv(ids,'n293_blbid_scanid.csv',row.names = FALSE,quote = FALSE)
pncit1$age=pncit1$ageAtScan/12
write.csv(pncit1,'n293_demographics.csv',row.names = FALSE,quote = FALSE)
```

### 1. CWAS-MDMR

The computation of  CWASMDMR was  done with  `cwasmdr` singularity image (`/cbica/projects/pncitc/cwasmdmr.simg`). We used packages from the connectir project at [https://github.com/czarrar/connectir](https://github.com/czarrar/connectir)

Distance matrix was first computed with the following script: 

```
#!/bin/bash
#$ -l h_vmem=320G #QSUB, can take some hours
#$ -l tmpfree=200G
singimage=/cbica/projects/pncitc/cwasmdmr.simg 
scriptdir=/usr/local/bin
mdmrouput=/cbica/projects/pncitc/ignore/cwas293 #output directory
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc
imagelist=/cbica/projects/pncitc/subjectData/imageinput_rest3.csv #list of image in nifti # HAD TO RE-GENERATE THIS LIST AS THE FILE PATHS HAD CHANGED
rm  -rf $mdmrouput # remove previous run if exist 
metric=pearson # pearson correlation 
# compute distance matrix
singularity exec -e -B /cbica/projects/pncitc $singimage $scriptdir/Rscript $scriptdir/connectir_subdist.R $mdmrouput --infuncs1=$imagelist --ztransform --automask1  -c 3 -t 4 --brainmask1=$brainmask --method="$metric" --bg=$bgim  --overwrite --memlimit=200

```

The output of distance matrix: `/cbica/projects/pncitc/ignore/cwas293`
   
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
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/ignore/cwas293 -f 'logk+relMeanRMSmotion+sex+age' -m /cbica/projects/pncitc/samplerecreation/n293_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```

Memory and formatting were the main problems with these scripts not running well - the numbers/format left in were what worked for me. The output is at: `/cbica/projects/pncitc/ignore/cwas293/logk_motion_sex_age`

### 2. Significant clusters from mdmr
The cluster analysis was computed  with the script `scripts/grf_fslcluster.sh`, written based on  [FSL cluster analysis](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster) with  Gaussian Random Field (GRF) theory

The script `scripts/cluster.sh` called grf_fslcluster.sh, as listed in this repo under scripts, with `z=3.09`
Two clusters were obtained: one at the at the frontal region the order at the TPJ. 

cluster.sh reads as below:
```
 #!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/ignore/cwas293/logk_motion_sex_age/zstats_logk.nii.gz  -m ${dir}/ignore/cwas293/mask.nii.gz -t 3.09 -o ${dir}/ignore/cluster_output 
```

while grf_fslcluster.sh reads as: 
```
#!/usr/bin/env bash

###################################################################
###################################################################

###################################################################
# Combine all text file output
###################################################################

###################################################################
# Usage function
###################################################################
Usage(){
  echo ""; echo ""; echo ""
  echo "Usage: `basename $0`  grf_fslcluster.sh -i zstat -m mask -t threshold -o output"
  echo ""
  echo "Compulsory arguments:"
  echo "  -i : zstats: compulsory"
  echo "  -m: mask"
  echo "  -o : Output file name"
  echo "       "
  exit 2
}

###################################################################
# Parse arguments
###################################################################
while getopts "i:t:m:o:h" OPTION ; do
  case ${OPTION} in
    i)
      zstat=${OPTARG}
      ;;
    t)
      thresh=${OPTARG}
      ;;
    m)
      mask=${OPTARG}
      ;;
    o)
      outdir=${OPTARG}
      ;;
    h)
      Usage
      ;;
    *)
      Usage
      ;;
    esac
done

###################################################################
# Ensure that all compulsory arguments have been defined
###################################################################
[[ -z ${outdir} ]] && Usage
[[ -z ${zstat} ]] && Usage
[[ -z ${mask} ]] && Usage

###################################################################
# Now run through each file that we find and append it to the output file
###################################################################
 
if [[ -z ${thresh} ]]; then 
   thresh=2.3
   echo "voxel threshold is 2.3 (default)"
fi 

echo " find d and v " 
dv=$(smoothest -z ${zstat} -m ${mask})

id0=$(echo $dv |cut -d' ' -f2)
id1=$(echo $dv |cut -d' ' -f4)
echo " the dlh is ${id0}"
echo "                  "
echo " the number of volume: ${id1}"
echo $thresh
echo $dv
mkdir -p ${outdir}/cluster_Z${thresh}

cluster -i ${zstat} -d ${id0} --volume=${id1} -t ${thresh} -p 0.05 \
   -o  ${outdir}/cluster_Z${thresh}/cluster_Z${thresh} >  \
    ${outdir}/cluster_Z${thresh}/cluster_Z${thresh}.csv

echo "done"
```

The output of cluster masks is at: `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09`. 

Numbers obtained from the CSV slightly different than before, but ultimately a close replication for the second cluster(likely due to changes in software version): 
| Cluster Index | Voxels | P | -log10(P) | MAX | MAX X (vox) | MAX Y (vox) | MAX Z (vox) | COG X (vox) | COG Y (vox) | COG Z (vox) |
| ---- | ---- | ---- | ------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|1|	12|	0.000103|	3.99|	3.54|	30|	44|	30|	30.9|	43.8|	30.3|

The first cluster was not replicated, but this had never survived sensitivity analyses. It was inconsistently generated in different locations, and sometimes not at all - we excluded it from final analysis. 


However, the output image of clusters produced was an .img and .hdr which needed to be turned into a nifti format, I used: 

 ```fslchfiletype NIFTI_GZ cluster_Z3.09.img cluster_Z3.09.nii```

Notably, the images were very sparse and hard to see - if viewed in mrpeek, appeared to be empty images with "NaN" scales - in Mango, only one cluster was visible. Ultimately downloaded MRICRON to view images locally - this worked much better. I used a template of `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/pnc_template_brain_2mm.nii.gz`



### 3. Seed-based correlation 
A mask was generated from the cluster_z3_09.nii, using fslmath at `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1` :

(See [https://mandymejia.com/fsl-maths-commands/](https://mandymejia.com/fsl-maths-commands/))

```
fslmaths cluster_Z3.09.nii.gz -thr 1 -uthr 2 mask1/mask1.nii.gz #CLUSTER1
```

Mask generated wsa again in .hdr and .img format, so I used fslchfiletype to turn them into niftis.

The mask was upsampled from 4mm to 2mm and was used as a seed for seed-based correlation.

This upsampling was done using the pnc_2mm template, found at `/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/pnc_template_brain_2mm.nii.gz`

The code used for this was: 
```
3dresample -master pnc_template_brain_2mm.nii.gz -input mask1.nii.gz -prefix mask1_2mm.nii.gz
```

The path to the seed is: 
`~/cluster_Z3.09/mask1/mask1_2mm.nii.gz `

The seed-based correlation was computed with the following script, and the `xcpengine.simg` file under `/cbica/projects/pncitc/ignore`:

```
#!/bin/bash #DEFINITELY QSUB THIS
cd /cbica/projects/pncitc/ignore
XCPEDIR=xcpEngine
seedpoint1=/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mm.nii.gz

bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv
image=/cbica/projects/pncitc/subjectData/rest293/
outputdir=/cbica/projects/pncitc/ignore/seedcorrmaps

mkdir -p ${outputdir}

cat $bblid | while IFS="," read -r a b ; 

do
     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/ignore /cbica/projects/pncitc/ignore/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint1 -o $outputdir -p ${a},${b} -k 6 -n mask1\
     #rm -rf $outputdir/seed/mask1/${a}_${b}_connectivity_mask1_seed.nii.gz

done
```

*NOTE: For this script, the output displays in the error log for some reason. Also note, the demographic csvs say "scandid", not "scanid"*

### 4. Linear regression with FSL `flameo` 

Flameo regression computation requires `design`,`contrast` and `group` text files. The script `scripts/makeflameodesig.R`:

```
# script to make the design matrices for flameo

library(pracma)
demogr=read.csv('~/n293_demographics.csv') 
#logk+relMeanRMSmotion+age+sex 
desigmatlogkonly=cbind(rep(1,307),demogr$logk,demogr$sex,demogr$relMeanRMSmotion,demogr$age)

grp=ones(293,1) # only one group

contrast4=zeros(5,5); 
diag(contrast4)=1; 


write.table(desigmatlogkonly,'~/desigmatlogkonly.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(contrast4,'~/contrast4.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(grp,'~/grp.txt',sep=' ',quote = FALSE,row.names = FALSE,col.names = FALSE)
```

 was used to generate these, but the package was out of date in R on CBICA, so I ran it locally and copied the needed files over into `/cbica/projects/pncitc/ignore/seedcorrmaps`


I also converted the designlogkmat from .txt to .mat and did the same for grp.txt, using: 

```
Text2Vest desigmatlogkonly.txt desigmatlogkonly.mat
Text2Vest grp.txt grp.grp
Text2Vest contrast4.txt contrast4.con
```
I also created the flameo csvs pointing to the mask1Z_sm6.nii.gz niftis in the same directory,`mask1.csv`  (`/cbica/projects/pncitc/ignore/seedcorrmaps`). 

The flameo linear regression was computed with this script:

```
#!/bin/bash
#$ -l h_vmem=200G
#$ -l tmpfree=200G # this will throw an error about ": $PATH does not agree with $PATH_modshare counter.", but it shouldn't stop the process. 
bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv #FIRST HALF RUNS QUICKLY, WOULD QSUB BECAUSE THE SECOND PART TAKES A BIT LONGER
imagedir=/cbica/projects/pncitc/ignore/seedcorrmaps/seed
scriptdir=/cbica/projects/pncitc/ignore/seedcorrmaps
outputdir=/cbica/projects/pncitc/ignore/seedcorrmaps/seed
demogdir=/cbica/projects/pncitc/ignore/seedcorrmaps


imagelist1=$scriptdir/mask1.csv

rm -rf $imagelist1


cat $bblid | while IFS="," read -r a b ; 

do 
     img1=$(ls -f $imagedir/mask1/${a}_${b}_connectivity_mask1Z_sm6.nii.gz)
     
     echo $img1 >> $imagelist1

done 


mask=/cbica/projects/pncitc/subjectData/PNCgrey2mm.nii.gz

fslmerge -t ${outputdir}/4Dcopeseed1.nii.gz $(cat $imagelist1)

flameo --copefile=${outputdir}/4Dcopeseed1.nii.gz   --mask=${mask}   --dm=${demogdir}/desigmatlogkonly.mat  --tc=${demogdir}/contrast4.con  --cs=${demogdir}/grp.grp --runmode=flame1 --ld=$outputdir/mask1/logk #SECOND PART, WHICH TAKES LONGER

```



The outputs of flameo regression: 

`/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk`


In the directory, there are zvalues: 

  `zstat1 : average `

  `zstat2 : logk `

  `zstat3 : sex`

  `zstat4 : motion`

  `zstat5 : age`


The zstats  were FDR corrected with this script `scripts/flameoutputfdrcorrection` in R on CBICA. For this, 

```
rm(list = ls())
library(RNifti, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
# for mask1 or seed1
for (i in 1:5 ){  
    mask=readNifti('mask1/logk/mask.nii.gz')
    z1=readNifti(paste0('mask1/logk/zstat',i,'.nii.gz'))
    Z=z1[mask==1]
    p <- 2*pnorm((-abs(Z)))
    p1=p.adjust(p, method = 'fdr')
    zvals = qnorm(1 - (p1/2)) 
    zvals[zvals==Inf]=10
    Z[Z>0]=1; Z[Z<0]=-1
    zm=mask
    zm[mask==1]=zvals*Z
    writeNifti(zm,paste0('mask1/logk/zfdr',i,'.nii.gz'),template = mask)
}


```


The outputs are located here:

`/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1`

FDR corrected z-values. 

  `zfdr1 : average `

  `zfdr2 : logk `

  `zfdr3 : sex`

  `zfdr4 : motion`

  `zfdr5 : age`
  


### 5. Vizualisation of Results - iPython in CBICA
All computations were done in PNC template. For vizualisation, all the nifti files  were tranformed to MNI before as below: 

```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]
flame1dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/'
zstats=['zfdr1','zfdr2']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz'
    at.run()

    
clusterdirectory = '/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09'
zstats=['/mask1/mask1_2mm','/mask1/mask1']
for i in range(len(zstats)):
    at.inputs.input_image = clusterdirectory + zstats[i]+'.nii.gz'
    at.inputs.output_image = clusterdirectory + zstats[i]+'MNI.nii.gz'
    at.run()
```

 a. for clusters and mean of seed-based correlation: 
_Note: had to use `flchfiletype` on the `copeseed` images before running the script, and move the `.hdr` and `.img` files out of the directory/ remove them altogether - having the nifti and img in the same directory can cause an error._
```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms
from nipype.interfaces.fsl import MultiImageMaths,maths, MeanImage

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

output_image = '/cbica/projects/pncitc/ignore/cluster_output/cluster_Z3.09/mask1/mask1_2mmMNI.nii.gz'
#put it on surface 
img1=img.load_img(output_image)
v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0.001,vmax=.5,title='clusters',cmap='RdYlBu_r') 

# plot of mean of seed-based correlation 

# average of all subject 
seedbasedir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/'
corrtm=['4Dcopeseed1'] # make sure to change to nifti and remove .img
viewim=[]
meanimage=MeanImage()
for i in range(len(corrtm)):
    meanimage.inputs.in_file=seedbasedir + corrtm[i]+ '.nii.gz'
    meanimage.inputs.dimension='T'
    meanimage.inputs.out_file=seedbasedir +corrtm[i] + 'mean.nii.gz' 
    meanimage.run()
    os.system('fslchfiletype NIFTI_GZ 4Dcopeseed1mean.img 4Dcopeseed1mean.nii.gz')
    os.system('rm -rf 4Dcopeseed1mean.hdr  4Dcopeseed1mean.img')
    at.inputs.input_image = meanimage.inputs.out_file
    at.inputs.output_image = seedbasedir +corrtm[i] + 'meanMNI.nii.gz'
    at.run()
    img1=img.load_img(at.inputs.output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0.1,vmax=0.5,title='mean of seed-based correlation' + ': mask' +str(i + 1 ),cmap='RdYlBu_r') 
    viewim.append(v)
    
ii = 0    
for x in viewim:
  ii += 1
  x.save_as_html("meanseedbasedcorr" + str(ii) + ".html")
 

```

  b. for mask1 (ran in iPython on CBICA): `notebook/flameomask1.ipynb`
  

  ```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")


import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

flame1dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/'
zstats=['zfdr1','zfdr2']
label=['mean','logk']
viewim=[]
for i in range(len(zstats)):
    output_image = flame1dir + zstats[i]+'MNI.nii.gz' 
    img1=img.load_img(output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0,vmax=5,title='zstat :'+label[i],cmap='RdYlBu_r') 
    viewim.append(v)

ii = 0
for x in viewim:
    ii+=1
    x.save_as_html("cluster"+str(ii)+".html")
  ```


Images I generated were saved in `/cbica/projects/pncitc/ignore/seedcorrmaps/seed` in the .html format. 

### 6. Regional plot of significant regions of logk
 
 I did this step locally as the R libraries on CBICA were very out of date. For N293, I redid this entire script as I found an error in the masks - it seems positive and negative masks were switched, and the negative masks were multiplied by -1. 
 
 Additionally, I pulled out all non-zero values for corrdata rather than values equal to 1, as was done originally. 
 
 Finally, I generated the insets seen in the manuscript (these were the niftis written out in the script below, I then projected them to the surface using similar code as in `/notebook/flameomask1.ipynb` on this Github.)
 
 These surface projections are in .html within `.../ignore/Vis`. My code for this part is as below:

```
library(RNifti, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(pracma, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(ggplot2, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(visreg, lib.loc = '/cbica/projects/pncitc/mehtareplicate')
library(Matrix, lib.loc = '/cbica/projects/pncitc/mehtareplicate') # and really any other packages that give you issues - as this can happen
# make the masks
mask1=readNifti('/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/zfdr2.nii.gz') # the fdr corrected flameo output for logk 

#get the  postive masks
p_m1=mask1; p_m1[p_m1<1.64]=0

#get the negative masks
n_m1=(-1)*mask1;  n_m1[n_m1<1.64]=0; n_m1 = (-1)*n_m1


writeNifti(p_m1, '/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/p_m1.nii.gz', template = NULL, datatype = "auto", version = 1)
writeNifti(n_m1, '/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/n_m1.nii.gz', template = NULL, datatype = "auto", version = 1)

b=read.csv('/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv',header=FALSE)

#make table 

corrdata=zeros(293,4)

for (i in 1:293) {
  img1=readNifti(paste0('/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/',b[i,1],'_',b[i,2],'_connectivity_mask1Z_sm6.nii.gz')) # flameo output
  datap1=img1[p_m1!=0]
  datam1=img1[n_m1!=0]
  corrdata[i,]=c(b[i,1],b[i,2],mean(datap1),mean(datam1))
}

colnames(corrdata)=c('bblid','scanid','mask1pos','mask1neg')

write.csv(corrdata,'/cbica/projects/pncitc/demographics/n293_meanseedcorr.csv',quote = FALSE,row.names = FALSE)

# merge CSV
x = read.csv('/cbica/projects/pncitc/demographics/n293_meanseedcorr.csv')
y = read.csv('/cbica/projects/pncitc/demographics/n307_demographics.csv') # demographics are right, when merged the n307 will become n293
z = read.csv('/cbica/projects/pncitc/demographics/n2416_RestQAData_20170714.csv')

final1=merge(x,y, by=c('bblid','scanid')) # merge by Ids  
final2=merge(final1,z, by=c('bblid','scanid')) # merge by Ids 
write.csv(final2,'n293_data.csv',quote = FALSE,row.names = FALSE)

# write as .rds
saveRDS(final2, file = "/cbica/projects/pncitc/demographics/my_data.RDS") 

#start plotting
ddata=readRDS('/cbica/projects/pncitc/demographics/my_data.RDS')
poscluster1mask_nologk=lm(mask1pos~age+sex+relMeanRMSmotion,data=ddata)
poscluster1mask=lm(mask1pos~logk+age+sex+relMeanRMSmotion,data=ddata)

ylab<-"Correlation (z(r))"

ddata$poscluster1resid<-poscluster1mask_nologk$residuals+mean(ddata$mask1pos)
ggplot(ddata,aes(x=logk,y=poscluster1resid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))
ggsave('/cbica/projects/pncitc/ignore/cluster1pos.png')

ddata$negcluster1resid<-negcluster1mask_nologk$residuals+mean(ddata$mask1neg)

ggplot(ddata,aes(x=logk,y=negcluster1resid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ggsave('/cbica/projects/pncitc/ignore/cluster1neg.png')
```
Insets were finally visualized on CBICA via iPython using: 

```
# import all the requirements and hide warnings
import warnings
warnings.filterwarnings("ignore")
import os

import nilearn.plotting as plott
import nilearn.image as img
from nilearn import datasets,surface
import matplotlib.pyplot as plt
from nipype.interfaces.ants import ApplyTransforms

big_fsaverage = datasets.fetch_surf_fsaverage('fsaverage') # for viz 

#registration paramteters

ref='/cbica/projects/pncitc/subjectData/PNC_transforms/MNI152_T1_2mm_brain.nii.gz'
transform1='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_0Warp.nii.gz'
transform2='/cbica/projects/pncitc/subjectData/PNC_transforms/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

dir='/cbica/projects/pncitc/ignore/seedcorrmaps/seed/mask1/logk/zfdrmask1/'
masks=['n_m1', 'p_m1']
mask = masks
for i in range(len(masks)):
    at.inputs.input_image = dir + masks[i]+'.nii.gz'
    at.inputs.output_image = dir + masks[i]+'MNI.nii.gz'
    img2 = mask[i]+'.img'
    nii = mask[i]+'.nii.gz'
    os.system('fslchfiletype NIFTI_GZ ' + img2 + ' ' + nii)
    os.system('rm -rf *img')
    at.run()

viewim=[]
label = ['neg','pos']

for i in range(len(masks)):
  output_image = dir + masks[i]+'MNI.nii.gz'
  img1=img.load_img(output_image)
  v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0,vmax=5,title='zstat :'+label[i],cmap='RdYlBu_r') 
  viewim.append(v)

ii = 0
for x in viewim:
  ii+=1
  x.save_as_html('inset'+str(ii)+'.html')

```

Finally, for the cluster image, I locally ran this on the mask1 in MNI: 

```
from nilearn import plotting
import nibabel as nib

stat_map_img = '/Users/kahinim/Desktop/mask1MNI.nii.gz'
plotting.plot_glass_brain(stat_map_img, output_file='/Users/kahinim/Desktop/test.png', display_mode='ortho', colorbar=False, figure=None, axes=None, title=None, threshold='auto', annotate=True, black_bg=False, cmap=None, alpha=0.7, vmin=None, vmax=None, plot_abs=True, symmetric_cbar='auto', resampling_interpolation='continuous')
```
------------------------------------------------------------------------------------------------------------------------------------------------------


