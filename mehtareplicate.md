---
title: "PNC-ITC Replication"
authors: Kahini
output: html_document
---
#### Copied and modified from pncitc.md

Sample replication: n307 did not exclude those with the `health_exclude` criteria. Analyses were re-run on n293, using Pehlivanova et al's n427 sample and then running restQA exclusions on them (all this information was from Pehlivanova .csvs). The code and .csvs for this are available at `cbica/projects/pncitc/samplerecreation`. All subsequent analyses were run in `cbica/projects/pncitc/mehtareplicaten293`, using the same steps as below in folder `mehtareplicaten293` - additionally, any new .csvs should be pointed to in the scripts in that directory. Any departure from the protocol for n307 will be noted here.  I also moved the bblid_scanid .csv to demographics, and created a folder within subjectData called rest293 for the n = 293 replication.
**Results were similar in N293 for the second cluster, things changed for the first cluster. Visualizations are available in the .html format within mehtareplicaten293**

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
write.csv(ids,'n294_blbid_scanid.csv',row.names = FALSE,quote = FALSE)
pncit1$age=pncit1$ageAtScan/12
write.csv(pncit1,'n294_demographics.csv',row.names = FALSE,quote = FALSE)
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
mdmrouput=/cbica/projects/pncitc/mehtareplicaten293/cwas293 #output directory
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc
imagelist=/cbica/projects/pncitc/subjectData/imageinput_rest3.csv #list of image in nifti # HAD TO RE-GENERATE THIS LIST AS THE FILE PATHS HAD CHANGED
rm  -rf $mdmrouput # remove previous run if exist 
metric=pearson # pearson correlation 
# compute distance matrix
singularity exec -e -B /cbica/projects/pncitc $singimage $scriptdir/Rscript $scriptdir/connectir_subdist.R $mdmrouput --infuncs1=$imagelist --ztransform --automask1  -c 3 -t 4 --brainmask1=$brainmask --method="$metric" --bg=$bgim  --overwrite --memlimit=200

```

The output of distance matrix: `/cbica/projects/pncitc/mehtareplicaten293/cwas293`
   
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
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/mehtareplicaten293/cwas293 -f 'logk+relMeanRMSmotion+sex+age' -m /cbica/projects/pncitc/samplerecreation/n293_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```

Memory and formatting were the main problems with these scripts not running well - the numbers/format left in were what worked for me. The output is at: `/cbica/projects/pncitc/mehtareplicate/cwas293/logk_motion_sex_age`

### 2. Significant clusters from mdmr
The cluster analysis was computed  with the script `scripts/grf_fslcluster.sh`, written based on  [FSL cluster analysis](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster) with  Gaussian Random Field (GRF) theory

The script `scripts/cluster.sh` called grf_fslcluster.sh, as listed in this repo under scripts, with `z=3.09`
Two clusters were obtained: one at the at the frontal region the order at the TPJ. 

cluster.sh reads as below:
```
 #!/bin/bash # NO NEED TO QSUB
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/mehtareplicaten293/cwas293/logk_motion_sex_age/zstats_logk.nii.gz  -m ${dir}/mehtareplicaten293/cwas293/mask.nii.gz -t 3.09 -o ${dir}/mehtareplicaten293/cluster_output 
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

The output of cluster masks is at: `/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09`. 

Numbers obtained from the CSV slightly different than before, but ultimately a close replication (likely due to changes in software version): 
| Cluster Index | Voxels | P | -log10(P) | MAX | MAX X (vox) | MAX Y (vox) | MAX Z (vox) | COG X (vox) | COG Y (vox) | COG Z (vox) |
| ---- | ---- | ---- | ------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 2 | 11 | 9.27E-05 | 4.03 | 3.54 | 30 | 44 | 30 | 30.9 | 43.9 | 30.4 |
| 1 | 5 | 0.0228 | 1.64	 | 3.54 | 13 | 28 | 24 | 13.4 | 28.6 | 24.2 |

However, the output image of clusters produced was an .img and .hdr which needed to be turned into a nifti format, I used: 

 ```fslchfiletype NIFTI_GZ cluster_Z3.09.img cluster_Z3.09.nii```

Notably, the images were very sparse and hard to see - if viewed in mrpeek, appeared to be empty images with "NaN" scales - in Mango, only one cluster was visible. Ultimately downloaded MRICRON to view images locally - this worked much better. 



### 3. Seed-based correlation 
Two different masks were generated from the cluster_z3_09.nii, using fslmath at `/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09/mask*` :

(See [https://mandymejia.com/fsl-maths-commands/](https://mandymejia.com/fsl-maths-commands/))

```
fslmaths cluster_Z3.09.nii.gz -thr 1 -uthr 2 mask1/mask1.nii.gz #CLUSTER1
fslmaths cluster_Z3.09.nii.gz -thr 2 -uthr 2 mask2/mask2.nii.gz #CLUSTER2
```

Masks generated were again in .hdr and .img format, so I used fslchfiletype to turn them into niftis.

The two masks were upsampled from 4mm to 2mm and were used as seeds for seed-based correlation.

This upsampling was done using the pnc_2mm template, found at `/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09/pnc_template_brain_2mm.nii.gz`

The code used for this was: 
```
3dresample -master pnc_template_brain_2mm.nii.gz -input mask1.nii.gz -prefix mask1_2mm.nii.gz
3dresample -master pnc_template_brain_2mm.nii.gz -input mask2.nii.gz -prefix mask2_2mm.nii.gz
```

The path to the two seeds is: 
`~/cluster_Z3.09/mask1/mask1_2mm.nii.gz `
`~/cluster_Z3.09/mask2/mask2_2mm.nii.gz` 

The seed-based correlation was computed with the following script, and the `xcpengine.simg` file under `/cbica/projects/pncitc/mehtareplicate`:

```
#!/bin/bash #DEFINITELY QSUB THIS
cd /cbica/projects/pncitc/mehtareplicaten293
XCPEDIR=xcpEngine
seedpoint1=/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09/mask1/mask1_2mm.nii.gz
seedpoint2=/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09/mask2/mask2_2mm.nii.gz

bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv
image=/cbica/projects/pncitc/subjectData/rest293/
outputdir=/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps

mkdir -p ${outputdir}

cat $bblid | while IFS="," read -r a b ; 

do
     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/mehtareplicaten293 /cbica/projects/pncitc/mehtareplicaten293/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint1 -o $outputdir -p ${a},${b} -k 6 -n mask1\
     #rm -rf $outputdir/seed/mask1/${a}_${b}_connectivity_mask1_seed.nii.gz


     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/mehtareplicaten293 /cbica/projects/pncitc/mehtareplicaten293/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint2 -o $outputdir -p ${a},${b} -k 6 -n mask2\
     #rm -rf $outputdir/seed/mask2/${a}_${b}_connectivity_mask2_seed.nii.gz
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

 was used to generate these, but the package was out of date in R on CBICA, so I ran it locally and copied the needed files over into `/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps`


I also converted the designlogkmat from .txt to .mat and did the same for grp.txt, using: 

```
Text2Vest desigmatlogkonly.txt desigmatlogkonly.mat
Text2Vest grp.txt grp.grp
Text2Vest contrast4.txt contrast4.con
```
I also created the flameo csvs pointing to the mask1Z_sm6.nii.gz niftis in the same directory, naming the csvs `mask1.csv` and `mask2.csv` respectively (`/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps`). 

The flameo linear regression was computed with this script:

```
#!/bin/bash
#$ -l h_vmem=200G
#$ -l tmpfree=200G # this will throw an error about ": $PATH does not agree with $PATH_modshare counter.", but it shouldn't stop the process. 
bblid=/cbica/projects/pncitc/demographics/n293_bblid_scandid.csv #FIRST HALF RUNS QUICKLY, WOULD QSUB BECAUSE THE SECOND PART TAKES A BIT LONGER
imagedir=/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed
scriptdir=/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps
outputdir=/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed
demogdir=/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps


imagelist1=$scriptdir/mask1.csv
imagelist2=$scriptdir/mask2.csv


rm -rf $imagelist1
rm -rf $imagelist2


cat $bblid | while IFS="," read -r a b ; 

do 
     img1=$(ls -f $imagedir/mask1/${a}_${b}_connectivity_mask1Z_sm6.nii.gz)
     img2=$(ls -f $imagedir/mask2/${a}_${b}_connectivity_mask2Z_sm6.nii.gz)
     
     echo $img1 >> $imagelist1
     echo $img2 >> $imagelist2

done 


mask=/cbica/projects/pncitc/subjectData/PNCgrey2mm.nii.gz

fslmerge -t ${outputdir}/4Dcopeseed1.nii.gz $(cat $imagelist1)
fslmerge -t ${outputdir}/4Dcopeseed2.nii.gz $(cat $imagelist2)

flameo --copefile=${outputdir}/4Dcopeseed1.nii.gz   --mask=${mask}   --dm=${demogdir}/desigmatlogkonly.mat  --tc=${demogdir}/contrast4.con  --cs=${demogdir}/grp.grp --runmode=flame1 --ld=$outputdir/mask1/logk #SECOND PART, WHICH TAKES LONGER

flameo --copefile=${outputdir}/4Dcopeseed2.nii.gz   --mask=${mask}   --dm=${demogdir}/desigmatlogkonly.mat  --tc=${demogdir}/contrast4.con  --cs=${demogdir}/grp.grp --runmode=flame1 --ld=$outputdir/mask2/logk
```



The outputs of flameo regression: 

`/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask1/logk`
`/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask2/logk`

In each directory, there are zvalues: 

  `zstat1 : average `

  `zstat2 : logk `

  `zstat3 : sex`

  `zstat4 : motion`

  `zstat5 : age`


The zstats  were FDR corrected with this script `scripts/flameoutputfdrcorrection`. For this, 
1. I opened 'R' in CBICA
2. Changed mask.img in each of the `logk` files to .nii as before, using `fslchtype`
3. Did:


```
setwd('~/seedcorrmaps')
.libPaths("/gpfs/fs001/cbica/home/mehtaka/R/x86_64-conda_cos6-linux-gnu-library/3.6")
```

Then chose to install [RNifti](https://rdrr.io/cran/RNifti/f/README.md) in a personal library, using mirror #78. Steps involved :
1. install.packages("remotes")
2. install.packages("devtools")
3. install_github("jonclayden/RNifti")

Then I ran the rest of the script:

```
rm(list = ls())
library(RNifti)
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

# for mask2 or seed2
for (i in 1:5 ){  
    mask=readNifti('mask2/logk/mask.nii.gz')
    z1=readNifti(paste0('mask2/logk/zstat',i,'.nii.gz'))
    Z=z1[mask==1]
    p <- 2*pnorm((-abs(Z)))
    p1=p.adjust(p, method = 'fdr')
    zvals = qnorm(1 - (p1/2)) 
    zvals[zvals==Inf]=10
    Z[Z>0]=1; Z[Z<0]=-1
    zm=mask
    zm[mask==1]=zvals*Z
    writeNifti(zm,paste0('mask2/logk/zfdr',i,'.nii.gz'),template = mask)
}

```


The outputs are located here:

`/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask1/logk/zfdrmask1`
`/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask2/logk/zfdrmask2`

FDR corrected z-values. 

  `zfdr1 : average `

  `zfdr2 : logk `

  `zfdr3 : sex`

  `zfdr4 : motion`

  `zfdr5 : age`
  

**Changes: ran locally for n293 due to library difficulties on CBICA**

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
flame1dir='/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask1/logk/zfdrmask1/'
zstats=['zfdr1','zfdr2']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz'
    at.run()

flame2dir='/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask2/logk/zfdrmask2/'
zstats=['zfdr1','zfdr2']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame2dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame2dir + zstats[i]+'MNI.nii.gz'
    at.run()
    
cluster_directory = '/cbica/projects/pncitc/mehtareplicaten293/cluster_output/cluster_Z3.09'
zstats=['/mask1/mask1_2mm','/mask2/mask2_2mm']
for i in range(len(zstats)):
    at.inputs.input_image = clusterdirectory + zstats[i]+'.nii.gz'
    at.inputs.output_image = clusterdirectory + zstats[i]+'MNI.nii.gz'
    at.run()
```

 a. for clusters and mean of seed-based correlation: 
_Note: had to use `flchfiletype` on the `copeseed` images, `mask1and2_2mm.nii`, `mean copeseed` images, and move the `.hdr` and `.img` files out of the directory/ remove them altogether - having the nifti and img in the same directory can cause an error._
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

ref='~/MNI152_T1_2mm_brain.nii.gz'
transform1='~/PNC-MNI_0Warp.nii.gz'
transform2='~/PNC-MNI_1Affine.mat'
at = ApplyTransforms()
at.inputs.dimension = 3
at.inputs.reference_image = ref
at.inputs.interpolation = 'NearestNeighbor'
at.inputs.default_value = 0
at.inputs.transforms = [transform1, transform2]
at.inputs.invert_transform_flags = [False, False]

# combine cluster mask1 and mask2

#multiply second cluster by 2 to differentiate from cluster 1 and add both cluster together 
mask2=MultiImageMaths()
mask2.inputs.in_file='~/mask2_2mm.nii.gz'
mask2.inputs.op_string='-mul -1 -add %s'
mask2.inputs.out_file='~/mask1and2_2mm.nii.gz' # make sure to change to nifti and remove .img
mask2.inputs.operand_files='/~/mask1_2mm.nii.gz'
mask2.cmdline
mask2.run()

#register it to MNI template 
at.inputs.input_image = mask2.inputs.out_file
at.inputs.output_image = '~/mask1and2_2mmMNI.nii.gz'
at.run()

#put it on surface 
img1=img.load_img(at.inputs.output_image)
v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0.001,vmax=.5,title='clusters',cmap='RdYlBu_r') 

# plot of mean of seed-based correlation of seed1 and seed2

# average of all subject 
seedbasedir='/Users/adebimpe/Box/projects/ITC2/output/seedcorrmaps/'
corrtm=['4Dcopeseed1','4Dcopeseed2'] # make sure to change to nifti and remove .img
viewim=[]
meanimage=MeanImage()
for i in range(len(corrtm)):
    meanimage.inputs.in_file=seedbasedir + corrtm[i]+ '.nii.gz'
    meanimage.inputs.dimension='T'
    meanimage.inputs.out_file=seedbasedir +corrtm[i] + 'mean.nii.gz' # make sure to change to nifti and remove .img
    meanimage.run()
    at.inputs.input_image = meanimage.inputs.out_file
    at.inputs.output_image = seedbasedir +corrtm[i] + 'meanMNI.nii.gz'
    at.run()
    img1=img.load_img(at.inputs.output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0.1,vmax=0.5,title='mean of seed-based correlation' + ': mask' +str(i + 1 ),cmap='RdYlBu_r') 
    viewim.append(v)
  

ii = 4
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
 

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

flame1dir='~/mask1/logk/'
zstats=['zfdr1','zfdr2']
label=['mean','logk']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz' # this should already be done! If not remember to change to .nii from .img, and remove .img
    at.run()
    img1=img.load_img(at.inputs.output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0,vmax=5,title='zstat :'+label[i],cmap='RdYlBu_r') 
    viewim.append(v)

  ii = 0
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
  ```

  c. for mask2 (iPython) : `notebook/flameomask2.ipynb`
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

flame1dir='~/mask2/logk/'
zstats=['zfdr1','zfdr2']
label=['mean','logk']
viewim=[]
for i in range(len(zstats)):
    at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz' # this should already be done! If not remember to change to .nii from .img, and remove .img
    at.run()
    img1=img.load_img(at.inputs.output_image)
    v= plott.view_img_on_surf(img1, surf_mesh='fsaverage',threshold=0,vmax=5,title='zstat :'+label[i],cmap='RdYlBu_r') 
    viewim.append(v)

  ii = 2
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
  ```

Images I generated were saved in `/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed` in the .html format. 

### 6. Regional plot of significant regions of logk
 
 I did this step locally as the R libraries on CBICA were very out of date. For N293, I redid this entire script as I found an error in the masks - it seems positive and negative masks were switched, and the negative masks were multiplied by -1. 
 
 Additionally, I pulled out all non-zero values for corrdata rather than values equal to 1, as was done originally. 
 
 Finally, I generated the insets seen in the manuscript (these were the niftis written out in the script below, I then projected them to the surface using similar code as in `/notebook/flameomask1.ipynb` on this Github.)
 
 These surface projections are in .html within `.../mehtareplicaten293/Vis`. My code for this part is as below, and all the files are in the `dropbox` under final_visualization**

```
library(RNifti)
library(pracma)

setwd('/Users/kahinim/Desktop')

# make the masks
mask1=readNifti('1Zfiles/zfdr2.nii.gz') # the fdr corrected flameo output for logk 
mask2=readNifti('2Zfiles/zfdr2.nii.gz')

#get the  postive masks
p_m1=mask1; p_m1[p_m1<1.64]=0
p_m2=mask2; p_m2[p_m2<1.64]=0

#get the negative masks
n_m1=(-1)*mask1;  n_m1[n_m1<1.64]=0; n_m1 = (-1)*n_m1
n_m2=(-1)*mask2;  n_m2[n_m2<1.64]=0; n_m2 = (-1)*n_m2

writeNifti(p_m1, 'p_m1.nii.gz', template = NULL, datatype = "auto", version = 1)
writeNifti(n_m1, 'n_m1.nii.gz', template = NULL, datatype = "auto", version = 1)
writeNifti(n_m2, 'n_m2.nii.gz', template = NULL, datatype = "auto", version = 1)
writeNifti(p_m2, 'p_m2.nii.gz', template = NULL, datatype = "auto", version = 1)

b=read.csv('n293_bblid_scanid.csv',header=FALSE)

#make table 

corrdata=zeros(293,6)

for (i in 1:293) {
  img1=readNifti(paste0('1Zfiles/',b[i,1],'_',b[i,2],'_connectivity_mask1Z_sm6.nii.gz')) # flameo output
  img2=readNifti(paste0('2Zfiles/',b[i,1],'_',b[i,2],'_connectivity_mask2Z_sm6.nii.gz'))
  datap1=img1[p_m1!=0]; datap2=img2[p_m2!=0]
  datam1=img1[n_m1!=0]; datam2=img2[n_m2!=0]
  corrdata[i,]=c(b[i,1],b[i,2],mean(datap1),mean(datam1),mean(datap2),mean(datam2))
}

colnames(corrdata)=c('bblid','scanid','mask1pos','mask1neg','mask2pos','mask2neg')

write.csv(corrdata,'n293_meanseedcorr.csv',quote = FALSE,row.names = FALSE)

# merge CSV
x = read.csv('n293_meanseedcorr.csv')
y = read.csv('n307_demographics.csv') # demographics are right, when merged the n307 will become n293
z = read.csv('n2416_RestQAData_20170714.csv')

final1=merge(x,y, by=c('bblid','scanid')) # merge by Ids  
final2=merge(final1,z, by=c('bblid','scanid')) # merge by Ids 
write.csv(final2,'n293_data.csv',quote = FALSE,row.names = FALSE)

# write as .rds
saveRDS(final2, file = "my_data.RDS") 

#start plotting
ddata=readRDS('/Users/kahinim/Desktop/my_data.rds')
library(visreg);
poscluster1mask_nologk=lm(mask1pos~age+sex+relMeanRMSmotion,data=ddata)
poscluster1mask=lm(mask1pos~logk+age+sex+relMeanRMSmotion,data=ddata)

#summary(postpjmask)
#svg("/Users/kahinim/Desktop/mask1pos.svg", width = 8, height = 8)
imageplot<-visreg(poscluster1mask, "logk", 
                  xlab="Log K", 
                  ylab="Correlation",line=list(col="red",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "red"),
                  fill=list(col=adjustcolor("red", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,ylim=c(-0.5,0.5),xlim=c(-9,-1) )

negcluster1mask=lm(mask1neg~age+sex+relMeanRMSmotion+logk,data=ddata)
negcluster1mask_nologk=lm(mask1neg~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
#svg("/Users/kahinim/Desktop/mask1neg.svg", width = 8, height = 8)
imageplot<-visreg(negcluster1mask, "logk", 
                  xlab="Log K", 
                  ylab="Correlation k",line=list(col="blue",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "blue"),
                  fill=list(col=adjustcolor("blue", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )

poscluster2mask=lm(mask2pos~age+sex+relMeanRMSmotion+logk,data=ddata)
poscluster2mask_nologk=lm(mask2pos~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
#svg("/Users/kahinim/Desktop/mask2pos.svg", width = 8, height = 8)
imageplot<-visreg(poscluster2mask, "logk", 
                  xlab="Log K", 
                  ylab="Correlation",line=list(col="red",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "red"),fill=list(col=adjustcolor("red", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )

negcluster2mask=lm(mask2neg~age+sex+relMeanRMSmotion+logk,data=ddata)
negcluster2mask_nologk=lm(mask2neg~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
#svg("/Users/kahinim/Desktop/mask2neg.svg", width = 8, height = 8)
imageplot<-visreg(negcluster2mask, "logk", 
                  xlab="Log K", 
                  ylab="Correlation",line=list(col="blue",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "blue"), fill=list(col=adjustcolor("blue", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )


library(ggplot2)

ylab<-"Correlation (z(r))"

ddata$poscluster1resid<-poscluster1mask_nologk$residuals+mean(ddata$mask1pos)
ggplot(ddata,aes(x=logk,y=poscluster1resid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))
ggsave('/Users/kahinim/Desktop/cluster1pos.png')

ddata$negcluster1resid<-negcluster1mask_nologk$residuals+mean(ddata$mask1neg)

ggplot(ddata,aes(x=logk,y=negcluster1resid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ggsave('/Users/kahinim/Desktop/cluster1neg.png')

ddata$poscluster2resid<-poscluster2mask_nologk$residuals+mean(ddata$mask2pos)
ggplot(ddata,aes(x=logk,y=poscluster2resid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ggsave('/Users/kahinim/Desktop/cluster2pos.png')

ddata$negcluster2resid<-negcluster2mask_nologk$residuals+mean(ddata$mask2neg)
ggplot(ddata,aes(x=logk,y=negcluster2resid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1)+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))
ggsave('/Users/kahinim/Desktop/cluster2neg.png')

```
 _Note: please ignore the 'mean_val_csvs', 'lthr.sh', 'uthr.sh', 'mean_val.sh' in `mehtareplicaten293/...` if you see them - these were part of a different approach we did not end up going with._
 
 ### 7. Final visualizations - missing from original pncitc.md
 
 1. Used this code on `cbica` to transform all niftis I would want to use for visualization to MNI space (as done in step 6): 
```
    ...: zstats=['zfdr1','zfdr2']
    ...: viewim=[]
    ...: for i in range(len(zstats)):
    ...:     at.inputs.input_image = flame1dir + zstats[i]+'.nii.gz'
    ...:     at.inputs.output_image = flame1dir + zstats[i]+'MNI.nii.gz'
    ...:     at.run()
    ...: 
    ...: flame2dir='/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/
    ...: mask2/logk/zfdrmask2/'
    ...: zstats=['zfdr1','zfdr2']
    ...: viewim=[]
    ...: for i in range(len(zstats)):
    ...:     at.inputs.input_image = flame2dir + zstats[i]+'.nii.gz'
    ...:     at.inputs.output_image = flame2dir + zstats[i]+'MNI.nii.gz'
    ...:     at.run()
    ...: 
    ...: clusterdirectory = '/cbica/projects/pncitc/mehtareplicaten293/cluster_
    ...: output/cluster_Z3.09'
    ...: zstats=['/mask1/mask1_2mm','/mask2/mask2_2mm']
    ...: for i in range(len(zstats)):
    ...:     at.inputs.input_image = clusterdirectory + zstats[i]+'.nii.gz'
    ...:     at.inputs.output_image = clusterdirectory + zstats[i]+'MNI.nii.gz'
    ...:     at.run()
```
2. Downloaded pysurfer locally (this was extremely painful, but apparently setup differs on each device. I still wasn't able to view images via the renderer, and had to save them out instead)
3. Downloaded all the MNI space niftis converted in step 1... 
4. Used the code below for visualization (changed minimum and maximum, but you can see the values I assigned at the bottom of each image. These values are based off what Azeez had in his manuscript)

```
import os
from surfer import Brain, project_volume_data #pysurfer imports as surfer

print(__doc__)

"""
Bring up the visualization window.
"""
brain = Brain("fsaverage", "lh", "inflated", subjects_dir='/Users/kahinim/Desktop/freesurfer/subjects')

"""
Get a path to the volume file.
"""
volume_file = "/Users/kahinim/Desktop/zfdr2.nii.gz"

"""
There are two options for specifying the registration between the volume and
the surface you want to plot on. The first is to give a path to a
Freesurfer-style linear transformation matrix that will align the statistical
volume with the Freesurfer anatomy.

Most of the time you will be plotting data that are in MNI152 space on the
fsaverage brain. For this case, Freesurfer actually ships a registration matrix
file to align your data with the surface.
"""
reg_file = "/Users/kahinim/Desktop/freesurfer/average/mni152.register.dat"
zstat = project_volume_data(volume_file, "lh", reg_file)

"""
Note that the contours of the fsaverage surface don't perfectly match the
MNI brain, so this will only approximate the location of your activation
(although it generally does a pretty good job). A more accurate way to
visualize data would be to run the MNI152 brain through the recon-all pipeline.

Alternatively, if your data are already in register with the Freesurfer
anatomy, you can provide project_volume_data with the subject ID, avoiding the
need to specify a registration file.

By default, 3mm of smoothing is applied on the surface to clean up the overlay
a bit, although the extent of smoothing can be controlled.
"""
#zstat = project_volume_data(volume_file, "lh",
                            #subject_id="fsaverage", subjects_dir = "/Users/kahinim/Desktop/freesurfer",smooth_fwhm=0.5)

"""
Once you have the statistical data loaded into Python, you can simply pass it
to the `add_overlay` method of the Brain object.
"""

brain.add_overlay(zstat, min=4, max=9)
brain.show_view('med')
brain.save_image('/Users/kahinim/Desktop/test1.png')
brain.show_view('lat')
brain.save_image('/Users/kahinim/Desktop/test2.png')
brain.show_view('ros')
brain.save_image('/Users/kahinim/Desktop/test3.png')
brain.show_view('caud')
brain.save_image('/Users/kahinim/Desktop/test4.png')
#brain.save_imageset(brain, ['med', 'lat', 'ros', 'caud'], '.png')
"""
It can also be a good idea to plot the inverse of the mask that was used in the
analysis, so you can be clear about areas that were not included.

It's good to change some parameters of the sampling to account for the fact
that you are projecting binary (0, 1) data.
"""


```

5. For the clusters, I used: 

```
import os
from surfer import Brain, project_volume_data

print(__doc__)

"""
Bring up the visualization window.
"""
brain = Brain("fsaverage", "rh", "inflated", subjects_dir='/Users/kahinim/Desktop/freesurfer/subjects')

"""
Get a path to the volume file.
"""
volume_file = "/Users/kahinim/Desktop/mask1/mask1MNI.nii.gz"

"""
There are two options for specifying the registration between the volume and
the surface you want to plot on. The first is to give a path to a
Freesurfer-style linear transformation matrix that will align the statistical
volume with the Freesurfer anatomy.

Most of the time you will be plotting data that are in MNI152 space on the
fsaverage brain. For this case, Freesurfer actually ships a registration matrix
file to align your data with the surface.
"""
reg_file = "/Users/kahinim/Desktop/freesurfer/average/mni152.register.dat"
zstat = project_volume_data(volume_file, "rh", reg_file)

"""
Note that the contours of the fsaverage surface don't perfectly match the
MNI brain, so this will only approximate the location of your activation
(although it generally does a pretty good job). A more accurate way to
visualize data would be to run the MNI152 brain through the recon-all pipeline.

Alternatively, if your data are already in register with the Freesurfer
anatomy, you can provide project_volume_data with the subject ID, avoiding the
need to specify a registration file.

By default, 3mm of smoothing is applied on the surface to clean up the overlay
a bit, although the extent of smoothing can be controlled.
"""
"""
Once you have the statistical data loaded into Python, you can simply pass it
to the `add_overlay` method of the Brain object.
"""
mask_file = volume_file
mask = project_volume_data(mask_file, "rh", subject_id="fsaverage",
                           smooth_fwhm=1, projsum="max").astype(bool)
mask = ~mask
#mask = ~(mask)
stat = mask 

brain.add_overlay(stat, min = 0, max=1)
brain.show_view('med')
brain.save_image('/Users/kahinim/Desktop/test1.png')
brain.show_view('ros')
brain.save_image('/Users/kahinim/Desktop/test2.png')
brain.show_view('caud')
brain.save_image('/Users/kahinim/Desktop/test3.png')
brain.show_view('lat')
brain.save_image('/Users/kahinim/Desktop/test4.png')
#brain.save_imageset(brain, ['med', 'lat', 'ros', 'caud'], '.png')
"""
It can also be a good idea to plot the inverse of the mask that was used in the
analysis, so you can be clear about areas that were not included.

It's good to change some parameters of the sampling to account for the fact
that you are projecting binary (0, 1) data.
"""
```
where the clusters were the masks I'd found in `cluster_Z309`.

I uploaded all these images to the `dropbox` under `manuscript_figures`
