---
title: "PNC-ITC Replication"
authors: Kahini
output: html_document
---
#### Copied and modified from pncitc.md


 
# 1. CWASMDMR

The computation of  CWASMDMR was  done with  `cwasmdr` singularity image (`/cbica/projects/pncitc/cwasmdmr.simg`). We used packages from the connectir project at [https://github.com/czarrar/connectir](https://github.com/czarrar/connectir)

Distance matrix was first computed with the following script: 

```
#!/bin/bash
#$ -l h_vmem=320G #QSUB, can take some hours
#$ -l tmpfree=200G

singimage=/cbica/projects/pncitc/cwasmdmr.simg 
scriptdir=/usr/local/bin


mdmrouput=/cbica/projects/pncitc/mehtareplicate/cwas307 #output directory
brainmask=/cbica/projects/pncitc/subjectData/PNCgrey.nii.gz # greymatter mask from pnc   
bgim=/cbica/projects/pncitc/subjectData/PNCbrain.nii.gz # pnc template from pnc

imagelist=/cbica/projects/pncitc/subjectData/imageinput_rest2.csv # list of image in nifti # HAD TO RE-GENERATE THIS LIST AS THE FILE PATHS HAD CHANGED


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
#!/bin/bash #QSUB, can take some hours
#$ -l h_vmem=300G
#$ -l tmpfree=300G
singularity exec -e -B /cbica/projects/pncitc  \
/cbica/projects/pncitc/cwasmdmr.simg \
/usr/local/bin/Rscript /usr/local/bin/connectir_mdmr.R -i /cbica/projects/pncitc/mehtareplicate/cwas307 -f 'logk+relMeanRMSmotion+sex+age' -m /cbica/projects/pncitc/demographics/n307_demographics.csv --factors2perm='logk' --save-perms -c 5 -t 5  --ignoreprocerror --memlimit=300 logk_motion_sex_age
```

Memory and formatting were the main problems with these scripts not running well - the numbers/format left in were what worked for me. The output is at: `/cbica/projects/pncitc/mehtareplicate/cwas307/logk_motion_sex_age`

### 2. Significant clusters from mdmr
The cluster analysis was computed  with the script `scripts/grf_fslcluster.sh`, written based on  [FSL cluster analysis](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Cluster) with  Gaussian Random Field (GRF) theory. 

The script `scripts/cluster.sh` called grf_fslcluster.sh, as listed in this repo under scripts, with `z=3.09`
Two clusters were obtained: one at the at the frontal region the order at the TPJ. 

cluster.sh reads as below:
```
 #!/bin/bash #NO NEED TO QSUB THIS
dir=/cbica/projects/pncitc
bash grf_fslcluster.sh -i ${dir}/mehtareplicate/test.nii  -m ${dir}/mehtareplicate/cwas307/mask.nii.gz -t 3.09 -o ${dir}/mehtareplicate/cluster_output
```


The output of cluster masks is at: `/cbica/projects/pncitc/mehtareplicate/cluster_output/cluster_Z3.09`. 

Numbers obtained from the CSV slightly different than before, but ultimately a close replication (likely due to changes in software version): 
| Cluster Index | Voxels | P | -log10(P) | MAX | MAX X (vox) | MAX Y (vox) | MAX Z (vox) | COG X (vox) | COG Y (vox) | COG Z (vox) |
| ---- | ---- | ---- | ------ | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| 2 | 11 | 9.27E-05 | 4.03 | 3.54 | 30 | 44 | 30 | 30.9 | 43.9 | 30.4 |
| 1 | 5 | 0.0228 | 1.64	 | 3.54 | 13 | 28 | 24 | 13.4 | 28.6 | 24.2 |

However, the output image of clusters produced was an .img and .hdr which needed to be turned into a nifti format. Used matlab with [this](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) toolbox and these commands on the output images (the .img and .hdr files had to be in the same directory for this to work): 

```
nii=load_nii('cluster_output/cluster_Z3.09/cluster_Z3.09.hdr')
save_nii(nii, 'cluster_Z3.09.nii');
```
*NOTE: ```fslchfiletype NIFTI_GZ cluster_Z3.09.img cluster_Z3.09.nii``` achieves the same thing and is faster*


Notably, the images were very sparse and hard to see - if viewed in mrpeek, appeared to be empty images with "NaN" scales - in Mango, only one cluster was visible. Ultimately downloaded MRICRON to view images locally - clusters appeared in basically the same location they did in previous replications. 



### 3. Seed-based correlation 
Two different masks were generated from the cluster_z3_09.nii, using fslmath at/cbica/projects/pncitc/mehtareplicate/cluster_output/cluster_Z3.09/masks:

(See [https://mandymejia.com/fsl-maths-commands/](https://mandymejia.com/fsl-maths-commands/))

```
fslmaths cluster_Z3.09.nii.gz -thr 1 -uthr 2 mask1.nii.gz #CLUSTER1
fslmaths cluster_Z3.09.nii.gz -thr 2 -uthr 2 mask2.nii.gz #CLUSTER2
```

Masks generated were again in .hdr and .img format, so I used fslchfiletype to turn them into niftis.

The two  masks were upsampled from 4mm to 2mm and were used as seeds for seed-based correlation.

This upsampling was done using the pnc_2mm template, found at `/cbica/projects/pncitc/mehtareplicate/cluster_output/cluster_Z3.09/masks/pnc_template_brain_2mm.nii.gz`

The code used for this was: 
```
3dresample -master pnc_template_brain_2mm.nii.gz -input mask1.nii.gz -prefix mask1_2mm.nii.gz
3dresample -master pnc_template_brain_2mm.nii.gz -input mask2.nii.gz -prefix mask2_2mm.nii.gz
```

The path to the two seeds is: 
`/cbica/projects/GURLAB/projects/pncitc/output/cluster_Z3.09/mask1_2mm.nii.gz `
`/cbica/projects/GURLAB/projects/pncitc/output/cluster_Z3.09/mask2_2mm.nii.gz` 

The seed-based correlation was computed with the following script, and the `xcpengine.simg` file under `/cbica/projects/pncitc/mehtareplicate`:

```
#!/bin/bash # DEFINITELY QSUB THIS ONE
cd /cbica/projects/pncitc/mehtareplicate
XCPEDIR=xcpEngine
seedpoint1=/cbica/projects/pncitc/mehtareplicate/cluster_output/cluster_Z3.09/masks/mask1_2mm.nii.gz
seedpoint2=/cbica/projects/pncitc/mehtareplicate/cluster_output/cluster_Z3.09/masks/mask2_2mm.nii.gz

bblid=/cbica/projects/pncitc/demographics/n307_bblid_scandid.csv
image=/cbica/projects/pncitc/subjectData/rest/
outputdir=/cbica/projects/pncitc/mehtareplicate/seedcorrmaps

mkdir -p ${outputdir}

cat $bblid | while IFS="," read -r a b ; 

do
     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/mehtareplicate /cbica/projects/pncitc/mehtareplicate/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint1 -o $outputdir -p ${a},${b} -k 6 -n mask1\
     #rm -rf $outputdir/seed/mask1/${a}_${b}_connectivity_mask1_seed.nii.gz


     img=$(ls -f $image/${a}_${b}_rest.nii.gz)
     singularity exec --cleanenv -B /cbica/projects/pncitc/mehtareplicate /cbica/projects/pncitc/mehtareplicate/xcpengine.simg /xcpEngine/utils/seedconnectivity -i $img -s $seedpoint2 -o $outputdir -p ${a},${b} -k 6 -n mask2\
     #rm -rf $outputdir/seed/mask2/${a}_${b}_connectivity_mask2_seed.nii.gz
done
```
*NOTE: For this script, the output displays in the error log for some reason*
*Also note, the demographic csvs say "scandid", not "scanid"

### 4. Linear regression with FSL `flameo` 

Flameo regression computation requires `design`,`contrast` and `group` text files. The script `scripts/makeflameodesig.R` was used to generate these, but the package was out of date in R on CBICA, so I ran it locally and copied the needed files over into `cbica/projects/pncitic/mehtareplicate/regression`

_Note: You can install the required packages in a required PERSONAL library and thus run the script on CUBIC_

I also converted the designlogkmat from .txt to .mat and did the same for grp.txt, using: 

```
Text2Vest desigmatlogkonly.txt desigmatlogkonly.mat
Text2Vest grp.txt grp.grp
Text2Vest contrast4.txt contrast4.con
```

The flameo linear regression was computed with this script:

```
#!/bin/bash
#$ -l h_vmem=200G
#$ -l tmpfree=200G # this will throw an error about ": $PATH does not agree with $PATH_modshare counter.", but it shouldn't stop the process. 
bblid=/cbica/projects/pncitc/demographics/n307_bblid_scandid.csv #FIRST HALF RUNS QUICKLY, WOULD QSUB BECAUSE THE SECOND PART TAKES A BIT LONGER
imagedir=/cbica/projects/pncitc/mehtareplicate/seedcorrmaps/seed/
scriptdir=/cbica/projects/pncitc/mehtareplicate/regression
outputdir=/cbica/projects/pncitc/mehtareplicate/regression
demogdir=/cbica/projects/pncitc/mehtareplicate/regression


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

I also created the flameo csvs pointing to the mask1Z_sm6.nii.gz niftis in the same directory, naming the csvs `mask1.csv` and `mask2.csv` respectively. 

The outputs of flameo regression: 

`/cbica/projects/pncitc/mehtareplicate/regression/mask1`
`/cbica/projects/pncitc/mehtareplicate/regression/mask2`

In each directory, there are zvalues: 

  `zstat1 : average `

  `zstat2 : logk `

  `zstat3 : sex`

  `zstat4 : motion`

  `zstat5 : age`

# to be updated

The zstats  were FDR corrected with this script `scripts/flameoutputfdrcorrection`. The outputs are located here:

`/cbica/projects/pncitc/mehtareplicate/seedcorrmaps/mask1/logk`
`/cbica/projects/pncitc/mehtareplicate/seedcorrmaps/mask2/logk`

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





 
