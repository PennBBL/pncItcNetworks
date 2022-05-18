---
title: "PNC-ITC Replication"
authors: Kahini
output: html_document
---
#### Copied and modified from pncitc.md

Sample replication: n307 did not exclude those with the `health_exclude` criteria. Analyses were re-run on n293, using Pehlivanova et al's n427 sample and then running restQA exclusions on them (all this information was from Pehlivanova .csvs). The code and .csvs for this are available at `cbica/projects/pncitc/samplerecreation`. All subsequent analyses were run in `cbica/projects/pncitc/mehtareplicaten293`, using the same steps as below in folder `mehtareplicaten293` - additionally, any new .csvs should be pointed to in the scripts in that directory. Any departure from the protocol for n307 will be noted here.  I also moved the bblid_scanid .csv to demographics, and created a folder within subjectData called rest293 for the n = 293 replication.
**Results were similar in N293 for the second cluster, things changed for the first cluster. Visualizations are available in the .html format within mehtareplicaten293**

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
bash grf_fslcluster.sh -i ${dir}/mehtareplicate/test.nii #this was the "zstats_logk.nii.gz" within the log k file, just moved out and renamed
-m ${dir}/mehtareplicate/cwas307/mask.nii.gz -t 3.09 -o ${dir}/mehtareplicate/cluster_output
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
`~/cluster_Z3.09/masks/mask1_2mm.nii.gz `
`~/cluster_Z3.09/masks/mask2_2mm.nii.gz` 

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


The zstats  were FDR corrected with this script `scripts/flameoutputfdrcorrection`. For this, 
1. I opened 'R' in CBICA
2. Did:


```
setwd('/cbica/projects/pncitc/mehtareplicate/regression/')
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

`/cbica/projects/pncitc/mehtareplicate/regression/mask1/logk`
`/cbica/projects/pncitc/mehtareplicate/regression/mask2/logk`

FDR corrected z-values. 

  `zfdr1 : average `

  `zfdr2 : logk `

  `zfdr3 : sex`

  `zfdr4 : motion`

  `zfdr5 : age`
  

**Changes: ran locally for n293 due to library difficulties on CBICA**

### 5. Vizualisation of Results - results pretty similar to first iteration, based off images in the notebooks in this repo
All computations were done in PNC template. For vizualisation, all the nifti files  were tranformed to MNI before. 

  a. for clusters and mean of seed-based correlation (I ran this one last, as I had trouble locating the script (it's now `cbica/projects/pncitc/dropbox/seed-basedcorrelation.ipynb`), hence the high value for ii) : 
At the end of the script, I added these lines to write out each output to .html: 
  ```
  ii = 4
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
  ```
  _Note: had to use flchfiletype on the copeseed images, mask1and2_2mm.nii, mean copeseed images, and move the .hdr and .img files out/ remove them altogether - having the nifti and img in the same directory can cause an error._
  
  b. for mask1 (ran in iPython on CBICA): `notebook/flameomask1.ipynb`
  
  At the end of the script, I added these lines to write out each output to .html: 
  ```
  ii = 0
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
  ```

  c. for mask2 (iPython) : `notebook/flameomask2.ipynb`

 At the end of the script, I added these lines to write out each output to .html: 
  ```
  ii = 2
for x in viewim:
    ii+=1
    x.save_as_html("visualization"+str(ii)+".html")
  ```
_Note: The test_vis folder uses the same code but with threshold 0 as those images are more similar to the ones in the manuscript_

Images I generated were saved in `cbica/projects/pncitc/mehtareplicate/KMVis/full_cortical_vis` as AA(produced by Azees)or KM(me). They are numbered to correspond.

### 6. Regional plot of significant regions of logk - results pretty similar to first iteration
 
 I did this step locally as the R libraries on CBICA were very out of date. 
 
The positive and negative zvalues of seed-based correlation regression with logk was extracted with the script `scripts/extractsignificantcluster.R` for both seed masks. For some reason, changing the loop to go from 2 to 308 (rather than 1 to 307) fixed errors I had, but the last subject kept getting left out and the first row was always populated with zeros, so I simply ran the last participant by itself and copy-pasted their results into the first row:

Other than changing the paths in general, I made this change to write out the one subject: 
Original script: 
```
for (i in 1:307 ) {
     img1=readNifti(paste0('/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/seed/mask1/',b[i,1],'_',b[i,2],'_connectivity_mask1Z_sm6.nii.gz'))
     img2=readNifti(paste0('/cbica/projects/GURLAB/projects/pncitc/output/seedcorrmaps/seed/mask2/',b[i,1],'_',b[i,2],'_connectivity_mask2Z_sm6.nii.gz'))
     datap1=img1[p_m1==1]; datap2=img2[p_m2==1]
     datam1=img1[n_m1==1]; datam2=img2[n_m2==1]
     corrdata[i,]=c(b[i,1],b[i,2],mean(datap1),mean(datam1),mean(datap2),mean(datam2))
}
```

Changes: 

```
img1=readNifti(#path to their nifti)
img2=readNifti(#path to their nifti)
datap1=img1[p_m1==1]; datap2=img2[p_m2==1]
datam1=img1[n_m1==1]; datam2=img2[n_m2==1]
corrdata[1,]=c(#BBLID,#SCANID],mean(datap1),mean(datam1),mean(datap2),mean(datam2)
```


The results were vizualised with `notebook/meanseedcorrelationplot.Rmd`. I also ran this locally, and compared my results (labelled KM) with Azeez's (labelled AA). These can be found at the path `cbica/projects/pncitc/mehtareplicate/KMVis`. They are numbered to correspond. 
 
 **Changes- For N293, I realized the final manuscript used different code, so did not reproduce this step.**

### 5. Vizualisation of Results 

1. I had to use the `zfdr2.nii.gz` as mask on all seed connectivity maps, eg: `121085_7602_connectivity_mask1Z_sm6.nii.gz` for both clusters. Then, I found the average of the z(r) values within those maps and plot them versus the delay discounting log(k) parameter available in hte `demographics.csv`.
2. Azeez did positive and negative values separately, so first I thresholded all negative values/positive values into separate connectivity maps for both clusters. I navigated into `/cbica/projects/pncitc/mehtareplicaten293/seedcorrmaps/seed/mask1` for cluster 1 (or `mask2` for cluster 2). 
3. I ran the following scripts: 
```
#/bin/bash
FILES=$(ls *Z*)
for f in $FILES
do
	fslmaths $f -uthr 0 negthresh_${f}
	trimmed=$(basename negthresh_${f} .nii.gz)
	fslchfiletype NIFTI_GZ ${trimmed}.img ${trimmed}.nii.gz
	rm -rf *img* *hdr*
	echo "Processed negthresh_${f}"  

done
```
and 
```
#/bin/bash
FILES=$(ls *Z*)
for f in $FILES
do
	fslmaths $f -thr 0 posthresh_${f}
	trimmed=$(basename posthresh_${f} .nii.gz)
	fslchfiletype NIFTI_GZ ${trimmed}.img ${trimmed}.nii.gz
	rm -rf *img* *hdr*
	echo "Processed posthresh_${f}"  

done
```
to produce images in the format `posthresh_116812_7087_connectivity_mask1Z_sm6.nii.gz` and `negthresh_112633_7573_connectivity_mask1Z_sm6.nii.gz`. Note that running the second script after the first accidentally generated posthresh_negthresh niftis, which I removed via `rm -rf`. 
3. I then used the script below to pipe mean connectivity values to a text file:
```
general=/data/joy/BBL/tutorials/exampleData/AMICO_NODDI/raw/*/*

for i in $general;do 
	bblIDs=$(echo ${i}|cut -d'/' -f9 |sed s@'/'@' '@g);
	SubDate_and_ID=$(echo ${i}|cut -d'/' -f10|sed s@'/'@' '@g|sed s@'x'@'x'@g)
	filepath=$(echo ${i} | rev | cut -d'/' -f4- | rev )
3dROIstats -mask ~/templates/pnc_wm_prior_bin2mm.nii.gz -1DRformat -nomeanout -nzmean ${filepath}/Processed_Data/${bblIDs}/${SubDate_and_ID}/norm/		*_ODI_Std.nii.gz >>~/ODI_mean_wm.txt
```
This generated the values I needed. 
4. Then, using the `demographics.csv` as well as resting-state QA data from Pehlivanova et al in the `samplerecreation` folder, I was able to regenerate the graphs from the manuscript by adapting the following R code from Adam: 
```
ddata=readRDS('~/Desktop/ITC/my_data.rds'). # replaced with .csv containing necessary info, i.e: age, bblid, scanid, relRMS, sex, age, logK, posmask values for cluster 1, posmask values for cluster 2, negativemask values for cluster 1, negative mask values for cluster2
library(visreg);
postpjmask_nologk=lm(mask1pos~age+sex+relMeanRMSmotion,data=ddata)
postpjmask=lm(mask1pos~logk+age+sex+relMeanRMSmotion,data=ddata)

#summary(postpjmask)
pdf("mask1pos.pdf", width = 8, height = 8)
imageplot<-visreg(postpjmask, "logk", 
                   xlab="Log K", 
                  ylab="Correlation",line=list(col="red",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "red"),
                  fill=list(col=adjustcolor("red", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,ylim=c(-0.5,0.5),xlim=c(-9,-1) )

negtpjmask=lm(mask1neg~age+sex+relMeanRMSmotion+logk,data=ddata)
negtpjmask_nologk=lm(mask1neg~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
pdf("mask1neg.pdf", width = 8, height = 8)
imageplot<-visreg(negtpjmask, "logk", 
                   xlab="Log K", 
                  ylab="Correlation k",line=list(col="blue",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "blue"),
                  fill=list(col=adjustcolor("blue", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )

posfrmask=lm(mask2pos~age+sex+relMeanRMSmotion+logk,data=ddata)
posfrmask_nologk=lm(mask2pos~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
pdf("mask2pos.pdf", width = 8, height = 8)
imageplot<-visreg(posfrmask, "logk", 
                   xlab="Log K", 
                  ylab="Correlation",line=list(col="red",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "red"),fill=list(col=adjustcolor("red", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )

negfrmask=lm(mask2neg~age+sex+relMeanRMSmotion+logk,data=ddata)
negfrmask_nologk=lm(mask2neg~age+sex+relMeanRMSmotion,data=ddata)
#summary(postpjmask)
pdf("mask2neg.pdf", width = 8, height = 8)
imageplot<-visreg(negfrmask, "logk", 
                   xlab="Log K", 
                  ylab="Correlation",line=list(col="blue",lwd=4),overlay=TRUE,rug = FALSE,points.par = list(pch =16, cex = 1, col = "blue"), fill=list(col=adjustcolor("blue", alpha.f = 0.5)), cex.axis=1.5,cex.lab=1.5,xlim=c(-10,0),ylim=c(-.8,.8) )


library(ggplot2)

ylab<-"Correlation (ρ)"

ddata$postpjresid<-postpjmask_nologk$residuals+mean(ddata$mask1pos)
ggplot(ddata,aes(x=logk,y=postpjresid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8) +ylim(c(-0.438,0.3))+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))


ddata$negtpjresid<-negtpjmask_nologk$residuals+mean(ddata$mask1neg)

ggplot(ddata,aes(x=logk,y=negtpjresid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1) +ylim(c(-0.5,0.3))+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ddata$posfrresid<-posfrmask_nologk$residuals+mean(ddata$mask2pos)
ggplot(ddata,aes(x=logk,y=posfrresid)) + geom_smooth(method = 'lm', colour=('#b40101'), fill = "#ef1212",size=2,alpha=.8) +ylim(c(-0.32,0.64))+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))

ddata$negfrresid<-negfrmask_nologk$residuals+mean(ddata$mask2neg)
ggplot(ddata,aes(x=logk,y=negfrresid)) + geom_smooth(method = 'lm', colour=('#0c3e6d'), fill = "#69abde",size=2,alpha=1) +ylim(c(-0.65,0.43))+xlim(c(-8.75,-1))+ geom_point() + xlab("Discount Rate (logK)") +ylab(ylab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'black', size = 2), axis.line.y = element_line(colour = 'black', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(face="bold",size=20), axis.title = element_text(size=26), axis.title.y = element_text(margin = margin(t = 0, r = 27, b = 0, l = 0)))
```

 
