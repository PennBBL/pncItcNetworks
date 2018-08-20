# multivariate analysis with flameo12 

# design  and contrast file 
projectdir=/data/jux/BBL/projects/pncItcNetworks
design=$projectdir/scripts/dti
maskfile=/data/joy/BBL/studies/pnc/template/pnc_template_brain_2mm_MASK.nii.gz

mkdir $projectdir/output/dti


outputdir1=$projectdir/output/dti

rm -rf $outputdir1/stats
rm -rf $outputdir1/stats_dememean

#merge the files to 4D

fslmerge -t $outputdir1/dti4D.nii.gz  $projectdir/output/dtidensity/*tractdensity.nii.gz

flameo --copefile=$outputdir1/dti4D.nii.gz  --mask=$maskfile  --dm=$design/design.mat --tc=$design/contrast.con --cs=$design/grp.grp --runmode=flame1 --ld=$outputdir1/stats 



flameo --copefile=$outputdir1/dti4D.nii.gz  --mask=$maskfile  --dm=$design/design_demean.mat --tc=$design/contrast.con --cs=$design/grp.grp --runmode=flame1 --ld=$outputdir1/stats_demean



