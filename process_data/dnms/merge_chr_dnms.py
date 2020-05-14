import os
import glob

sourcedir = "/net/snowwhite/home/jedidiah/"
vcfdir = sourcedir + "topmed_freeze5_dnms/"
outdir = vcfdir + "merged/"
pythonpath = sourcedir + "anaconda3/envs/anno/bin/python "
bcftoolspath = sourcedir + "anaconda3/envs/anno/bin/bcftools "


for i in range(1,23):
    chrom = "chr" + str(i)
    # files = " <(ls " + vcfdir + chrom + ".*) > "
    files = " ".join(glob.glob(vcfdir + chrom + '.*'))
    # if i in [1,16]:
    cmd = bcftoolspath + "concat " + files + " > " + outdir + chrom + ".freeze5.dnms.vcf"
    os.system(cmd)
    