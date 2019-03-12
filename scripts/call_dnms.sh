 bcftools view -Ou -c 1 -C 1 -v snps -I -f PASS /net/topmed8/working/call_sets/out.2017_05/release/minDP10/freeze65k.chr22.pass_and_fail.anno.gtonly.minDP10.vcf.gz | bcftools view -c 1 -S ~/topmed_freeze5_trio_children_nr.txt > ~/topmed_freeze5_trio_children_dnms.vcf
 
 # bcftools view -Ou -S ~/topmed_freeze5_trio_children_nr.txt -v snps -I -f PASS /net/topmed8/working/call_sets/out.2017_05/release/minDP10/freeze65k.chr22.pass_and_fail.anno.gtonly.minDP10.vcf.gz | bcftools view -I -c 1 -C 1  > ~/topmed_freeze5_trio_children_dnms_v2.vcf
 
 # remove duplicates from ped file
 sed -E "s/(,[A-Z0-9,]+*?)//" /net/topmed8/working/call_sets/out.2017_05/index/freeze5.64960.4606.ped > ~/freeze5.nc.ped
 
 
 # get global ancestry of children in pedigrees
 grep -Fwf ~/topmed_freeze5_trio_children_nr.txt /net/topmed2/working/klemmerr/global_list_b38all.txt > ~/topmed_freeze5_trio_children_nr_ancestry.txt