#rule multicov:
#    input:
#        bam='results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam',
#        index='results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam.bai',
#        bed='input/{sample1}_{sample2}_mutations.bed'
#    output:
#        'results/depths/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}_depths.tsv'
#    envmodules:
#        'bedtools/2.30.0'
#    threads:
#        1
#    resources:
#        mem_mb=32768,
#        runtime='1-0:0:0'
#    shell:
#        'bedtools multicov -D -F -bams {input.bam} -bed {input.bed} > {output}'

rule mpileup:
    input:
        bam='results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam',
        index='results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam.bai',
        bed='input/{sample1}_{sample2}_mutations.bed'
    output:
        'results/pileup/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}_pileup.vcf'
    envmodules:
        'gcc/8.3.0',
        'bcftools/1.9'
    threads:
        1
    resources:
        mem_mb=32768,
        runtime='1-0:0:0'
    shell:
        '''
        bcftools mpileup -B -A -d{config[maxdepth]} -f {config[ref]} \
            -q{config[baseQ]} -Q{config[mapQ]} -R {input.bed} -a AD {input.bam} > {output}
        '''
