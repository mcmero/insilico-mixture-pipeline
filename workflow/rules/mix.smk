rule preprocess:
    input:
        'input/{sample}.bam'
    output:
        'results/preproc/{sample}_{pseed}.{norm}.bam'
    threads:
        8
    envmodules:
        'gcc/8.3.0',
        'samtools/1.9'
    shell:
        '''
        if [ {wildcards.norm} -eq 0 ]; then
            ln -s ../../{input} {output} ;
        else
            samtools view -@ {threads} -s {wildcards.pseed}.{wildcards.norm} -b {input} > {output} ;
        fi
        '''

rule subsample:
    input:
        'results/preproc/{sample}_{pseed}.{norm}.bam'
    output:
        'results/subsampled/{sample}_{pseed}.{norm}_{seed}.{prop}.bam',
    threads:
        8
    envmodules:
        'gcc/8.3.0',
        'samtools/1.9'
    shell:
        'samtools view -@ {threads} -s {wildcards.seed}.{wildcards.prop} -b {input} > {output}'

rule sort:
    input:
        'results/subsampled/{sample}_{pseed}.{norm}_{seed}.{prop}.bam',
    output:
        'results/subsampled-sorted/{sample}_{pseed}.{norm}_{seed}.{prop}.bam'
    threads:
        8
    envmodules:
        'gcc/8.3.0',
        'samtools/1.9'
    shell:
        'samtools sort -@ {threads} -m {config[mem]} {input} -o {output}'

rule index:
    input:
        'results/subsampled-sorted/{sample}_{pseed}.{norm}_{seed}.{prop}.bam'
    output:
        'results/subsampled-sorted/{sample}_{pseed}.{norm}_{seed}.{prop}.bam.bai'
    threads:
        8
    envmodules:
        'gcc/8.3.0',
        'samtools/1.9'
    shell:
        'samtools index -@ {threads} {input}'

rule merge:
    input:
        bams=[
            'results/subsampled-sorted/{sample}_{pseed}.{norm}_{seed}.{prop}.bam'.format(
                sample=sample, pseed=pseed, norm=norm, seed=seed, prop=prop
             ) for sample, pseed, norm, seed, prop in zip(samples, pseeds, norms, seeds, props)
        ],
        indexes=[
            'results/subsampled-sorted/{sample}_{pseed}.{norm}_{seed}.{prop}.bam.bai'.format(
                sample=sample, pseed=pseed, norm=norm, seed=seed, prop=prop
             ) for sample, pseed, norm, seed, prop in zip(samples, pseeds, norms, seeds, props)
        ],
    output:
        'results/merged/{sample1}_{seed1}.{prop1}_{sample2}_{seed2}.{prop2}.bam'
    threads:
        8
    envmodules:
        'gcc/8.3.0',
        'samtools/1.9'
    shell:
        'samtools merge -n -@ {threads} {output} {input.bams}'
