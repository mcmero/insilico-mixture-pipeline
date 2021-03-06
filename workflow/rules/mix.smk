rule preprocess:
    input:
        'input/{sample}.bam'
    output:
        'results/preproc/{sample}_{norm_seed}.{norm}.bam'
    log:
        'logs/preprocess_{sample}_{norm_seed}.{norm}.log'
    threads:
        8
    resources:
        mem_mb=65536,
        runtime='1-0:0:0'
    envmodules:
        'samtools/1.9'
    shell:
        '''
        if [ {wildcards.norm} -eq 0 ]; then
            ln -s ../../{input} {output} ;
        else
            samtools view -@ {threads} -s {wildcards.norm_seed}.{wildcards.norm} -b {input} > {output} ;
        fi
        '''

rule subsample:
    input:
        'results/preproc/{sample}_{norm_seed}.{norm}.bam'
    output:
        'results/subsampled/{norm_seed}.{norm}_{seed}/{sample}_0.{prop}.bam',
    log:
        'logs/subsample_{norm_seed}.{norm}_{seed}_{sample}_0.{prop}.log'
    threads:
        8
    resources:
        mem_mb=65536,
        runtime='1-0:0:0'
    envmodules:
        'samtools/1.9'
    shell:
        'samtools view -@ {threads} -s {wildcards.seed}.{wildcards.prop} -b {input} > {output}'

rule merge:
    input:
        bams=[
            'results/subsampled/{norm_seed}.{norm}_{seed}/{sample}_0.{prop}.bam'.format(
                sample=sample, norm_seed=norm_seed, norm=norm, seed=seed, prop=prop
             ) for sample, norm_seed, norm, seed, prop in zip(samples, norm_seeds, norms, seeds, props)
        ]
    output:
        'results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam'
    log:
        'logs/merge_{norm_seed1}.{norm1}_{norm_seed2}.{norm2}_{sample1}_0.{prop1}_{sample2}_0.{prop2}.log'
    threads:
        8
    resources:
        mem_mb=65536,
        runtime='2-0:0:0'
    envmodules:
        'samtools/1.9'
    params:
        mem_per_thread='7G'
    shell:
        '''
        samtools cat \
            results/subsampled/{wildcards.norm_seed1}.{wildcards.norm1}_{wildcards.norm_seed1}/{wildcards.sample1}_0.{wildcards.prop1}.bam \
            results/subsampled/{wildcards.norm_seed2}.{wildcards.norm2}_{wildcards.norm_seed2}/{wildcards.sample2}_0.{wildcards.prop2}.bam \
            | samtools sort -m {params.mem_per_thread} -@ {threads} - -o {output}
        '''

rule index_merged:
    input:
        'results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam'
    output:
        'results/merged/{norm_seed1}.{norm1}_{norm_seed2}.{norm2}/{sample1}_0.{prop1}_{sample2}_0.{prop2}.bam.bai'
    log:
       'logs/index_merged_{norm_seed1}.{norm1}_{norm_seed2}.{norm2}_{sample1}_0.{prop1}_{sample2}_0.{prop2}.log'
    threads:
        8
    resources:
        mem_mb=65536,
        runtime='1-0:0:0'
    envmodules:
        'samtools/1.9'
    shell:
        'samtools index -@ {threads} {input}'
