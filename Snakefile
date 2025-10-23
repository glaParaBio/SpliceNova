import pandas

STRAND = ['forward', 'reverse']

def get_fastq(wc):
    fq = [ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r1.iloc[0]]
    if not ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r2.isna().values.any():
        fq2 = ssfq[ssfq.fastq_base == wc.fastq_base].fastq_r2.iloc[0]
        fq.append(fq2)
    return fq

species = pandas.read_csv(config['species'], sep= '\t', comment= '#')

ssfq = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#')
ss = ssfq.drop(['fastq_r1', 'fastq_r2'], axis= 1).drop_duplicates()
assert len(ss.sample_id) == len(set(ss.sample_id))

for fq in ['fastq_r1', 'fastq_r2']:
    ssfq[fq] = [os.path.join(config['fastqdir'], x) for x in ssfq[fq]]

ssfq['fastq_base'] = [re.sub('\.fastq\.gz$|\.fq.\.gz', '', os.path.basename(x)) for x in ssfq.fastq_r1]
assert len(ssfq.fastq_base) == len(set(ssfq.fastq_base))

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in ss.sample_id]),
    fastq_base= '|'.join([re.escape(x) for x in ssfq.fastq_base]),
    species=  '|'.join([re.escape(x) for x in ss.species.unique()]),
    STRAND=  '|'.join([re.escape(x) for x in STRAND]),

rule all:
    input:
        expand('{species}/bigwig/{sample}.forward.bw', zip, species= ss.species, sample= ss.sample_id),
        expand('{species}/bigwig/{sample}.reverse.bw', zip, species= ss.species, sample= ss.sample_id),
        expand('{species}/edger/junction_counts.tsv.gz', species= ss.species),
        expand('{species}/edger/differential_junctions.tsv.gz', species= ss.species),
        # 'multiqc/fastqc_report.html',


include: 'workflows/prep_ref.smk'
include: 'workflows/quant_splice_junctions.smk'


rule sra_download:
    output:
         temp('sra/{sra_id}_1.fastq.gz'),
         temp(touch('sra/{sra_id}_2.fastq.gz')),
    shell:
        r"""
        parallel-fastq-dump --gzip --split-files --outdir `dirname {output[0]}` \
            --threads 8 --tmpdir . --sra-id {wildcards.sra_id}
        """


rule hisat2:
    input:
        fq= get_fastq,
        idx= lambda wc: 'ref/{species}.8.ht2',
    output:
        bam= temp('{species}/hisat2/{fastq_base}.bam'),
        hlog= '{species}/hisat2/{fastq_base}.log',
        cutlog= '{species}/cutadapt/{fastq_base}.log',
    params:
        max_intron_len= config['max_intron_len'],
    run:
        if len(input.fq) == 1:
            cutadapt = f'cutadapt --quality-cutoff 15 --minimum-length 10 --cores 4 -a AGATCGGAAGAGC -o {output.bam}.R1.fq {input.fq[0]} > {output.cutlog}'
            args= f'-U {output.bam}.R1.fq'
        elif len(input.fq) == 2:
            cutadapt = f'cutadapt --quality-cutoff 15 --minimum-length 10 --cores 4 -a AGATCGGAAGAGC -o {output.bam}.R1.fq -A AGATCGGAAGAGC -p {output.bam}.R2.fq {input.fq[0]} {input.fq[1]} > {output.cutlog}'
            args= f'-1 {output.bam}.R1.fq -2 {output.bam}.R2.fq'
        else:
            raise Exception('Unexpected number of fastq files')
        
        idx= re.sub('\.8\.ht2$', '', input.idx)

        shell(f"""
        {cutadapt}

        hisat2 --summary-file {output.hlog} --new-summary --fr --rna-strandness RF \
           --max-intronlen {params.max_intron_len} --threads 4 -x {idx} {args} \
        | samtools view -u -@ 4 \
        | samtools sort -@ 8 > {output.bam}

        rm {output.bam}.R[12].fq
        """)


rule merge_hisat2:
    input:
        bam= lambda wc: expand('{{species}}/hisat2/{fastq_base}.bam', fastq_base= ssfq[ssfq.sample_id == wc.sample].fastq_base),
    output:
        bam= '{species}/hisat2/{sample}.bam',       
    run:
        if len(input.bam) == 1:
            shell("mv {input.bam} {output.bam}")
        else:
            shell("samtools merge {output.bam} {input.bam}")
        shell("samtools index -@ 4 {output.bam}")


rule bigwig:
    input:
        bam= '{species}/hisat2/{sample}.bam',       
    output:
        bw= '{species}/bigwig/{sample}.{strand}.bw',
    shell:
        r"""
        bamCoverage --filterRNAstrand {wildcards.strand} -b {input.bam} -o {output} \
            --binSize 25 \
            --minMappingQuality 5 \
            --normalizeUsing BPM \
            --numberOfProcessors 4
        """


base_fq= {}
for i, row in ssfq.iterrows():
    key = os.path.basename(re.sub('\.gz$', '', row['fastq_r1']))
    key = re.sub('\.fastq$|\.gz$', '', key)
    base_fq[key] = row['fastq_r1']
    if not pandas.isna(row['fastq_r2']):
        key = os.path.basename(re.sub('\.gz$', '', row['fastq_r2']))
        key = re.sub('\.fastq$|\.gz$', '', key)
        base_fq[key] = row['fastq_r2']

rule fastqc:
    priority: -10
    input:
        lambda wc: base_fq[wc.base_fq],
    output:
        'fastqc/{base_fq}_fastqc.zip',
    shell:
        r"""
        fastqc -o fastqc {input}
        """

rule multiqc_fastqc:
    priority: -10
    input:
        expand('fastqc/{base_fq}_fastqc.zip', base_fq= base_fq.keys()),
    output:
        'multiqc/fastqc_report.html',
    shell:
        r"""
        multiqc --force --outdir `dirname {output}` \
            --filename `basename {output}` {input}
        """
