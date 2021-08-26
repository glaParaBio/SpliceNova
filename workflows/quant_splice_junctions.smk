def featureCountsStrand(wc):
    """Numeric code used by featureCounts to identify strandness of protocol
    for RNAseq
    """
    strand= ss[ss.sample_id == wc.sample].strand
    assert len(strand) == 1
    strand= strand.iloc[0]
    code= None
    if strand == 'FR':
        code= 1
    elif strand == 'RF':
        code= 2
    elif strand == 'unstranded':
        code= 0
    else:
        sys.stderr.write('\n\nInvalid strand identifier found in sample sheet: %s\n' % strand)
        sys.exit(1)
    return code

rule pre_filter_bam:
    input:
        bam= lambda wc: '{species}/hisat2/{sample}.bam',
    output:
        temp('{species}/hisat2/{sample}.pre.bam'),
    shell:
        r"""
        samtools view -@ 2 -u -q 10 -F 2828 {input.bam} > {output}
        """

rule filter_reversestrand:
    input:
        '{species}/hisat2/{sample}.pre.bam',
    output:
        bam= temp('{species}/hisat2/{sample}.reverse.bam'),
        bai= temp('{species}/hisat2/{sample}.reverse.bam.bai'),
    run:
        if ssfq[ssfq.sample_id == wildcards.sample].fastq_r2.isna().any():
            shell(r"""
            samtools view -u -@ 4 -F 16 {input} > {output.bam}
            """)
        else:
            shell(r"""
            samtools view -u -@ 4 -f 64 -F 16 {input} > {output.bam}.1
            samtools view -u -@ 4 -f 144 {input} > {output.bam}.2
            
            samtools merge {output.bam} {output.bam}.1 {output.bam}.2
            rm {output.bam}.1 {output.bam}.2
            """)
        shell("samtools index -@ 4 {output.bam}")


rule filter_forwardstrand:
    input:
        '{species}/hisat2/{sample}.pre.bam',
    output:
        bam= temp('{species}/hisat2/{sample}.forward.bam'),
        bai= temp('{species}/hisat2/{sample}.forward.bam.bai'),
    run:
        if ssfq[ssfq.sample_id == wildcards.sample].fastq_r2.isna().any():
            shell(r"""
            samtools view -u -@ 4 -f 16 {input} > {output.bam}
            """)
        else:
            shell(r"""
            samtools view -u -@ 4 -f 80 {input} > {output.bam}.1
            samtools view -u -@ 4 -f 128 -F 16 {input} > {output.bam}.2
            
            samtools merge {output.bam} {output.bam}.1 {output.bam}.2
            rm {output.bam}.1 {output.bam}.2
            """)
        shell("samtools index -@ 4 {output.bam}")


rule split_strand:
    input:
        gff= 'ref/{species}.gff',
    output:
        plus= temp('ref/{species}.plus.gff'),
        minus= temp('ref/{species}.minus.gff'),
    shell:
        r"""
        awk '$7 == "+"' {input.gff} > {output.plus}
        awk '$7 == "-"' {input.gff} > {output.minus}
        """

jcmd = r"featureCounts -a {input.gff} -o {output.cnt} -g gene_id -s {params.strand} -J -G {input.fa} -p -T 4 {input.bam}"

rule featureCounts_plus:
    input:
        bam= '{species}/hisat2/{sample}.forward.bam',
        gff= 'ref/{species}.plus.gff',
        fa= 'ref/{species}.fasta',
    output:
        cnt= temp('{species}/featureCounts/{sample}.counts.plus'),
        jx= temp('{species}/featureCounts/{sample}.counts.plus.jcounts'),
        sumry= temp('{species}/featureCounts/{sample}.counts.plus.summary'),
    params:
        strand= featureCountsStrand,
    shell:
        jcmd


rule featureCounts_minus:
    input:
        bam= '{species}/hisat2/{sample}.reverse.bam',
        gff= 'ref/{species}.minus.gff',
        fa= 'ref/{species}.fasta',
    output:
        cnt= temp('{species}/featureCounts/{sample}.counts.minus'),
        jx= temp('{species}/featureCounts/{sample}.counts.minus.jcounts'),
        sumry= temp('{species}/featureCounts/{sample}.counts.minus.summary'),
    params:
        strand= featureCountsStrand,
    shell:
        jcmd


rule merge_split_counts:
    input:
        plus= lambda wc: [f'{wc.species}/featureCounts/{x}.counts.plus.jcounts' for x in sorted(ss[ss.species == wc.species].sample_id)],
        minus= lambda wc: [f'{wc.species}/featureCounts/{x}.counts.minus.jcounts' for x in sorted(ss[ss.species == wc.species].sample_id)],
    output:
        jx= '{species}/featureCounts/counts.txt.jcounts',
    run:
        dts = []
        for x in input.plus + input.minus:
            dt = pandas.read_csv(x, sep= '\t')
            fid = dt.columns[len(dt.columns)-1]
            dt.rename(columns= {fid: 'count'}, inplace= True)
            sample_id = re.sub('\.counts\.(plus|minus)\.jcounts$', '', os.path.basename(x))
            dt['sample_id'] = sample_id
            assert sample_id in fid

            if x in input.plus:
                dt = dt[(dt.Site1_strand == '+') & (dt.Site2_strand == '+')]
            elif x in input.minus:
                dt = dt[(dt.Site1_strand == '-') & (dt.Site2_strand == '-')]
            else:
                raise Exception
            
            dts.append(dt)

        sites = pandas.concat(dts)
        sites.to_csv(output.jx, sep= '\t', index= False, na_rep= 'NA')

# At this point you have all the splice junctions quantified according to strand
# We need to count the non-split reads spanning these junctions

rule junctions_to_gff:
    input:
        jx= '{species}/featureCounts/counts.txt.jcounts',
        gff= 'ref/{species}.gff',
        utils= os.path.join(workflow.basedir, 'scripts/utils.R'),
    output:
        gff= '{species}/featureCounts/junctions.gff',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)
source('{input.utils}')

jx <- fread('{input.jx}', select= c('PrimaryGene', 'Site1_chr', 'Site1_strand', 'Site1_location', 'Site2_chr', 'Site2_strand', 'Site2_location'))
jx <- unique(jx)
jx <- jx[!is.na(PrimaryGene) & !is.na(Site1_strand)]
stopifnot(identical(jx$Site1_chr, jx$Site2_chr))
stopifnot(identical(jx$Site1_strand, jx$Site2_strand))

ref <- fread('{input.gff}')
ref[, gene_id := gff_attribute(V9, 'gene_id')]
ref <- unique(ref[gene_id != '', list(gene_id, strand= V7)])

# Keep junctions on the same strand as the underlying gene
jx <- merge(jx, ref, by.x= c('PrimaryGene', 'Site1_strand'), by.y= c('gene_id', 'strand'))

gff <- jx[, list(
    chrom= Site1_chr,
    source= '.',
    feature= 'junction',
    start= Site1_location,
    end= Site2_location,
    v1= '.',
    strand= Site1_strand,
    v2= '.',
    attr= sprintf('gene_id "%s";', PrimaryGene)
)]
stopifnot(nrow(gff) == nrow(unique(gff)))

write.table(x= gff, file= '{output.gff}', col.names= FALSE, row.names= FALSE, quote= FALSE, sep= '\t')

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


rule count_nonsplit_reads:
    input:
        gff= '{species}/featureCounts/junctions.gff',
        bam= '{species}/hisat2/{sample}.pre.bam',
    output:
        cnt= temp('{species}/featureCounts/{sample}.junctions.all_nonsplit'),
        smry= temp('{species}/featureCounts/{sample}.junctions.all_nonsplit.summary'),
    params:
        strand= featureCountsStrand,
    shell:
        r"""
        featureCounts --nonSplitOnly -O -f -a {input.gff} -o {output.cnt} -t junction -g gene_id -s {params.strand} -p -T 2 {input.bam}
        """


rule count_fully_contained:
    input:
        gff= '{species}/featureCounts/junctions.gff',
        bam= '{species}/hisat2/{sample}.pre.bam',
    output:
        cnt= temp('{species}/featureCounts/{sample}.junctions.contained'),
        smry= temp('{species}/featureCounts/{sample}.junctions.contained.summary'),
    params:
        strand= featureCountsStrand,
    shell:
        r"""
        featureCounts --fracOverlap 0.98 --nonSplitOnly -O -f -a {input.gff} -o {output.cnt} -t junction -g gene_id -s {params.strand} -p -T 2 {input.bam}
        """


rule count_intron_supporting_junctions:
    input:
        contained= '{species}/featureCounts/{sample}.junctions.contained',
        full= '{species}/featureCounts/{sample}.junctions.all_nonsplit',
    output:
        cnt= temp('{species}/featureCounts/{sample}.junctions.nonsplit'),
    run:
        with open(input.full) as full, open(input.contained) as contained, open(output.cnt, 'w') as fout:
            for f, c in zip(full, contained):
                if f.startswith('#'):
                    fout.write('#\n')
                elif f.startswith('Geneid\tChr\t'):
                    assert f == c
                    fout.write(f)
                else:
                    f = f.strip().split('\t')
                    c = c.strip().split('\t')
                    assert f[0:6] == c[0:6]
                    assert len(f) == 7 and len(c) == 7
                    out = f
                    cnt = int(f[6]) - int(c[6])
                    assert cnt >= 0
                    f[6] = str(cnt)
                    fout.write('\t'.join(f) + '\n')
        fout.close()


rule merge_nonsplit_counts:
    input:
        cnt= lambda wc: [f'{wc.species}/featureCounts/{x}.junctions.nonsplit' for x in sorted(ss[ss.species == wc.species].sample_id)],
    output:
        cnt= '{species}/featureCounts/junctions.nonsplit',
    run:
        dts = []
        for x in input.cnt:
            dt = pandas.read_csv(x, sep= '\t', comment= '#')
            fid = dt.columns[len(dt.columns)-1]
            dt.rename(columns= {fid: 'count'}, inplace= True)
            sample_id = re.sub('\.junctions\.nonsplit$', '', os.path.basename(x))
            dt['sample_id'] = sample_id
            assert sample_id in fid
            dts.append(dt)

        sites = pandas.concat(dts)
        sites.to_csv(output.cnt, sep= '\t', index= False, na_rep= 'NA')


# rule count_summary:
#     input:
#         nonsplit= expand('{{species}}/hisat2/{sample}.junctions.nonsplit.summary', sample= ss.sample_id),
#         jminus= expand('{{species}}/hisat2/{sample}.counts.minus.summary', sample= ss.sample_id), 
#         jplus= expand('{{species}}/hisat2/{sample}.counts.plus.summary', sample= ss.sample_id), 
#     output:
#         smry= '{species}/hisat2/counts.summary',
#     run:
#         with open(output.smry, 'w') as fout:
#             fout.write('\t'.join(['sample_id', 'junction_type', 'count_of', 'count']) + '\n')
#             for x in input.nonsplit + input.jminus + input.jplus:
#                 fn = os.path.basename(x).split('.')
#                 sample_id, species, jx= [fn[i] for i in [0, 1, 3]]
#                 assert len(ss[(ss.sample_id == sample_id) & (ss.species == species)]) == 1
#                 assert jx in ['nonsplit', 'minus', 'plus']
#                 dat = open(x).readlines()
#                 for x in dat:
#                     if not x.startswith('Status'):
#                         fout.write('\t'.join([sample_id, jx, x.strip()]) + '\n')
# 
#         shell(r"""
# cat <<'EOF' > {rule}.$$.tmp.R 
# 
# library(data.table)
# 
# smry <- fread('{output.smry}')
# tot <- smry[, list(tot= sum(count)), by= list(sample_id, junction_type)]        
# smry <- merge(smry, tot, by= c('sample_id', 'junction_type'))
# smry[, pct := count/tot * 100]
# smry[, junction_type := factor(junction_type, levels= c('plus', 'minus', 'nonsplit'))]
# smry[, sample_id := factor(sample_id, levels= unique(smry[junction_type == 'nonsplit' & count_of == 'Unassigned_NoFeatures'][order(pct)]$sample_id))]
# write.table(file= '{output.smry}', x= smry, sep= '\t', row.names= FALSE, quote= FALSE)
# 
# EOF
# Rscript {rule}.$$.tmp.R
# rm {rule}.$$.tmp.R
#         """)


rule junction_counts:
    input:
        split= '{species}/featureCounts/counts.txt.jcounts',
        nonsplit= '{species}/featureCounts/junctions.nonsplit',
    output:
        cnt= '{species}/edger/junction_counts.tsv.gz',
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

options(scipen= 9) # Prevent genomic coordinates with round numbers to be written with scientific notation

split_cnt <- fread('{input.split}')
split_cnt <- split_cnt[!is.na(PrimaryGene)]
setnames(split_cnt, c('PrimaryGene', 'Site1_chr', 'Site1_location', 'Site2_location'), c('gene_id', 'chrom', 'jx_start', 'jx_end'))
split_cnt[, jx_type := 'split']

nonsplit_cnt <- fread('{input.nonsplit}')
setnames(nonsplit_cnt, c('Geneid', 'Chr', 'Start', 'End'), c('gene_id', 'chrom', 'jx_start', 'jx_end'))
nonsplit_cnt[, jx_type := 'nonsplit']

stopifnot(identical(sort(unique(split_cnt$jx_id)), sort(unique(nonsplit_cnt$jx_id))))

cnt <- rbind(split_cnt[, list(sample_id, chrom, jx_start, jx_end, gene_id, jx_type, count)], 
          nonsplit_cnt[, list(sample_id, chrom, jx_start, jx_end, gene_id, jx_type, count)])
cnt[, jx_start := jx_start - 1] # 0-based coord

ocount <- dcast.data.table(data= cnt, sample_id + gene_id + chrom + jx_start + jx_end ~ jx_type, value.var= 'count', fill= 0) 
gz <- gzfile('{output.cnt}', 'w')
write.table(file= gz, x= ocount, sep= '\t', row.names= FALSE, quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

rule differential_junctions:
    input:
        contrasts= config['contrasts'], 
        cnt= '{species}/edger/junction_counts.tsv.gz',
    output:
        bed= temp('{species}/edger/differential_junctions.{contrast}.bed.gz'),
    params:
        min_count= config['min_count'],
        min_pct_pass= config['min_pct_pass'],
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R

library(data.table)
library(edgeR)

options(scipen= 9) # Prevent genomic coordinates with round numbers to be written with scientific notation

contrasts <- fread(cmd= 'grep -v "^#" {input.contrasts}')
contrasts <- contrasts[species == '{wildcards.species}' & contrast == '{wildcards.contrast}']
contrasts[, species := NULL]
contrasts[, contrast := NULL]

cnt <- fread('{input.cnt}')
cnt <- melt(data= cnt, id.vars= c('sample_id', 'gene_id', 'chrom', 'jx_start', 'jx_end'), variable.name= 'jx_type', value.name= 'count')

stopifnot(contrasts$sample_id %in% cnt$sample_id)

cnt <- cnt[sample_id %in% contrasts$sample_id]
cnt[, jx_id := paste(gene_id, chrom, jx_start, jx_end, sep= '::')]
cnt <- cnt[, list(sample_id, jx_id, jx_type, count)]

cnt <- merge(cnt, contrasts[, list(sample_id, side)], by= 'sample_id', sort= FALSE)

# Filter sites with low counts. We are not intersted in differential
# expression so want a decent number of counts in both lhs and rhs
tot <- cnt[, list(count= sum(count)), by= list(sample_id, jx_id, side)]
tot <- tot[, list(n_libs= .N, n_libs_pass= sum(count > {params.min_count})), by= list(jx_id, side)]
tot[, pct := n_libs_pass/n_libs]
keep <- tot[, list(n_cmp_pass= sum(pct >= {params.min_pct_pass})), by= jx_id][n_cmp_pass >= 2]
cnt <- cnt[jx_id %in% keep$jx_id]

cat(sprintf('%s %s\n', '{wildcards.contrast}', length(unique(cnt$jx_id))))

# DGE #
# Using strategy for BSseq data. See
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5747346/
mat <- dcast.data.table(data= cnt, jx_id ~ sample_id + side + jx_type, value.var= 'count', sep= '::', fill= 0)
mat <- as.matrix(mat, rownames= 'jx_id')

sample_id <- as.factor(sapply(strsplit(colnames(mat), split= '::'), function(x) x[[1]]))
side <- as.factor(sapply(strsplit(colnames(mat), split= '::'), function(x) x[[2]]))
jx_type <- as.factor(sapply(strsplit(colnames(mat), split= '::'), function(x) x[[3]]))

design.samples <- model.matrix(~ 0 + sample_id)
colnames(design.samples) <- gsub('-', '_', colnames(design.samples))
colnames(design.samples) <- sub('sample_id', '', colnames(design.samples))

design.side <- model.matrix(~ 0 + side)
colnames(design.side) <- sub('side', '', colnames(design.side))
design  <- cbind(design.samples, (jx_type == "split") * design.side)
colnames(design) <- ifelse(grepl('^\\d', colnames(design)), paste0('x', colnames(design)), colnames(design))
contr <- makeContrasts(contr= lhs - rhs, levels= design)

libsize <- cnt[, list(libsize= sum(count)), by= sample_id]
size <- libsize$libsize
names(size) <- libsize$sample_id
sizes <- rep(size, each= 2)
sizes <- sizes[match(names(sizes), sub('::.*', '', colnames(mat)))]
stopifnot(identical(sub('::.*', '', colnames(mat)), names(sizes)))

y <- DGEList(mat, lib.size= sizes)
y <- estimateDisp(y, design= design)

## NB: glmQLFit is faaaar more conservative than glmFit
fit <- glmFit(y, design, robust= TRUE)
lrt <- glmTreat(fit, contrast= contr, lfc= log2(1.5))

detable <- as.data.table(topTags(lrt, sort.by= 'none', n= Inf)$table, keep.rownames= 'jx_id')
detable[, unshrunk.logFC := NULL]
detable[, contrast := '{wildcards.contrast}']

side_avg <- dcast.data.table(data= cnt[, list(tot= sum(count)), by= list(jx_id, jx_type, side)], jx_id ~ side + jx_type, value.var= 'tot', fill= 0)
side_avg <- side_avg[, list(jx_id, lhs= 100 * lhs_split/(lhs_nonsplit + lhs_split + 0.125), 
                                   rhs= 100 * rhs_split/(rhs_nonsplit + rhs_split + 0.125))]
side_avg[, lhs_rhs_diff := lhs - rhs]

detable <- merge(detable, side_avg, by= 'jx_id', sort= FALSE)
detable[, c('gene_id', 'chrom', 'jx_start', 'jx_end') := tstrsplit(detable$jx_id, '::', fixed= TRUE)]
detable[, jx_start := as.numeric(jx_start)]
detable[, jx_end := as.numeric(jx_end)]

detable[, jx_id := NULL]
detable[, species := '{wildcards.species}']

gz <- gzfile('{output.bed}', 'w')
write.table(detable, file= gz, row.names= FALSE, sep= '\t', quote= FALSE)
close(gz)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """

def collect_differential_junctions(contrasts_tsv, species):
    cntr = pandas.read_csv(contrasts_tsv, sep= '\t', comment= '#')
    cntr = cntr[cntr.species == species]
    jx = []
    for contrast in cntr.contrast.unique():
        jx.append(f'{species}/edger/differential_junctions.{contrast}.bed.gz')
    return jx

rule cat_differential_junctions:
    input:
        contrasts= config['contrasts'], 
        gff= 'ref/{species}.gff',
        bed= lambda wc: collect_differential_junctions(config['contrasts'], wc.species),
    output:
        tsv= '{species}/edger/differential_junctions.tsv.gz',
    run:
        sys.path.append(workflow.basedir)
        
        from scripts import utils
        from urllib.parse import unquote

        tsv = []
        for x in input.bed:
            tsv.append(pandas.read_csv(x, sep= '\t'))
        tsv = pandas.concat(tsv)

        gff = pandas.read_csv(input.gff, sep= '\t', header= None)
        gff = gff[gff[2] == 'gene'].copy()
        gff['gene_id'] = [utils.gff_attribute(x, 'ID') for x in gff[8]]
        assert(len(gff[gff.gene_id == '']) == 0)
        gff['description'] = [utils.gff_attribute(x, 'description') for x in gff[8]]
        gff['description'] = [unquote(x) for x in gff.description]
        gff.rename(columns= {6: 'strand'}, inplace= True) 
        gff = gff[['gene_id', 'strand', 'description']].drop_duplicates()
        
        tsv = tsv.merge(gff, on= 'gene_id')
        
        col_order = ['chrom', 'jx_start', 'jx_end', 'gene_id', 'strand', 'contrast', 'lhs', 'rhs', 'lhs_rhs_diff', 'logFC', 'logCPM', 'PValue', 'FDR', 'species', 'description']
        assert(sorted(col_order) == sorted(tsv.columns.tolist()))
        tsv = tsv[col_order]

        tsv.to_csv(output.tsv, compression= 'gzip', sep= '\t', index= False)
