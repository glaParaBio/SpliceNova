rule download_fasta_reference:
    output:
        genome= 'ref/{species}.fasta',
        fai= 'ref/{species}.fasta.fai',
    params:
        url= lambda wc: species[species.species== wc.species].fasta.iloc[0],
    shell:
        r"""
        curl -L {params.url} > {output.genome}
        samtools faidx {output.genome}
        """


rule download_gff_reference:
    output:
        gff= 'ref/{species}.gff',
    params:
        url= lambda wc: species[species.species== wc.species].gff.iloc[0],
    shell:
        r"""
        curl  -L {params.url} \
        | awk -v FS='\t' -v OFS='\t' '$1 !~ "^#" {{
                                        biotype = ";biotype="$3;
                                        if($3 == "protein_coding_gene") {{$3 = "gene"}}
                                        if($3 == "ncRNA_gene") {{$3 = "gene"}}
                                        if($3 == "pseudogene") {{$3 = "gene"}}
                                        print $0 biotype}}' > {output.gff}
        """


rule hisat_index:
    input:
        species= 'ref/{species}.fasta',
    output:
        idx= 'ref/{species}.8.ht2',
    shell:
        r"""
        hisat2-build -p 4 --seed 1234 -f {input.species} ref/{wildcards.species}
        """

