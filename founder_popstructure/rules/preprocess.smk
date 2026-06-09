# workflow/rules/preprocess.smk

rule angsd_bamlist:
    input:
        expand("mapping/final/{sample}.bam", sample=config["global"]["sample-names"])
    output:
        "angsd/individuals.bamlist"
    log:
        "logs/angsd/bamlist.log"
    run:
        os.makedirs("angsd", exist_ok=True)
        with open(output[0], "w") as out:
            for bam in input:
                out.write(str(bam) + "\n")

rule generate_popmap:
    input:
        "angsd/individuals.bamlist"
    output:
        "angsd/popmap.txt"
    log:
        "logs/angsd/popmap.log"
    shell:
        """
        (while read p; do
            name=$(basename $p .bam)
            echo -e "${{name}}\tColony"
        done < {input} > {output}) 2> {log}
        """

rule generate_genome_mask:
    input:
        ref = config["data"]["reference-genome"],
        repeats = REPEAT_GFF
    output:
        bed = "refs/genome_mask.bed.gz"
    log:
        "logs/refs/genome_mask.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        r"""
        (seqkit locate -p "[Nn]+" {input.ref} --bed > refs/gaps.tmp.bed
        awk '!/^#/ {{print $1, $4-1, $5}}' {input.repeats} > refs/repeats.tmp.bed
        cat refs/gaps.tmp.bed refs/repeats.tmp.bed \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - \
            | gzip > {output.bed}
        rm refs/gaps.tmp.bed refs/repeats.tmp.bed) &> {log}
        """