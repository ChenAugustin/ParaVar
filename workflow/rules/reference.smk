rule hisat2_index:
    output:
        indexes = expand(RES_DIR_RULE + "hisat2_index/{{REFGENOME}}.{index}.ht2", index=[1, 2, 3, 4, 5, 6, 7, 8])
    conda:
        "../envs/hisat2.yml"
    log:
        out = LOG_DIR_RULE + "hisat2_index/{REFGENOME}.out",
        err = LOG_DIR_RULE + "hisat2_index/{REFGENOME}.err"
    params:
        out_dir = RES_DIR_RULE + "hisat2_index"
    shell:
        """
        hisat2-build -p {threads} \
                     -f {REFPATH} \
                     {params.out_dir}/{REFGENOME} \
                     1>> {log.out} \
                     2>> {log.err}
        """

rule bwa_mem2_index:
    output:
        indexes = expand(RES_DIR_RULE + "bwa_mem2_index/{{REFGENOME}}.fasta.{index}", index=["0123", "amb", "ann", "bwt.2bit.64", "pac"])
    conda:
        "../envs/bwa_mem2.yml"
    log:
        out = LOG_DIR_RULE + "bwa_mem2_index/{REFGENOME}.out",
        err = LOG_DIR_RULE + "bwa_mem2_index/{REFGENOME}.err"
    params:
        out_dir = RES_DIR_RULE + "bwa_mem2_index"
    shell:
        # To avoid copying the reference FASTA, create a symbolic link
        # And index the symbolic link to the reference FASTA
        """
        ln --symbolic --no-target-directory {REFPATH} {params.out_dir}/{REFGENOME}.fasta
        bwa-mem2 index {params.out_dir}/{REFGENOME}.fasta \
                       1>> {log.out} \
                       2>> {log.err}
        """

rule minimap2_index:
    output:
        indexes = expand(RES_DIR_RULE + "minimap2_index/{{REFGENOME}}.{index}", index=["mmi"])
    conda:
        "../envs/minimap2.yml"
    log:
        out = LOG_DIR_RULE + "minimap2_index/{REFGENOME}.out",
        err = LOG_DIR_RULE + "minimap2_index/{REFGENOME}.err"
    shell:
        """
        minimap2 -d {output.indexes} \
                 {REFPATH}
        """

rule samtools_dict_faidx:
    output:
        dict = RES_DIR_RULE + "samtools_dict_faidx/{REFGENOME}.dict",
        fai = RES_DIR_RULE + "samtools_dict_faidx/{REFGENOME}.fasta.fai"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "samtools_dict_faidx/{REFGENOME}.out",
        err = LOG_DIR_RULE + "samtools_dict_faidx/{REFGENOME}.err"
    params:
        out_dir = RES_DIR_RULE + "samtools_dict_faidx"
    shell:
        # To avoid copying the reference FASTA, create a symbolic link
        # And index the symbolic link to the reference FASTA
        """
        ln --symbolic --no-target-directory {REFPATH} {params.out_dir}/{REFGENOME}.fasta
        samtools dict {params.out_dir}/{REFGENOME}.fasta \
                      --output {output.dict} \
                      1>> {log.out} \
                      2>> {log.err}
        samtools faidx {params.out_dir}/{REFGENOME}.fasta \
                       1>> {log.out} \
                       2>> {log.err}
        """
