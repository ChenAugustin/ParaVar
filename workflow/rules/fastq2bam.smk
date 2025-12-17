rule hisat2:
    input:
        unpack(get_reads_to_map),
        indexes = expand(RES_DIR_RULE + "hisat2_index/{refGenome}.{index}.ht2", refGenome={REFGENOME}, index=[1, 2, 3, 4, 5, 6, 7, 8])
    output:
        bam = RES_DIR_RULE + "hisat2/{sample_id}/{sample_id}.bam",
        bai = RES_DIR_RULE + "hisat2/{sample_id}/{sample_id}.bam.bai",
        sum = RES_DIR_RULE + "hisat2/{sample_id}/{sample_id}.flagstat.txt"
    conda:
        "../envs/hisat2.yml"
    log:
        out = LOG_DIR_RULE + "hisat2/{sample_id}/{sample_id}.out",
        err = LOG_DIR_RULE + "hisat2/{sample_id}/{sample_id}.err"
    params:
        out_dir = RES_DIR_RULE + "hisat2/{sample_id}",
        idx_dir = RES_DIR_RULE + "hisat2_index"
    shell:
        # If >1 r1 file, then concat them before mapping (same for r2)
        # Output is in SAM format by default so, use samtools view to convert it to bam format
        # Use same number of threads for all tools. Note that samtool's multithreading is incompatible with piping it to sambamba sort
        # Use touch to create the INT_DIR_STRUCT directory, otherwise sambamba throws "no such file or directory"
        # Index sorted BAM file required for variant calling
        # Get some performace statistics with samtools flagstat
        """
        n_r1=$(echo "{input.r1}" | awk '{{print NF}}')
        if [ "$n_r1" -gt 1 ];
        then
            zcat {input.r1} > {params.out_dir}/{wildcards.sample_id}_1P.fastq;
            gzip {params.out_dir}/{wildcards.sample_id}_1P.fastq;
            zcat {input.r2} > {params.out_dir}/{wildcards.sample_id}_2P.fastq;
            gzip {params.out_dir}/{wildcards.sample_id}_2P.fastq;

            hisat2_input_r1="{params.out_dir}/{wildcards.sample_id}_1P.fastq.gz";
            hisat2_input_r2="{params.out_dir}/{wildcards.sample_id}_2P.fastq.gz";
        else
            hisat2_input_r1="{input.r1}";
            hisat2_input_r2="{input.r2}";
        fi

        touch {output.bam}

        hisat2 --threads {threads} \
               --time \
               --un-gz {params.out_dir}/{wildcards.sample_id}.un.fastq.gz \
               --al-gz {params.out_dir}/{wildcards.sample_id}.al.fastq.gz \
               --un-conc-gz {params.out_dir}/{wildcards.sample_id}.un_conc.fastq.gz \
               --al-conc-gz {params.out_dir}/{wildcards.sample_id}.al_conc.fastq.gz \
               --summary-file {params.out_dir}/{wildcards.sample_id}.summary.txt \
               --new-summary \
               --met-file {params.out_dir}/{wildcards.sample_id}.met.txt \
               --rg-id {wildcards.sample_id} \
               --no-temp-splicesite \
               --no-spliced-alignment \
               -q \
               -x {params.idx_dir}/{REFGENOME} \
               -1 $hisat2_input_r1 \
               -2 $hisat2_input_r2 \
               2>> {log.err} \
        | \
        samtools view /dev/stdin \
                      -b \
                      -h \
                      2>> {log.err} \
        | \
        sambamba sort /dev/stdin \
                      --nthreads {threads} \
                      --memory-limit "{resources.mem_mb}MB" \
                      --tmpdir {params.out_dir}/{wildcards.sample_id} \
                      --out {output.bam} \
                      1>> {log.out} \
                      2>> {log.err}

        samtools index -@ {threads} \
                       {output.bam}

        samtools flagstat -@ {threads} {output.bam} > {output.sum}
        """

rule bwa_mem2:
    input:
        unpack(get_reads_to_map),
        indexes = ancient(expand(RES_DIR_RULE + "bwa_mem2_index/{refGenome}.fasta.{index}", refGenome={REFGENOME}, index=["0123", "amb", "ann", "bwt.2bit.64", "pac"]))
    output:
        bam = RES_DIR_RULE + "bwa_mem2/{sample_id}/{sample_id}.bam",
        bai = RES_DIR_RULE + "bwa_mem2/{sample_id}/{sample_id}.bam.bai",
        sum = RES_DIR_RULE + "bwa_mem2/{sample_id}/{sample_id}.flagstat.txt"
    conda:
        "../envs/bwa_mem2.yml"
    log:
        out = LOG_DIR_RULE + "bwa_mem2/{sample_id}/{sample_id}.out",
        err = LOG_DIR_RULE + "bwa_mem2/{sample_id}/{sample_id}.err"
    params:
        out_dir = RES_DIR_RULE + "bwa_mem2/{sample_id}",
        idx_dir = RES_DIR_RULE + "bwa_mem2_index"
    shell:
        # If >1 r1 file, then concat them before mapping (same for r2)
        # Rune BWA-MEM2 using parameters recommended by Regier et al. 2018 (10.1038/s41467-018-06159-4) for reproducibility
        """
        n_r1=$(echo "{input.r1}" | awk '{{print NF}}')
        if [ "$n_r1" -gt 1 ];
        then
            zcat {input.r1} > {params.out_dir}/{wildcards.sample_id}_1P.fastq;
            gzip {params.out_dir}/{wildcards.sample_id}_1P.fastq;
            zcat {input.r2} > {params.out_dir}/{wildcards.sample_id}_2P.fastq;
            gzip {params.out_dir}/{wildcards.sample_id}_2P.fastq;

            bwa_mem2_input_r1="{params.out_dir}/{wildcards.sample_id}_1P.fastq.gz";
            bwa_mem2_input_r2="{params.out_dir}/{wildcards.sample_id}_2P.fastq.gz";
        else
            bwa_mem2_input_r1="{input.r1}";
            bwa_mem2_input_r2="{input.r2}";
        fi

        bwa-mem2 mem -t {threads} \
                     -K 100000000 \
                     -Y \
                     -R "@RG\\tID:{wildcards.sample_id}" \
                     {params.idx_dir}/{REFGENOME}.fasta \
                     $bwa_mem2_input_r1 \
                     $bwa_mem2_input_r2 \
                     2>> {log.err} \
        | \
        samtools view /dev/stdin \
                      -b \
                      -h \
                      2>> {log.err} \
        | \
        sambamba sort /dev/stdin \
                      --nthreads {threads} \
                      --memory-limit "{resources.mem_mb}MB" \
                      --tmpdir {params.out_dir}/{wildcards.sample_id} \
                      --out {output.bam} \
                      1>> {log.out} \
                      2>> {log.err}

        samtools index -@ {threads} \
                       {output.bam}

        samtools flagstat -@ {threads} {output.bam} > {output.sum}
        """

if config["read_type"] == "long_hifi":

    rule minimap2:
        input:
            unpack(get_reads_to_map),
            indexes = ancient(expand(RES_DIR_RULE + "minimap2_index/{refGenome}.{index}", refGenome={REFGENOME}, index=["mmi"]))
        output:
            bam = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.bam",
            bai = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.bam.bai",
            sum = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.flagstat.txt"
        conda:
            "../envs/minimap2.yml"
        log:
            out = LOG_DIR_RULE + "minimap2/{sample_id}/{sample_id}.out",
            err = LOG_DIR_RULE + "minimap2/{sample_id}/{sample_id}.err"
        params:
            out_dir = RES_DIR_RULE + "minimap2"
        shell:
            # If >1 r1 file, then concat them before mapping (same for r2, if present)
            # -a: Require output in SAM format instead of the default PAF format for variant calling
            # --MD: Get the MD tag, potentially relevant for downstream analyses
            """
            n_r1=$(echo "{input.r1}" | awk '{{print NF}}')
            if [ "$n_r1" -gt 1 ];
            then
                zcat {input.r1} > {params.out_dir}/{wildcards.sample_id}_1P.fastq;
                gzip {params.out_dir}/{wildcards.sample_id}_1P.fastq;

                minimap2_input_r1="{params.out_dir}/{wildcards.sample_id}_1P.fastq.gz";
            else
                minimap2_input_r1="{input.r1}";
            fi
    
            minimap2 -t {threads} \
                     -a \
                     -x "map-hifi" \
                     --MD \
                     {input.indexes} \
                     $minimap2_input_r1 \
                     2>> {log.err} \
            | \
            samtools view /dev/stdin \
                          -b \
                          -h \
                          2>> {log.err} \
            | \
            sambamba sort /dev/stdin \
                          --nthreads {threads} \
                          --memory-limit "{resources.mem_mb}MB" \
                          --tmpdir {params.out_dir}/{wildcards.sample_id} \
                          --out {output.bam} \
                          1>> {log.out} \
                          2>> {log.err}
    
            samtools index -@ {threads} \
                           {output.bam}
    
            samtools flagstat -@ {threads} {output.bam} > {output.sum}
            """

elif config["read_type"] == "short_pe":

    rule minimap2:
        input:
            unpack(get_reads_to_map),
            indexes = ancient(expand(RES_DIR_RULE + "minimap2_index/{refGenome}.{index}", refGenome={REFGENOME}, index=["mmi"]))
        output:
            bam = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.bam",
            bai = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.bam.bai",
            sum = RES_DIR_RULE + "minimap2/{sample_id}/{sample_id}.flagstat.txt"
        conda:
            "../envs/minimap2.yml"
        log:
            out = LOG_DIR_RULE + "minimap2/{sample_id}/{sample_id}.out",
            err = LOG_DIR_RULE + "minimap2/{sample_id}/{sample_id}.err"
        params:
            out_dir = RES_DIR_RULE + "minimap2"
        shell:
            # If >1 r1 file, then concat them before mapping (same for r2, if present)
            # -a: Require output in SAM format instead of the default PAF format for variant calling
            # --MD: Get the MD tag, potentially relevant for downstream analyses
            """
            n_r1=$(echo "{input.r1}" | awk '{{print NF}}')
            if [ "$n_r1" -gt 1 ];
            then
                zcat {input.r1} > {params.out_dir}/{wildcards.sample_id}_1P.fastq;
                gzip {params.out_dir}/{wildcards.sample_id}_1P.fastq;
                zcat {input.r2} > {params.out_dir}/{wildcards.sample_id}_2P.fastq;
                gzip {params.out_dir}/{wildcards.sample_id}_2P.fastq;
    
                minimap2_input_r1="{params.out_dir}/{wildcards.sample_id}_1P.fastq.gz";
                minimap2_input_r2="{params.out_dir}/{wildcards.sample_id}_2P.fastq.gz";
            else
                minimap2_input_r1="{input.r1}";
                minimap2_input_r2="{input.r2}";
            fi
    
            minimap2 -t {threads} \
                     -a \
                     -x "sr" \
                     --MD \
                     {input.indexes} \
                     $minimap2_input_r1 \
                     $minimap2_input_r2 \
                     2>> {log.err} \
            | \
            samtools view /dev/stdin \
                          -b \
                          -h \
                          2>> {log.err} \
            | \
            sambamba sort /dev/stdin \
                          --nthreads {threads} \
                          --memory-limit "{resources.mem_mb}MB" \
                          --tmpdir {params.out_dir}/{wildcards.sample_id} \
                          --out {output.bam} \
                          1>> {log.out} \
                          2>> {log.err}
    
            samtools index -@ {threads} \
                           {output.bam}
    
            samtools flagstat -@ {threads} {output.bam} > {output.sum}
            """
