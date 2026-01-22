# controlseq
#-----------

rule split_intervals_controlseq:
    input:
        dict = expand(RES_DIR_RULE + "samtools_dict_faidx/{refGenome}.dict", refGenome={REFGENOME}),
        fai = expand(RES_DIR_RULE + "samtools_dict_faidx/{refGenome}.fasta.fai", refGenome={REFGENOME})
    output:
        # The "-scattered" part of outfile name is added by SplitIntervals and cannot be changed
        intervals = expand(RES_DIR_RULE + "split_intervals_controlseq/{interval_controlseq}-scattered.interval_list", interval_controlseq = INTERVALS_CONTROLSEQ)
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "split_intervals_controlseq/split_intervals_controlseq.out",
        err = LOG_DIR_RULE + "split_intervals_controlseq/split_intervals_controlseq.err"
    params:
        inp_dir = RES_DIR_RULE + "samtools_dict_faidx",
        out_dir = RES_DIR_RULE + "split_intervals_controlseq"
    shell:
        """
        gatk SplitIntervals --reference {params.inp_dir}/{REFGENOME}.fasta \
                            --intervals {CONTROLSEQ} \
                            --scatter-count {N_INTERVALS_CONTROLSEQ} \
                            --interval-file-num-digits {N_DIGITS_CONTROLSEQ} \
                            --output {params.out_dir} \
                            1>> {log.out} \
                            2>> {log.err}
        """


# targetseq
#----------

rule split_intervals_targetseq:
    input:
        dict = expand(RES_DIR_RULE + "samtools_dict_faidx/{refGenome}.dict", refGenome={REFGENOME}),
        fai = expand(RES_DIR_RULE + "samtools_dict_faidx/{refGenome}.fasta.fai", refGenome={REFGENOME})
    output:
        # The "-scattered" part of outfile name is added by SplitIntervals and cannot be changed
        intervals = expand(RES_DIR_RULE + "split_intervals_targetseq/{interval_targetseq}-scattered.interval_list", interval_targetseq = INTERVALS_TARGETSEQ)
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "split_intervals_targetseq/split_intervals_targetseq.out",
        err = LOG_DIR_RULE + "split_intervals_targetseq/split_intervals_targetseq.err"
    params:
        inp_dir = RES_DIR_RULE + "samtools_dict_faidx",
        out_dir = RES_DIR_RULE + "split_intervals_targetseq"
    shell:
        """
        gatk SplitIntervals --reference {params.inp_dir}/{REFGENOME}.fasta \
                            --intervals {TARGETSEQ} \
                            --scatter-count {N_INTERVALS_TARGETSEQ} \
                            --interval-file-num-digits {N_DIGITS_TARGETSEQ} \
                            --output {params.out_dir} \
                            1>> {log.out} \
                            2>> {log.err}
        """

rule bcftools_call_allsites_controllib_targetseq:
    input:
        bam = get_controllib_bam_to_call_allsites_on_targetseq(),
        intervals = RES_DIR_RULE + "split_intervals_targetseq/{interval_targetseq}-scattered.interval_list"
    output:
        vcf = RES_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.raw.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.err"
    params:
        temp_interval = RES_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.txt",
        sample_ids_str = lambda wildcards: "\n".join(sample_ids_controllib),
        temp_samples = RES_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.samples.txt"
    shell:
        # Save in temporary file:
        #     - interval file after removing the header (lines starting with "@") as not supported by bcftools mpileup
        #     - sample IDs separated by newlines as required by bcftools reheader
        """
        awk 'BEGIN {{OFS="\\t"}};
             !/^@/ {{print $1, $2, $3}}' {input.intervals} \
            > {params.temp_interval}

        echo "{params.sample_ids_str}" > {params.temp_samples}

        bcftools mpileup --threads {threads} \
                         --fasta-ref {REFPATH} \
                         --max-depth 2147483647 \
                         -a FORMAT/DP \
                         --regions-file {params.temp_interval} \
                         {input.bam} \
                         -Ou \
                         2>> {log.err} \
        | \
        bcftools call --threads {threads} \
                      --multiallelic-caller \
                      --format-fields GQ \
                      --output-type z \
                      2>> {log.err} \
        | \
        bcftools reheader --threads {threads} \
                          --samples {params.temp_samples} \
                          --output {output.vcf} \
                          1>> {log.out} \
                          2>> {log.err}

        rm {params.temp_interval}
        rm {params.temp_samples}
        """

rule bcftools_call_allsites_targetlib_targetseq:
    input:
        bam = get_targetlib_bam_to_call_allsites_on_targetseq(),
        intervals = RES_DIR_RULE + "split_intervals_targetseq/{interval_targetseq}-scattered.interval_list"
    output:
        vcf = RES_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.raw.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.err"
    params:
        temp_interval = RES_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.txt",
        sample_ids_str = lambda wildcards: "\n".join(sample_ids_targetlib),
        temp_samples = RES_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.samples.txt"
    shell:
        # Save in temporary file:
        #     - interval file after removing the header (lines starting with "@") as not supported by bcftools mpileup
        #     - sample IDs separated by newlines as required by bcftools reheader
        """
        awk 'BEGIN {{OFS="\\t"}};
             !/^@/ {{print $1, $2, $3}}' {input.intervals} \
            > {params.temp_interval}

        echo "{params.sample_ids_str}" > {params.temp_samples}

        bcftools mpileup --threads {threads} \
                         --fasta-ref {REFPATH} \
                         --max-depth 2147483647 \
                         -a FORMAT/DP \
                         --regions-file {params.temp_interval} \
                         {input.bam} \
                         -Ou \
                         2>> {log.err} \
        | \
        bcftools call --threads {threads} \
                      --multiallelic-caller \
                      --format-fields GQ \
                      --output-type z \
                      2>> {log.err} \
        | \
        bcftools reheader --threads {threads} \
                          --samples {params.temp_samples} \
                          --output {output.vcf} \
                          1>> {log.out} \
                          2>> {log.err}

        rm {params.temp_interval}
        rm {params.temp_samples}
        """
