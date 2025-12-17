rule fastqc:
    input:
        unpack(get_reads_to_trim)
    output:
        html_R1 = RES_DIR_RULE + "fastqc/{fastq_id}_R1_fastqc.html",
        zip_R1 = RES_DIR_RULE + "fastqc/{fastq_id}_R1_fastqc.zip",
        html_R2 = RES_DIR_RULE + "fastqc/{fastq_id}_R2_fastqc.html",
        zip_R2 = RES_DIR_RULE + "fastqc/{fastq_id}_R2_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    log:
        out = LOG_DIR_RULE + "fastqc/{fastq_id}.out",
        err = LOG_DIR_RULE + "fastqc/{fastq_id}.err"
    params:
        out_dir = RES_DIR_RULE + "fastqc"
    shell:
        # If FASTQ base name does not match {fastq_id}, then output files won't have the correct name and the pipeline will stall forever
        # To avoid this, stream data from stdin so that output file name can be customized to get the correct name
        """
        gunzip --stdout {input.fq1} \
        | \
        fastqc --outdir {params.out_dir} \
               stdin:{wildcards.fastq_id}_R1 \
               1>> {log.out} \
               2>> {log.err}

        gunzip --stdout {input.fq2} \
        | \
        fastqc --outdir {params.out_dir} \
               stdin:{wildcards.fastq_id}_R2 \
               1>> {log.out} \
               2>> {log.err}
        """

rule multiqc:
    input:
        zip_R1 = expand(RES_DIR_RULE + "fastqc/{fastq_id}_R1_fastqc.zip", fastq_id = fastq_ids),
        zip_R2 = expand(RES_DIR_RULE + "fastqc/{fastq_id}_R2_fastqc.zip", fastq_id = fastq_ids)
    output:
        html = RES_DIR_RULE + "multiqc/multiqc_report.html",
        data = directory(RES_DIR_RULE + "multiqc/multiqc_report_data")
    conda:
        "../envs/fastqc.yml"
    log:
        out = LOG_DIR_RULE + "multiqc/multiqc.out",
        err = LOG_DIR_RULE + "multiqc/multiqc.err"
    params:
        out_dir = RES_DIR_RULE + "multiqc"
    shell:
        # Let multiQC search the fastqc.zip files and generate an html report based on any log files that it recognises
        # --filename: name the html report file with the same name as {output.html}
        # --data-dir: along with the html report, create a parsed data directory "multiqc_data/" containing tab-delimited data files with additional information
        # --verbose: get more details of what MultiQC does (e.g., overwritting files)
        """
        multiqc {input.zip_R1} \
                {input.zip_R2} \
                --outdir {params.out_dir} \
                --filename "multiqc_report.html" \
                --data-dir \
                --verbose \
                1>> {log.out} \
                2>> {log.err}
        """

rule trimmomatic:
    input:
        unpack(get_reads_to_trim)
    output:
        r1p = RES_DIR_RULE + "trimmomatic/{fastq_id}_1P.fastq.gz",
        r1u = RES_DIR_RULE + "trimmomatic/{fastq_id}_1U.fastq.gz",
        r2p = RES_DIR_RULE + "trimmomatic/{fastq_id}_2P.fastq.gz",
        r2u = RES_DIR_RULE + "trimmomatic/{fastq_id}_2U.fastq.gz",
        sum = RES_DIR_RULE + "trimmomatic/{fastq_id}_summary.txt"
    conda:
        "../envs/trimmomatic.yml"
    log:
        out = LOG_DIR_RULE + "trimmomatic/{fastq_id}.out",
        err = LOG_DIR_RULE + "trimmomatic/{fastq_id}.err"
    params:
        out_dir = RES_DIR_RULE + "trimmomatic"
    shell:
        """
        # If output directory does not exist, create it
        if [[ -d {params.out_dir} ]]
        then
            echo "Ok, output directory exists"
        else
            mkdir -p {params.out_dir}
        fi;

        # Trim adapters!
        trimmomatic PE \
                    -threads {threads} \
                    {input.fq1} \
                    {input.fq2} \
                    -baseout "{params.out_dir}/{wildcards.fastq_id}.fastq.gz" \
                    {TRIMMOMATIC_SETTING} \
                    -summary {output.sum} \
                    1>> {log.out} \
                    2>> {log.err}
        """

rule fastqc_post:
    input:
        r1p = RES_DIR_RULE + "trimmomatic/{fastq_id}_1P.fastq.gz",
        r2p = RES_DIR_RULE + "trimmomatic/{fastq_id}_2P.fastq.gz"
    output:
        html_R1 = RES_DIR_RULE + "fastqc_post/{fastq_id}_1P_fastqc.html",
        zip_R1 = RES_DIR_RULE + "fastqc_post/{fastq_id}_1P_fastqc.zip",
        html_R2 = RES_DIR_RULE + "fastqc_post/{fastq_id}_2P_fastqc.html",
        zip_R2 = RES_DIR_RULE + "fastqc_post/{fastq_id}_2P_fastqc.zip"
    conda:
        "../envs/fastqc.yml"
    log:
        out = LOG_DIR_RULE + "fastqc_post/{fastq_id}.out",
        err = LOG_DIR_RULE + "fastqc_post/{fastq_id}.err"
    params:
        out_dir = RES_DIR_RULE + "fastqc_post"
    shell:
        """
        fastqc --outdir {params.out_dir} \
               {input.r1p} \
               1>> {log.out} \
               2>> {log.err}

        fastqc --outdir {params.out_dir} \
               {input.r2p} \
               1>> {log.out} \
               2>> {log.err}
        """

rule multiqc_post:
    input:
        zip_R1 = expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_1P_fastqc.zip", fastq_id = fastq_ids),
        zip_R2 = expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_2P_fastqc.zip", fastq_id = fastq_ids)
    output:
        html = RES_DIR_RULE + "multiqc_post/multiqc_report.html",
        data = directory(RES_DIR_RULE + "multiqc_post/multiqc_report_data")
    conda:
        "../envs/fastqc.yml"
    log:
        out = LOG_DIR_RULE + "multiqc_post/multiqc.out",
        err = LOG_DIR_RULE + "multiqc_post/multiqc.err"
    params:
        out_dir = RES_DIR_RULE + "multiqc_post"
    shell:
        # Let multiQC search the fastqc.zip files and generate an html report based on any log files that it recognises
        # --filename: name the html report file with the same name as {output.html}
        # --data-dir: along with the html report, create a parsed data directory "multiqc_data/" containing tab-delimited data files with additional information
        # --verbose: get more details of what MultiQC does (e.g., overwritting files)
        """
        multiqc {input.zip_R1} \
                {input.zip_R2} \
                --outdir {params.out_dir} \
                --filename "multiqc_report.html" \
                --data-dir \
                --verbose \
                1>> {log.out} \
                2>> {log.err}
        """
