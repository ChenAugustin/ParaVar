rule awk_remove_indels_controllib_targetseq:
    """
    Remove indels from controllib VCF file as problematic for downstream steps (e.g., creates duplicates in .gdepth).
    """
    input:
        vcf = RES_DIR_RULE + "bcftools_call_allsites_controllib_targetseq/{interval_targetseq}.raw.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "awk_remove_indels_controllib_targetseq/{interval_targetseq}.nindel.vcf.gz",
        tbi = RES_DIR_RULE + "awk_remove_indels_controllib_targetseq/{interval_targetseq}.nindel.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "awk_remove_indels_controllib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "awk_remove_indels_controllib_targetseq/{interval_targetseq}.err"
    shell:
        """
        zcat {input.vcf} \
             2>> {log.err} \
        | \
        awk '$0 !~/INDEL/ {{print}}' \
            2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule awk_remove_indels_targetlib_targetseq:
    """
    Remove indels from targetlib VCF file as problematic for downstream steps (e.g., creates duplicates in .gdepth).
    """
    input:
        vcf = RES_DIR_RULE + "bcftools_call_allsites_targetlib_targetseq/{interval_targetseq}.raw.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.nindel.vcf.gz",
        tbi = RES_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.nindel.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        zcat {input.vcf} \
             2>> {log.err} \
        | \
        awk '$0 !~/INDEL/ {{print}}' \
            2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule vcftools_get_inv_targetlib_targetseq:
    """
    Separate invariant and variant sites as subject to different filtering.
    """
    input:
        vcf = RES_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.nindel.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "vcftools_get_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv.vcf.gz",
        tbi = RES_DIR_RULE + "vcftools_get_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_get_inv_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_get_inv_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --max-maf 0 \
                 --recode \
                 --recode-INFO-all \
                 --stdout \
                 2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule vcftools_get_var_targetlib_targetseq:
    """
    Separate invariant and variant sites as subject to different filtering.
    """
    input:
        vcf = RES_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.nindel.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "vcftools_get_var_targetlib_targetseq/{interval_targetseq}.nindel_var.vcf.gz",
        tbi = RES_DIR_RULE + "vcftools_get_var_targetlib_targetseq/{interval_targetseq}.nindel_var.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_get_var_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_get_var_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --mac 1 \
                 --recode \
                 --recode-INFO-all \
                 --stdout \
                 2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule vcftools_filter_inv_targetlib_targetseq:
    """
    Apply filters to invariant sites.
    """
    input:
        vcf = RES_DIR_RULE + "vcftools_get_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "vcftools_filter_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv_filtered.vcf.gz",
        tbi = RES_DIR_RULE + "vcftools_filter_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv_filtered.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_filter_inv_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_filter_inv_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --max-missing {MIN_MISSING} \
                 --recode \
                 --recode-INFO-all \
                 --stdout \
                 2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule vcftools_filter_var_targetlib_targetseq:
    """
    Apply filters to variant sites.
    """
    input:
        vcf = RES_DIR_RULE + "vcftools_get_var_targetlib_targetseq/{interval_targetseq}.nindel_var.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "vcftools_filter_var_targetlib_targetseq/{interval_targetseq}.nindel_var_filtered.vcf.gz",
        tbi = RES_DIR_RULE + "vcftools_filter_var_targetlib_targetseq/{interval_targetseq}.nindel_var_filtered.vcf.gz.tbi"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_filter_var_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_filter_var_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --max-missing {MIN_MISSING} \
                 --minQ {MIN_QUAL} \
                 --recode \
                 --recode-INFO-all \
                 --stdout \
                 2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {output.vcf} \
              2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule bcftools_concat_filtered_targetlib_targetseq:
    """
    Concat the filtered invariant and variant sites VCF files.
    """
    input:
        vcf_inv = RES_DIR_RULE + "vcftools_filter_inv_targetlib_targetseq/{interval_targetseq}.nindel_inv_filtered.vcf.gz",
        vcf_var = RES_DIR_RULE + "vcftools_filter_var_targetlib_targetseq/{interval_targetseq}.nindel_var_filtered.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "bcftools_concat_filtered_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_concat_filtered_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_concat_filtered_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        bcftools concat --allow-overlaps \
                        {input.vcf_inv} \
                        {input.vcf_var} \
                        --output-type z \
                        --threads {threads} \
                        --output {output.vcf} \
                        2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """

rule bcftools_annotate_mismap_control_targetseq:
    """
    Annotate VCF file for filtering with bcftools.
    """
    input:
        mismap = RES_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.nindel.mismap",
        vcf = RES_DIR_RULE + "bcftools_concat_filtered_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "bcftools_annotate_mismap_control_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_annotate_mismap_control_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_annotate_mismap_control_targetseq/{interval_targetseq}.err"
    params:
        mismap_no_header = RES_DIR_RULE + "bcftools_annotate_mismap_control_targetseq/{interval_targetseq}.nindel.mismap_no_header.gz"
    shell:
        # Create temporary mismap file without header as required by bcftools annotate
        # Index temporary mismap file as required by bcftools annotate
        """
        tail -n +2 \
             {input.mismap} \
             2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {params.mismap_no_header}

        tabix --threads {threads} \
              --sequence 1 \
              --begin 2 \
              --end 2 \
              {params.mismap_no_header}
              1>> {log.out} \
              2>> {log.err}

        bcftools annotate --threads {threads} \
                          --annotations {params.mismap_no_header} \
                          --columns CHROM,POS,FORMAT/MMCT,INFO/MMMCT \
                          --header-lines {VCF_HEADER_MISMAP} \
                          --output-type z \
                          --output {output.vcf} \
                          {input.vcf} \
                          1>> {log.out} \
                          2>> {log.err}

        rm {params.mismap_no_header}
        rm "{params.mismap_no_header}.tbi"
        """

rule bcftools_annotate_gdepth_targetlib_targetseq:
    """
    Annotate VCF file for filtering with bcftools.
    """
    input:
        gdepth = RES_DIR_RULE + "r_apply_gdepth_filter_targetlib_targetseq/{interval_targetseq}.nindel.gdepth",
        vcf = RES_DIR_RULE + "bcftools_annotate_mismap_control_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "bcftools_annotate_gdepth_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_annotate_gdepth_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_annotate_gdepth_targetlib_targetseq/{interval_targetseq}.err"
    params:
        gdepth_no_header = RES_DIR_RULE + "bcftools_annotate_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth_no_header.gz"
    shell:
        # Create temporary gdepth file without header as required by bcftools annotate
        # Keep only CHROM, POS, and n_depth_pass columns as need to skip all the FORMAT/DP columns
        # Index temporary gdepth file as required by bcftools annotate
        """
        awk 'BEGIN {{OFS="\t"}}
             NR>1  {{print $1, $2, $NF}}' \
             {input.gdepth} \
             2>> {log.err} \
        | \
        bgzip --threads {threads} \
              --stdout \
              > {params.gdepth_no_header}

        tabix --threads {threads} \
              --sequence 1 \
              --begin 2 \
              --end 2 \
              {params.gdepth_no_header}
              1>> {log.out} \
              2>> {log.err}

        bcftools annotate --threads {threads} \
                          --annotations {params.gdepth_no_header} \
                          --columns CHROM,POS,INFO/NDP \
                          --header-lines {VCF_HEADER_GDEPTH} \
                          --output-type z \
                          --output {output.vcf} \
                          {input.vcf} \
                          1>> {log.out} \
                          2>> {log.err}

        rm {params.gdepth_no_header}
        rm "{params.gdepth_no_header}.tbi"
        """

rule bcftools_filter_targetlib_targetseq:
    """
    Filter VCF files based on MMMCT and NDP.
    """
    input:
        vcf = RES_DIR_RULE + "bcftools_annotate_gdepth_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    output:
        vcf = RES_DIR_RULE + "bcftools_filter_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "bcftools_filter_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "bcftools_filter_targetlib_targetseq/{interval_targetseq}.err"
    shell:
        """
        bcftools filter --threads {threads} \
                        --exclude 'INFO/MQ<{MIN_MQ} || INFO/MQ0F>{MAX_MQ0F} || INFO/MMMCT>{MAX_MMMCT} || INFO/NDP<{MIN_NDP}' \
                        --output-type z \
                        --output {output.vcf} \
                        {input.vcf} \
                        1>> {log.out} \
                        2>> {log.err}
        """

rule list_vcf_filtered_targetlib_targetseq:
    """
    List `.vcf.gz` files for merging.
    """
    input:
        vcf = expand(RES_DIR_RULE + "bcftools_filter_targetlib_targetseq/{interval_targetseq}.nindel_filtered.vcf.gz", interval_targetseq=INTERVALS_TARGETSEQ)
    output:
        txt = RES_DIR_RULE + "list_vcf_filtered_targetlib_targetseq/list_vcf.txt"
    run:
        f = open(output.txt, "w")
        for path in input.vcf:
            f.write(path + "\n")
        f.close()

rule gather_vcfs_filtered_targetlib_targetseq:
    """
    Merge `.vcf.gz` files.
    """
    input:
        txt = RES_DIR_RULE + "list_vcf_filtered_targetlib_targetseq/list_vcf.txt"
    output:
        vcf = RES_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/targetseq.nindel_filtered.vcf.gz",
        tbi = RES_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/targetseq.nindel_filtered.vcf.gz.tbi",
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/gather_vcfs_filtered_targetlib_targetseq.out",
        err = LOG_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/gather_vcfs_filtered_targetlib_targetseq.err"
    shell:
        """
        gatk GatherVcfs --INPUT {input.txt} \
                        --OUTPUT {output.vcf} \
                        1>> {log.out} \
                        2>> {log.err}

        tabix -p vcf {output.vcf} \
              1>> {log.out} \
              2>> {log.err}
        """
