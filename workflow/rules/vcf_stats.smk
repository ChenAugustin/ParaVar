# Mismapping of controllib on targetseq
#--------------------------------------

rule vcftools_stat_gdepth_controllib_targetseq:
    """
    Extract read depth per individual per site across controllib and targetseq.
    This will be used for the mismapping filter.
    """
    input:
        vcf = RES_DIR_RULE + "awk_remove_indels_controllib_targetseq/{interval_targetseq}.nindel.vcf.gz",
    output:
        gdepth = RES_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.nindel.gdepth"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.err"
    params:
        basename = RES_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.nindel"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --geno-depth \
                 --out {params.basename} \
                 1>> {log.out} \
                 2>> {log.err} \
        """

rule vcftools_stat_gdepth_targetlib_targetseq:
    """
    Extract read depth per individual per site across targetlib and targetseq.
    This will be used for the depth and mismapping filters.
    """
    input:
        vcf = RES_DIR_RULE + "awk_remove_indels_targetlib_targetseq/{interval_targetseq}.nindel.vcf.gz",
    output:
        gdepth = RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.err"
    params:
        basename = RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel"
    shell:
        """
        vcftools --gzvcf {input.vcf} \
                 --geno-depth \
                 --out {params.basename} \
                 1>> {log.out} \
                 2>> {log.err} \
        """

rule samtools_depth_lib_controlseq:
    """
    Extract read depth per individual per site across controllib + targetlib and controlseq.
    This will be used to compute median read depth of each library across controlseq for the mismapping filter.
    """
    input:
        bam = expand(RES_DIR_RULE + config["mapper"] + "/{sample_id}/{sample_id}.bam", sample_id = sample_ids_controllib + sample_ids_targetlib),
        intervals = RES_DIR_RULE + "split_intervals_controlseq/{interval_controlseq}-scattered.interval_list"
    output:
        tsv = RES_DIR_RULE + "samtools_depth_lib_controlseq/{interval_controlseq}.raw.tsv"
    conda:
        "../envs/bcftools_call.yml"
    log:
        out = LOG_DIR_RULE + "samtools_depth_lib_controlseq/{interval_controlseq}.out",
        err = LOG_DIR_RULE + "samtools_depth_lib_controlseq/{interval_controlseq}.err"
    params:
        temp_interval = RES_DIR_RULE + "samtools_depth_lib_controlseq/{interval_controlseq}.txt"
    shell:
        # Save in temporary file:
        #     - interval file after removing the header (lines starting with "@") as not supported by samtools depth
        """
        awk 'BEGIN {{OFS="\\t"}};
             !/^@/ {{print $1, $2, $3}}' {input.intervals} \
            > {params.temp_interval}

        samtools depth -b {params.temp_interval} \
                       -o {output.tsv} \
                       {input.bam} \
                       1>> {log.out} \
                       2>> {log.err}

        rm {params.temp_interval}
        """

rule shuf_samtools_depth_lib_controlseq:
    """
    Randomly samples N sites across each interval of controllib + targetlib and controlseq.
    This will be used to compute median read depth of each library across controlseq for the mismapping filter.
    """
    input:
        tsv = expand(RES_DIR_RULE + "samtools_depth_lib_controlseq/{interval_controlseq}.raw.tsv", interval_controlseq=INTERVALS_CONTROLSEQ)
    output:
        tsv = RES_DIR_RULE + "shuf_samtools_depth_lib_controlseq/controlseq.shuf.tsv"
    log:
        err = LOG_DIR_RULE + "shuf_samtools_depth_lib_controlseq/shuf_samtools_depth_lib_controlseq.err"
    params:
        sample_ids_str = lambda wildcards: " ".join(sample_ids_controllib + sample_ids_targetlib)
    shell:
        # Add column names as not printed out by samtools depth
        """

        echo "CHROM" "POS" {params.sample_ids_str} | \
        awk 'BEGIN {{OFS = "\t"}} {{for (i = 1; i <= NF; i++) printf "%s%s", $i, (i == NF ? ORS : OFS)}}' - \
        1> {output.tsv} 2>> {log.err}

        for file in {input.tsv};
        do
            shuf -n {SHUF_N_LINES_PER_INTERVAL} \
                 $file \
                 2>> {log.err}
        done 1>> {output.tsv} \
             2>> {log.err}
        """

rule r_stat_median_depth_lib_controlseq:
    """
    Compute median read depth of each targetlib across a random set of sites from controlseq for the mismapping filter.
    """
    input:
        tsv = RES_DIR_RULE + "shuf_samtools_depth_lib_controlseq/controlseq.shuf.tsv"
    output:
        tsv = RES_DIR_RULE + "r_stat_median_depth_lib_controlseq/controlseq.shuf_median.tsv"
    conda:
        "../envs/r_tidyverse.yml"
    log:
        out = LOG_DIR_RULE + "r_stat_median_depth_lib_controlseq/r_stat_median_depth_lib_controlseq.out",
        err = LOG_DIR_RULE + "r_stat_median_depth_lib_controlseq/r_stat_median_depth_lib_controlseq.err"
    params:
        script = "workflow/scripts/r_stat_median_depth_lib_controlseq.R"
    shell:
        """
        Rscript {params.script} \
                {input.tsv} \
                {output.tsv} \
                1>> {log.out} \
                2>> {log.err}
        """

if config["run_samtools_depth_lib_controlseq"]:

    rule r_stat_mismap_control_targetseq:
        """
        Compute per-individual and mean mismapping rate of control reads on targetseq.
        """
        input:
            gdepth_controllib = RES_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.nindel.gdepth",
            gdepth_targetlib = RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth",
            tsv = RES_DIR_RULE + "r_stat_median_depth_lib_controlseq/controlseq.shuf_median.tsv"
        output:
            mismap = RES_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.nindel.mismap"
        conda:
            "../envs/r_tidyverse.yml"
        log:
            out = LOG_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.out",
            err = LOG_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.err"
        params:
            script_precise = "workflow/scripts/r_stat_mismap_control_targetseq_precise.R",
        shell:
            """
            Rscript {params.script_precise} \
                    {input.gdepth_controllib} \
                    {input.gdepth_targetlib} \
                    {SAMPLE_SHEET_PATH} \
                    {input.tsv} \
                    {output.mismap} \
                    1>> {log.out} \
                    2>> {log.err}
            """

else:

    rule r_stat_mismap_control_targetseq:
        """
        Compute per-individual and mean mismapping rate of control reads on targetseq.
        """
        input:
            gdepth_controllib = RES_DIR_RULE + "vcftools_stat_gdepth_controllib_targetseq/{interval_targetseq}.nindel.gdepth",
            gdepth_targetlib = RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth"
        output:
            mismap = RES_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.nindel.mismap"
        conda:
            "../envs/r_tidyverse.yml"
        log:
            out = LOG_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.out",
            err = LOG_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.err"
        params:
            script_approx = "workflow/scripts/r_stat_mismap_control_targetseq_approx.R"
        shell:
            """
            Rscript {params.script_approx} \
                    {input.gdepth_controllib} \
                    {input.gdepth_targetlib} \
                    {SAMPLE_SHEET_PATH} \
                    {FRACTION_TARGETSEQ_IN_TARGETLIB} \
                    {output.mismap} \
                    1>> {log.out} \
                    2>> {log.err}
            """

# Depth filter of targetlib on targetseq
#---------------------------------------

rule get_gdepth_targetlib_targetseq_singlecopy:
    """
    Get gdepth across target sequences presumed to be in single-copy (= not repeated within the target sequence).
    """
    input:
        gdepth = expand(RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth", interval_targetseq=INTERVALS_TARGETSEQ)
    output:
        gdepth = RES_DIR_RULE + "get_gdepth_targetlib_targetseq_singlecopy/singlecopy.nindel.gdepth"
    conda:
        "../envs/bedtools.yml"
    log:
        out = LOG_DIR_RULE + "get_gdepth_targetlib_targetseq_singlecopy/singlecopy.out",
        err = LOG_DIR_RULE + "get_gdepth_targetlib_targetseq_singlecopy/singlecopy.err"
    shell:
        # Scan all gdepth files, and print sites within high-confidence single-copy target regions
        # For this, convert the gdepth files to BED format on the fly for bedtools
        # Then, convert the output back to gdepth format
        # -wa Keep all columns of entry A in output
        # -wb Keep all columns of entry B in output
        """
        # Initiate increment just to print header once
        n=1

        for gdepth in {input.gdepth};
        do
            awk -v n=$n \
                'FNR==NR && FNR==1 && n==1 {{printf("%s\\t", $0)}};
                 FNR!=NR && FNR==1 && n==1 {{printf("%s\\n", $0)}}' \
                 $gdepth {TARGETSEQ_SINGLECOPY}

            awk -v n=$n \
                'NR>1 {{printf("%s\\t%s\\t%s", $1, $2-1, $2);
                        for(i=3; i<=NF; i++) {{printf("\\t%s", $i)}};
                        printf("\\n")}}' \
                $gdepth \
                2>> {log.err} \
            | \
            bedtools intersect -a stdin \
                               -b {TARGETSEQ_SINGLECOPY} \
                               -wa \
                               -wb \
                               2>> {log.err};

            ((n++))

        done \
        | \
        awk 'NR==1 {{print}};
             NR>1  {{printf("%s\\t%s", $1, $3);
                     for(i=4; i<=NF; i++) {{printf("\\t%s", $i)}};
                     printf("\\n")}}' - \
            1>> {output.gdepth} \
            2>> {log.err}
        """

rule get_mismap_targetlib_targetseq_singlecopy:
    """
    Get mismap across target sequences presumed to be in single-copy (= not repeated within the target sequence).
    """
    input:
        mismap = expand(RES_DIR_RULE + "r_stat_mismap_control_targetseq/{interval_targetseq}.nindel.mismap", interval_targetseq=INTERVALS_TARGETSEQ)
    output:
        mismap = RES_DIR_RULE + "get_mismap_targetlib_targetseq_singlecopy/singlecopy.nindel.mismap"
    conda:
        "../envs/bedtools.yml"
    log:
        out = LOG_DIR_RULE + "get_mismap_targetlib_targetseq_singlecopy/singlecopy.out",
        err = LOG_DIR_RULE + "get_mismap_targetlib_targetseq_singlecopy/singlecopy.err"
    shell:
        # Scan all mismap files, and print sites within high-confidence single-copy target regions
        # Do it the same way as for gdepth files
        """
        # Initiate increment just to print header once
        n=1

        for mismap in {input.mismap};
        do
            awk -v n=$n \
                'FNR==NR && FNR==1 && n==1 {{printf("%s\\t", $0)}};
                 FNR!=NR && FNR==1 && n==1 {{printf("%s\\n", $0)}}' \
                 $mismap {TARGETSEQ_SINGLECOPY}

            awk -v n=$n \
                'NR>1 {{printf("%s\\t%s\\t%s", $1, $2-1, $2);
                        for(i=3; i<=NF; i++) {{printf("\\t%s", $i)}};
                        printf("\\n")}}' \
                $mismap \
                2>> {log.err} \
            | \
            bedtools intersect -a stdin \
                               -b {TARGETSEQ_SINGLECOPY} \
                               -wa \
                               -wb \
                               2>> {log.err};

            ((n++))

        done \
        | \
        awk 'NR==1 {{print}};
             NR>1  {{printf("%s\\t%s", $1, $3);
                     for(i=4; i<=NF; i++) {{printf("\\t%s", $i)}};
                     printf("\\n")}}' - \
            1>> {output.mismap} \
            2>> {log.err}
        """

rule r_stat_median_gdepth_targetlib_targetseq_singlecopy:
    """
    Compute median read depth of each targetlib across target sequences presumed to be in single-copy (= not repeated within the target sequence).
    Do it with/without excluding all sites with mismapping of control reads on targetseq.
    """
    input:
        gdepth = RES_DIR_RULE + "get_gdepth_targetlib_targetseq_singlecopy/singlecopy.nindel.gdepth",
        mismap = RES_DIR_RULE + "get_mismap_targetlib_targetseq_singlecopy/singlecopy.nindel.mismap"
    output:
        median_gdepth_0_mismap = RES_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/by_individual.0_mismap.tsv",
        plot_gdepth_by_individual = RES_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/by_individual.pdf",
        plot_gdepth_by_region = RES_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/by_region.pdf",
        plot_gdepth_by_individual_by_region = RES_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/by_individual_by_region.pdf"
    conda:
        "../envs/r_tidyverse.yml"
    log:
        out = LOG_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/r_stat_median_gdepth_targetlib_targetseq_singlecopy.out",
        err = LOG_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/r_stat_median_gdepth_targetlib_targetseq_singlecopy.err"
    params:
        script = "workflow/scripts/r_stat_median_gdepth_targetlib_targetseq_singlecopy.R"
    shell:
        """
        Rscript {params.script} \
                {input.gdepth} \
                {input.mismap} \
                {FACTOR_MAX_GDEPTH} \
                {output.median_gdepth_0_mismap} \
                {output.plot_gdepth_by_individual} \
                {output.plot_gdepth_by_region} \
                {output.plot_gdepth_by_individual_by_region} \
                1>> {log.out} \
                2>> {log.err}
        """

rule r_apply_gdepth_filter_targetlib_targetseq:
    """
    For each site, compute the number of individuals whose gdepth falls within the range to be likely a single-copy region.
    Will use this statistic to filter on gdepth (= loose consensus approach) and keep only sites likely to be in a single-copy region.
    """
    input:
        gdepth = RES_DIR_RULE + "vcftools_stat_gdepth_targetlib_targetseq/{interval_targetseq}.nindel.gdepth",
        median_gdepth_0_mismap = RES_DIR_RULE + "r_stat_median_gdepth_targetlib_targetseq_singlecopy/by_individual.0_mismap.tsv"
    output:
        gdepth = RES_DIR_RULE + "r_apply_gdepth_filter_targetlib_targetseq/{interval_targetseq}.nindel.gdepth"
    conda:
        "../envs/r_tidyverse.yml"
    log:
        out = LOG_DIR_RULE + "r_apply_gdepth_filter_targetlib_targetseq/{interval_targetseq}.out",
        err = LOG_DIR_RULE + "r_apply_gdepth_filter_targetlib_targetseq/{interval_targetseq}.err"
    params:
        script = "workflow/scripts/r_apply_gdepth_filter_targetlib_targetseq.R"
    shell:
        """
        Rscript {params.script} \
                {input.gdepth} \
                {input.median_gdepth_0_mismap} \
                {MIN_GDEPTH} \
                {output.gdepth} \
                1>> {log.out} \
                2>> {log.err}
        """

# Nucleotide diversity
#---------------------

rule pixy_nucleotide_diversity_filtered_targetlib_targetseq:
    """
    Use pixy to compute nucleotide diversity across filtered sites of targetlib mapped against targetseq.
    """
    input:
        vcf = RES_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/targetseq.nindel_filtered.vcf.gz",
    output:
        pixy_bed = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.pixy.bed",
        txt = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.nindel_filtered_pi.txt"
    conda:
        "../envs/pixy.yml"
    log:
        out = LOG_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/pixy_nucleotide_diversity_filtered_targetlib_targetseq.out",
        err = LOG_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/pixy_nucleotide_diversity_filtered_targetlib_targetseq.err"
    params:
        out_dir = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq",
        out_prefix = "targetseq.nindel_filtered"
    shell:
        # Force pixy to analyze the whole target sequence, including regions not represented in the filtered VCF (because all sites were filtered)
        """
        awk -v w={PIXY_WINDOW_SIZE} \
            '{{
                chrom=$1;
                len=$3;
                for (start=1; start<=len; start+=w) {{
                                                       end=start+w-1;
                                                       if (end>len) end=len;
                                                       print chrom"\t"start"\t"end
                                                    }}
             }}' \
             {TARGETSEQ} \
             1>> {output.pixy_bed} \
             2>> {log.err}

        pixy --n_cores {threads} \
             --stats pi \
             --bed_file {output.pixy_bed} \
             --vcf {input.vcf} \
             --populations {PIXY_POPULATION_FILE} \
             --output_folder {params.out_dir} \
             --output_prefix {params.out_prefix} \
             1>> {log.out} \
             2>> {log.err}
        """
