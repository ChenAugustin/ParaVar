rule circos_plot_stats_filtered_targetlib_targetseq:
    """
    Visualize various statistics along the target sequence using circos.
    """
    input:
        txt = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.nindel_filtered_pi.txt"
    output:
        svg = RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/circos.svg"
    conda:
        "../envs/circos.yml"
    log:
        out = LOG_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/circos_plot_stats_filtered_targetlib_targetseq.out",
        err = LOG_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/circos_plot_stats_filtered_targetlib_targetseq.err"
    params:
        tsv_scorable = RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/scorability.txt",
        tsv_pi = RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/pi.txt",
        config = RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/circos_config_targets",
        out_dir = RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq"
    shell:
        # Reformat input data for Circos
        """
        awk 'BEGIN {{OFS="\\t"}}
             NR>1  {{print $2, $3, $4, $6 > "{params.tsv_scorable}";
                     print $2, $3, $4, $5 > "{params.tsv_pi}"}}' \
            {input.txt}

        sed 's,tsv_scorable,'"{params.tsv_scorable}"',' \
            {CIRCOS_CONF_TARGETLIB_TARGETSEQ} \
            2>> {log.err} \
        | \
        sed 's,tsv_pi,'"{params.tsv_pi}"',' \
            - \
            2>> {log.err} \
        | \
        sed 's,out_dir,'"{params.out_dir}"',' \
            - \
            1> {params.config} \
            2>> {log.err};

        circos -conf {params.config} \
               -param karyotype={CIRCOS_KARYOTYPE_TARGETSEQ} \
               -outputdir {params.out_dir} \
               1>> {log.out} \
               2>> {log.err}
        """

rule rideogram_plot_stats_filtered_targetlib_targetseq:
    """
    Visualize various statistics along the target sequence using RIdeogram.
    """
    input:
        pixy_bed = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.pixy.bed",
        txt = RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.nindel_filtered_pi.txt"
    output:
        svg_scor = RES_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/scorability.svg",
        svg_pi = RES_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/pi.svg"
    conda:
        "../envs/rideogram.yml"
    log:
        out = LOG_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/rideogram_plot_stats_filtered_targetlib_targetseq.out",
        err = LOG_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/rideogram_plot_stats_filtered_targetlib_targetseq.err"
    params:
        script = "workflow/scripts/rideogram_plot_stats_filtered_targetlib_targetseq.R"
    shell:
        """
        Rscript {params.script} \
                {TARGETSEQ} \
                {input.pixy_bed} \
                {input.txt} \
                {output.svg_scor} \
                {output.svg_pi} \
                1>> {log.out} \
                2>> {log.err}
        """
