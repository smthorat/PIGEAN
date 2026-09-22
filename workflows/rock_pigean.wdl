version 1.0

workflow rock_pigean {
    input {
        File input_file
        String analysis_type = "positive-controls"
        String preset = "standard"
        String gene_set_profile = "default"
        String genome_build = "hg19"
        Array[File] custom_gene_set_files = []
        String custom_gene_set_format = "engine"
        String custom_gene_set_action = "supplement"
        File? background_gene_list

        # Phase 3: Column mapping for tabular evidence files
        String? gene_column
        String? score_column

        # Phase 3: Percentile-specific
        Boolean higher_is_better = false

        # Phase 3: Exome-specific column overrides
        String? exomes_gene_col
        String? exomes_p_col
        String? exomes_beta_col
        String? exomes_se_col
        String? exomes_n_col
        Float? exomes_n

        # Phase 4: GWAS column overrides
        String? gwas_chrom_col
        String? gwas_pos_col
        String? gwas_p_col
        String? gwas_beta_col
        String? gwas_se_col
        String? gwas_n_col
        Float? gwas_n
        String? gwas_freq_col
        String? gwas_locus_col
        String? gwas_filter_col
        String? gwas_filter_value

        # Phase 5: Convergence diagnostics
        Boolean enable_convergence_trace = false

        Int memory_gb = 8
        Int cpu = 4
        Int disk_gb = 50
    }
    call run_pigean {
        input:
            input_file = input_file,
            analysis_type = analysis_type,
            preset = preset,
            gene_set_profile = gene_set_profile,
            genome_build = genome_build,
            custom_gene_set_files = custom_gene_set_files,
            custom_gene_set_format = custom_gene_set_format,
            custom_gene_set_action = custom_gene_set_action,
            background_gene_list = background_gene_list,
            gene_column = gene_column,
            score_column = score_column,
            higher_is_better = higher_is_better,
            exomes_gene_col = exomes_gene_col,
            exomes_p_col = exomes_p_col,
            exomes_beta_col = exomes_beta_col,
            exomes_se_col = exomes_se_col,
            exomes_n_col = exomes_n_col,
            exomes_n = exomes_n,
            gwas_chrom_col = gwas_chrom_col,
            gwas_pos_col = gwas_pos_col,
            gwas_p_col = gwas_p_col,
            gwas_beta_col = gwas_beta_col,
            gwas_se_col = gwas_se_col,
            gwas_n_col = gwas_n_col,
            gwas_n = gwas_n,
            gwas_freq_col = gwas_freq_col,
            gwas_locus_col = gwas_locus_col,
            gwas_filter_col = gwas_filter_col,
            gwas_filter_value = gwas_filter_value,
            enable_convergence_trace = enable_convergence_trace,
            memory_gb = memory_gb,
            cpu = cpu,
            disk_gb = disk_gb,
    }
    output {
        File gs = run_pigean.gs
        File gss = run_pigean.gss
        File ggss = run_pigean.ggss
        File params = run_pigean.params
        File manifest = run_pigean.manifest
        File? gene_qc = run_pigean.gene_qc
        File? evidence_qc = run_pigean.evidence_qc
        File report = run_pigean.report
        File? report_html = run_pigean.report_html
        File? convergence_json = run_pigean.convergence_json
        File resolved_config = run_pigean.resolved_config
    }
}

task run_pigean {
    input {
        File input_file
        String analysis_type
        String preset
        String gene_set_profile
        String genome_build
        Array[File] custom_gene_set_files
        String custom_gene_set_format
        String custom_gene_set_action
        File? background_gene_list

        # Phase 3: Column mapping for tabular evidence files
        String? gene_column
        String? score_column

        # Phase 3: Percentile-specific
        Boolean higher_is_better

        # Phase 3: Exome-specific column overrides
        String? exomes_gene_col
        String? exomes_p_col
        String? exomes_beta_col
        String? exomes_se_col
        String? exomes_n_col
        Float? exomes_n

        # Phase 4: GWAS column overrides
        String? gwas_chrom_col
        String? gwas_pos_col
        String? gwas_p_col
        String? gwas_beta_col
        String? gwas_se_col
        String? gwas_n_col
        Float? gwas_n
        String? gwas_freq_col
        String? gwas_locus_col
        String? gwas_filter_col
        String? gwas_filter_value

        # Phase 5: Convergence diagnostics
        Boolean enable_convergence_trace

        Int memory_gb
        Int cpu
        Int disk_gb
    }
    runtime {
        docker: "gcr.io/nitrogenase-docker/rock-pigean:5.0.0"
        memory: memory_gb + " GB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
    }
    command <<<
        python3 -u /app/run_pigean.py \
            --analysis ~{analysis_type} \
            --input ~{input_file} \
            --gene-sets ~{gene_set_profile} \
            --genome-build ~{genome_build} \
            --preset ~{preset} \
            --output results \
            ~{if length(custom_gene_set_files) > 0 then "--custom-gene-set-files " + sep(" ", custom_gene_set_files) + " --custom-gene-set-format " + custom_gene_set_format + " --custom-gene-set-action " + custom_gene_set_action else ""} \
            ~{if defined(background_gene_list) then "--background " + background_gene_list else ""} \
            ~{if defined(gene_column) then "--gene-column " + gene_column else ""} \
            ~{if defined(score_column) then "--score-column " + score_column else ""} \
            ~{if higher_is_better then "--higher-is-better" else ""} \
            ~{if defined(exomes_gene_col) then "--exomes-gene-col " + exomes_gene_col else ""} \
            ~{if defined(exomes_p_col) then "--exomes-p-col " + exomes_p_col else ""} \
            ~{if defined(exomes_beta_col) then "--exomes-beta-col " + exomes_beta_col else ""} \
            ~{if defined(exomes_se_col) then "--exomes-se-col " + exomes_se_col else ""} \
            ~{if defined(exomes_n_col) then "--exomes-n-col " + exomes_n_col else ""} \
            ~{if defined(exomes_n) then "--exomes-n " + exomes_n else ""} \
            ~{if defined(gwas_chrom_col) then "--gwas-chrom-col " + gwas_chrom_col else ""} \
            ~{if defined(gwas_pos_col) then "--gwas-pos-col " + gwas_pos_col else ""} \
            ~{if defined(gwas_p_col) then "--gwas-p-col " + gwas_p_col else ""} \
            ~{if defined(gwas_beta_col) then "--gwas-beta-col " + gwas_beta_col else ""} \
            ~{if defined(gwas_se_col) then "--gwas-se-col " + gwas_se_col else ""} \
            ~{if defined(gwas_n_col) then "--gwas-n-col " + gwas_n_col else ""} \
            ~{if defined(gwas_n) then "--gwas-n " + gwas_n else ""} \
            ~{if defined(gwas_freq_col) then "--gwas-freq-col " + gwas_freq_col else ""} \
            ~{if defined(gwas_locus_col) then "--gwas-locus-col " + gwas_locus_col else ""} \
            ~{if defined(gwas_filter_col) then "--gwas-filter-col " + gwas_filter_col else ""} \
            ~{if defined(gwas_filter_value) then "--gwas-filter-value " + gwas_filter_value else ""} \
            ~{if enable_convergence_trace then "--enable-convergence-trace" else ""}
    >>>
    output {
        File gs = "results/gs.out"
        File gss = "results/gss.out"
        File ggss = "results/ggss.out"
        File params = "results/p.out"
        File manifest = "results/run_manifest.json"
        File? gene_qc = "results/input_gene_qc.tsv"
        File? evidence_qc = "results/input_evidence_qc.tsv"
        File report = "results/report.txt"
        File? report_html = "results/report.html"
        File? convergence_json = "results/convergence.json"
        File resolved_config = "results/resolved_config.json"
    }
}
