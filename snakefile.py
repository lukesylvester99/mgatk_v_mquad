import os

# --------- Load Config Variables ---------
configfile: "config.yaml"

# global variables
SAMPLE = config["sample"]
BASE_DIR = os.path.expanduser(os.path.expandvars(config["base_dir"]))
OUT_DIR = config["out_dir"].format(base_dir=BASE_DIR)
OUT_DIR = os.path.expanduser(os.path.expandvars(OUT_DIR))

# epianeufinder pipeline paths
BLACKLIST = os.path.expanduser(os.path.expandvars(config["blacklist"]))
CLEANED_FRAGS = os.path.expanduser(os.path.expandvars(config["cleaned_frags"]))

# mgatk paths
BAM = os.path.expanduser(os.path.expandvars(config["bam"]))

#vireo trunk analysis paths
SCRIPTS_DIR = os.path.abspath(os.path.join(BASE_DIR, "workflows", "scripts"))
VIREO_TRUNK_PY = os.path.join(SCRIPTS_DIR, "vireo_trunk.py")

# --------- Rules ---------
rule all:
    input:
        # --- run_epi outputs ---
        epi_cnv_calls=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results","cnv_calls.rds"
        ),
        epi_counts_gc_corrected=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "counts_gc_corrected.rds"
        ),
        epi_count_summary=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "count_summary.rds"
        ),
        epi_karyogram=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "Karyogram.png"
        ),
        epi_results_gc_corrected=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "results_gc_corrected.rds"
        ),
        epi_results_table=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "results_table.tsv"
        ),

        # --- get_epi_clones output ---
        epi_cancerous_cells=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results",f"{SAMPLE}_cells_for_mquad.tsv"
        ),

        # --- mgatk output ---
        mgatk_done=os.path.join(OUT_DIR, SAMPLE, "mgatk", ".mgatk_done"),

        # --- process_mgatk outputs ---
        mgatk_cells=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "cells_passing_depth.tsv"),
        mgatk_varfilt=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "variant_stats_filtered.tsv"),
        mgatk_hetfilt=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "cell_heteroplasmy_filtered.tsv"),
        mgatk_summary=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "variant_summary.tsv"),

        # --- cellsnp outputs ---
        cellsnp_barcodes=os.path.join(OUT_DIR, SAMPLE, "mquad", "barcodes_epi_x_depth.tsv"),
        cellsnp_vcf_keep=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz"
        ),
        cellsnp_vcf_keep_tbi=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz.tbi"
        ),
        cellsnp_done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".cellsnp_done"),

        # --- mquad output ---
        mquad_done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".mquad_done"),

        # --- vireo outputs ---
        vireo_done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".vireo_done"),
        vireo_fit_log=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo", "fit_mito_clones.log"
        ),
        vireo_clone_assign=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo", "clone_assignments.tsv"
        ),

        # --- vireo_trunk outputs ---
        trunk_mean_af=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis", "vireo", "mean_AF.png"
        ),
        trunk_heatmap=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis", "vireo", "AF_heatmap.png"
        ),
        trunk_clone_prob=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis", "vireo", "clone_id_prob.tsv"
        ),
        trunk_clone_assign=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis", "vireo", "clone_assignments.tsv"
        ),
        trunk_beta_mu=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis", "vireo", "clone_beta_mu_estimated_AF.tsv"
        )


rule run_epi:
    """Run epiAneufinder and require all expected output files."""
    input:
        frags=CLEANED_FRAGS,
        blacklist=BLACKLIST,
        r_script="../scripts/run_epi.R"
    params:
        outdir=os.path.join(OUT_DIR, SAMPLE, "epi_results")
    output:
        cnv_calls            = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "cnv_calls.rds"),
        counts_gc_corrected  = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "counts_gc_corrected.rds"),
        count_summary        = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "count_summary.rds"),
        karyogram            = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "Karyogram.png"),
        results_gc_corrected = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "results_gc_corrected.rds"),
        results_table        = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "results_table.tsv")
    threads: 16
    conda:
        "../../envs/epi.yml"
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "run_epi.log")
    shell:
        r"""
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        Rscript {input.r_script} \
          --sample "{SAMPLE}" \
          --fragments "{input.frags}" \
          --outdir "{params.outdir}" \
          --blacklist "{input.blacklist}" \
          --threads {threads} \
          2>&1 | tee {log}
        """


rule get_epi_clones:
    """Run epiAneufinder clone calling and export cancerous cell barcodes."""
    input:
        results_table=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", "results_table.tsv"
        ),
        r_script="../scripts/get_epi_clones.R"
    params:
        out_dir=os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results")
    output:
        cancerous_cells=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results",
            "epiAneufinder_results", f"{SAMPLE}_cells_for_mquad.tsv"
        )
    threads: 16
    conda:
        "../../envs/epi.yml"
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "get_epi_clones.log")
    shell:
        r"""
        mkdir -p {params.out_dir}
        mkdir -p $(dirname {log})

        Rscript {input.r_script} \
          --sample "{SAMPLE}" \
          --results_table "{input.results_table}" \
          --out_dir "{params.out_dir}" \
          2>&1 | tee {log}
        """

rule mgatk:
    """
    Process mitochondrial data with mgatk for the epi-selected (cancerous) cells.
    """
    conda:
        "../../envs/mgatk.yml"
    input:
        bam=BAM,
        barcodes=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "epiAneufinder_results", f"{SAMPLE}_cells_for_mquad.tsv")
    output:
        done=os.path.join(OUT_DIR, SAMPLE, "mgatk", ".mgatk_done"),
        depth=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.depthTable.txt"),
        var=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.variant_stats.tsv.gz"),
        het=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.cell_heteroplasmic_df.tsv.gz")
    params:
        outdir=os.path.join(OUT_DIR, SAMPLE, "mgatk"),
        sample_name=SAMPLE
    threads: 16
    resources:
        mem_mb=64000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "mgatk.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log}")"

        (
          echo "==== mgatk START $(date) ===="
          echo "BAM: {input.bam}"
          echo "BARCODES: {input.barcodes}"
          echo "OUTDIR: {params.outdir}"
          echo

          cd "{params.outdir}"

          mgatk tenx \
            -i "{input.bam}" \
            -b "{input.barcodes}" \
            -bt CB \
            -n "{params.sample_name}" \
            -o "{params.outdir}" \
            --ncores {threads} \
            --keep-temp-files

          echo "==== mgatk END $(date) ===="
        ) &> "{log}"

        touch "{output.done}"
        """


rule process_mgatk:
    """
    Run mgatk.py on mgatk outputs to produce:
      - cells_passing_depth.tsv
      - variant_stats_filtered.tsv
      - cell_heteroplasmy_filtered.tsv
      - variant_summary.tsv
    """
    conda:
        "../../envs/mgatk.yml"
    input:
        done=os.path.join(OUT_DIR, SAMPLE, "mgatk", ".mgatk_done"),
        depth=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.depthTable.txt"),
        var=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.variant_stats.tsv.gz"),
        het=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final", f"{SAMPLE}.cell_heteroplasmic_df.tsv.gz"),
        script=os.path.join("..", "scripts", "mgatk.py")  
    output:
        cells=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "cells_passing_depth.tsv"),
        varfilt=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "variant_stats_filtered.tsv"),
        hetfilt=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "cell_heteroplasmy_filtered.tsv"),
        summary=os.path.join(OUT_DIR, SAMPLE, "process_mgatk", "variant_summary.tsv")
    params:
        mgatk_final=os.path.join(OUT_DIR, SAMPLE, "mgatk", "final"),
        outdir=os.path.join(OUT_DIR, SAMPLE, "process_mgatk")
    threads: 6
    resources:
        mem_mb=16000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "process_mgatk.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname "{log}")"

        (
          echo "==== process_mgatk START $(date) ===="
          echo "MGATK_FINAL: {params.mgatk_final}"
          echo "OUTDIR: {params.outdir}"
          echo

          # quick sanity checks so failures are obvious
          ls -lh "{params.mgatk_final}" || true
          [[ -s "{input.depth}" ]] || {{ echo "Missing/empty: {input.depth}" >&2; exit 2; }}
          [[ -s "{input.var}"   ]] || {{ echo "Missing/empty: {input.var}" >&2; exit 2; }}
          [[ -s "{input.het}"   ]] || {{ echo "Missing/empty: {input.het}" >&2; exit 2; }}
          echo

          python -u "{input.script}" \
            --mgatk_dir "{params.mgatk_final}" \
            --out_dir "{params.outdir}" \
            --sample "{SAMPLE}"

          echo
          echo "==== process_mgatk END $(date) ===="
        ) &> "{log}"
        """


rule cellsnp:
    """
    Run cellSNP-lite (2b -> site discovery), intersect barcodes (epi ∩ depth),
    filter sites to mgatk passing variants, then cellSNP-lite (1a -> genotyping).
    """
    conda:
        "../../envs/mquad.yml"
    input:
        bam=BAM,
        epi_barcodes=os.path.join(
            OUT_DIR, SAMPLE,
            "epi_results", "epiAneufinder_results", "epiAneufinder_results",
            f"{SAMPLE}_cells_for_mquad.tsv"
        ),
        depth_cells=os.path.join(
            OUT_DIR, SAMPLE,
            "process_mgatk", "cells_passing_depth.tsv"
        ),
        var_filt=os.path.join(
            OUT_DIR, SAMPLE,
            "process_mgatk", "variant_stats_filtered.tsv"
        ),
        script="../scripts/cellsnp.sh"
    output:
        # key outputs we expect from this stage
        barcodes=os.path.join(OUT_DIR, SAMPLE, "mquad", "barcodes_epi_x_depth.tsv"),
        vcf_keep=os.path.join(OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz"),
        vcf_keep_tbi=os.path.join(OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz.tbi"),
        done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".cellsnp_done")
    params:
        out_base=os.path.join(OUT_DIR, SAMPLE, "mquad"),
        mt_contig="chrM"
    threads: 20
    resources:
        mem_mb=64000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "cellsnp.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_base}" "$(dirname "{log}")"

        (
          echo "==== cellsnp_lite START $(date) ===="
          echo "BAM: {input.bam}"
          echo "EPI_BARCODES: {input.epi_barcodes}"
          echo "DEPTH_CELLS: {input.depth_cells}"
          echo "VARIANT_FILTER: {input.var_filt}"
          echo "OUT_BASE: {params.out_base}"
          echo "THREADS: {threads}"
          echo "MT_CONTIG: {params.mt_contig}"
          echo

          bash "{input.script}" \
            "{input.bam}" \
            "{input.epi_barcodes}" \
            "{input.depth_cells}" \
            "{input.var_filt}" \
            "{params.out_base}" \
            "{threads}" \
            "{params.mt_contig}"

          echo
          echo "==== cellsnp_lite END $(date) ===="
        ) &> "{log}"

        # sentinel
        touch "{output.done}"
        """

rule mquad:
    """
    Run MQuad on cellSNP-lite outputs.
    """
    conda:
        "../../envs/mquad.yml"
    input:
        done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".cellsnp_done"),
        script="../scripts/mquad.sh"
    output:
        done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".mquad_done"),
        passed=os.path.join(OUT_DIR, SAMPLE, "mquad", "mquad_mt", "passed_variant_names.txt")
    params:
        cellsnp_out=os.path.join(OUT_DIR, SAMPLE, "mquad", "cellsnp_per_cell_mtSNP"),
        mquad_out=os.path.join(OUT_DIR, SAMPLE, "mquad", "mquad_mt")
    threads: 20
    resources:
        mem_mb=64000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "mquad.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.mquad_out}" "$(dirname "{log}")"

        (
          echo "==== mquad START $(date) ===="
          echo "CELL_SNP_OUT: {params.cellsnp_out}"
          echo "MQUAD_OUT: {params.mquad_out}"
          echo "THREADS: {threads}"
          echo

          bash "{input.script}" \
            "{params.cellsnp_out}" \
            "{params.mquad_out}" \
            "{threads}"

          echo
          echo "==== mquad END $(date) ===="
        ) &> "{log}"

        touch "{output.done}"
        """

rule vireo:
    """
    Run vireo/clone fitting (fit_mito_clones.py) inside the MQuad output dir.
    """
    conda:
        "../../envs/mquad.yml"
    input:
        done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".mquad_done"),
        script="../scripts/vireo.sh"
    output:
        done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".vireo_done"),
        vireo_log=os.path.join(OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo", "fit_mito_clones.log"),
        clone_assign=os.path.join(OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo", "clone_assignments.tsv")
    params:
        mquad_out=os.path.join(OUT_DIR, SAMPLE, "mquad", "mquad_mt"),
        vireo_py=os.path.join(BASE_DIR, "fit_mito_clones.py")
    threads: 1
    resources:
        mem_mb=16000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "vireo.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{log}")"

        (
          echo "==== vireo START $(date) ===="
          echo "MQUAD_OUT: {params.mquad_out}"
          echo "VIERO_PY: {params.vireo_py}"
          echo

          bash "{input.script}" \
            "{params.mquad_out}" \
            "{params.vireo_py}"

          echo
          echo "==== vireo END $(date) ===="
        ) &> "{log}"

        # sanity check that the expected log exists
        [[ -s "{output.vireo_log}" ]] || {{ echo "ERROR: missing vireo log: {output.vireo_log}" >&2; exit 2; }}

        touch "{output.done}"
        """

rule vireo_trunk:
    """
    Reuse existing vireo clone assignments (from mquad run) to visualize/emit
    trunk-candidate mtSNP patterns by:
      1) taking the mgatk-filtered bulk VCF (cellSNP.base.filtered_to_mgatk.vcf.gz)
      2) removing mquad-informative variants (passed_variant_names.txt)
      3) re-running cellSNP-lite mode 1a to generate new AD/DP matrices
      4) running vireo_trunk.py on those new matrices with fixed clone labels
    """
    conda:
        "../../envs/mquad.yml"
    input:
        bam=BAM,

        # barcode list used for per-cell genotyping (already produced by cellsnp rule)
        barcodes=os.path.join(OUT_DIR, SAMPLE, "mquad", "barcodes_epi_x_depth.tsv"),

        # mgatk-filtered bulk VCF produced by cellsnp rule
        vcf_mgatk=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz"
        ),
        vcf_mgatk_tbi=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "cellsnp_bulk_mtSNP", "cellSNP.base.filtered_to_mgatk.vcf.gz.tbi"
        ),

        # list of mquad-informative variants to remove (lines like 14783T>C)
        mquad_pass=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "passed_variant_names.txt"
        ),

        # fixed clone assignments from original vireo (mquad-based)
        clone_assign=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo", "clone_assignments.tsv"
        ),

        # scripts
        trunk_sh=os.path.join(SCRIPTS_DIR, "vireo_trunk.sh"),
        trunk_py=VIREO_TRUNK_PY,

        # ensure upstream clone-fitting finished (keeps ordering consistent)
        vireo_done=os.path.join(OUT_DIR, SAMPLE, "mquad", ".vireo_done")
    output:
        mean_af=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt",
            "vireo_trunk_analysis", "vireo", "mean_AF.png"
        ),
        heatmap=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt",
            "vireo_trunk_analysis", "vireo", "AF_heatmap.png"
        ),
        clone_prob=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt",
            "vireo_trunk_analysis", "vireo", "clone_id_prob.tsv"
        ),
        clone_assign=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt",
            "vireo_trunk_analysis", "vireo", "clone_assignments.tsv"
        ),
        beta_mu=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt",
            "vireo_trunk_analysis", "vireo", "clone_beta_mu_estimated_AF.tsv"
        )
    params:
        workdir=os.path.join(
            OUT_DIR, SAMPLE, "mquad", "mquad_mt", "vireo_trunk_analysis"
        ),
        mt_contig="chrM"
    threads: 8
    resources:
        mem_mb=16000
    log:
        os.path.join(OUT_DIR, SAMPLE, "logs", "vireo_trunk.log")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.workdir}" "$(dirname "{log}")"

        (
          echo "==== vireo_trunk START $(date) ===="
          echo "WORKDIR: {params.workdir}"
          echo "BAM: {input.bam}"
          echo "BARCODES: {input.barcodes}"
          echo "VCF_MGATK: {input.vcf_mgatk}"
          echo "MQUAD_PASS: {input.mquad_pass}"
          echo "CLONE_ASSIGN: {input.clone_assign}"
          echo "TRUNK_SH: {input.trunk_sh}"
          echo "TRUNK_PY: {input.trunk_py}"
          echo "THREADS: {threads}"
          echo "MT_CONTIG: {params.mt_contig}"
          echo

          bash "{input.trunk_sh}" \
            "{input.bam}" \
            "{input.barcodes}" \
            "{input.vcf_mgatk}" \
            "{input.mquad_pass}" \
            "{input.clone_assign}" \
            "{params.workdir}" \
            "{input.trunk_py}" \
            "{threads}" \
            "{params.mt_contig}"

          echo
          echo "==== vireo_trunk END $(date) ===="
        ) &> "{log}"
        """
