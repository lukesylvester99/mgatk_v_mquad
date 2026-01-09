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


# --------- Rules ---------
rule all:
    input:
        os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results",
            f"{SAMPLE}_cells_for_mquad.tsv"
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
        cnv_calls            = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "cnv_calls.rds"),
        counts_gc_corrected  = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "counts_gc_corrected.rds"),
        count_summary        = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "count_summary.rds"),
        karyogram            = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "Karyogram.png"),
        results_gc_corrected = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "results_gc_corrected.rds"),
        results_table        = os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "results_table.tsv")
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
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results", "results_table.tsv"
        ),
        r_script="../scripts/get_epi_clones.R"
    params:
        out_dir=os.path.join(OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results")
    output:
        cancerous_cells=os.path.join(
            OUT_DIR, SAMPLE, "epi_results", "epiAneufinder_results",
            f"{SAMPLE}_cells_for_mquad.tsv"
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
