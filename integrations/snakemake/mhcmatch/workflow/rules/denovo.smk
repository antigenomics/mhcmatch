"""The de novo arm: your mutation-window FASTA in, our epitope table out."""


def _alleles_str(w, input):
    """`--alleles` as a literal, from the config or from the file `mhcmatch alleles` wrote."""
    lit = config.get("alleles" if w.cls == "mhc1" else "alleles_mhc2")
    if lit:
        return f"--alleles '{lit}'"
    return f"--alleles \"$(cat {input.alleles[0]})\"" if input.alleles else ""


def _alleles_in(w):
    if config.get("alleles" if w.cls == "mhc1" else "alleles_mhc2"):
        return []
    return [f"{OUT}/alleles/{w.sample}.{w.cls}.txt"] if "typing" in SAMPLES_MAP[w.sample] else []


rule mhcmatch_predict:
    """Per-allele binding over the same windows, in the generic native TSV.

    `--native` is the output to read. `--scored-csv` beside it is a **legacy wide-CSV
    compatibility export** in a fixed column schema, kept for a caller whose downstream already
    reads that shape; nothing here consumes it and it is not a default deliverable.
    """
    input:
        fasta=lambda w: SAMPLES_MAP[w.sample][f"fasta_{w.cls}"],
        alleles=_alleles_in,
    output:
        scored=f"{OUT}/denovo/{{sample}}.{{cls}}.mhcmatch.scored.csv",
        native=f"{OUT}/denovo/{{sample}}.{{cls}}.mhcmatch.native.tsv",
    params:
        alleles=_alleles_str,
        tier=lambda w: opt("--tier", config["tier"]),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch predict {input.fasta} {params.alleles} --cls {wildcards.cls} {params.tier} "
        "--scored-csv {output.scored} --native {output.native}"


rule mhcmatch_rank:
    """The fitted EPIC aggregate over the same windows, one ordered table."""
    input:
        fasta=lambda w: SAMPLES_MAP[w.sample][f"fasta_{w.cls}"],
        alleles=_alleles_in,
    output:
        f"{OUT}/denovo/{{sample}}.{{cls}}.mhcmatch.ranked.tsv",
    params:
        alleles=_alleles_str,
        tier=lambda w: opt("--tier", config["tier"]),
        tumor=lambda w: opt("--tumor", config.get("tumor")),
        prev=lambda w: opt("--prevalence", config.get("prevalence")),
        thr=lambda w: opt("--threshold", None if config["rank_threshold"] == "none"
                          else config["rank_threshold"]),
        epitope=lambda w: opt("--epitope", config["rank_epitope"]),
    threads: 8
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch rank fasta {input.fasta} {params.alleles} --cls {wildcards.cls} "
        "{params.tier} {params.tumor} {params.prev} {params.thr} {params.epitope} --out {output}"
