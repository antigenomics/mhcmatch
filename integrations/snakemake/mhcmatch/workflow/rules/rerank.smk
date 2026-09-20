"""The rerank arm: your candidate table in, the same table plus an `mm_` block out."""


def _alleles_in(w):
    """The typing-derived allele file, unless one literal list was given for every sample."""
    if config.get("alleles" if w.cls == "mhc1" else "alleles_mhc2"):
        return []
    return [f"{OUT}/alleles/{w.sample}.{w.cls}.txt"] if "typing" in SAMPLES_MAP[w.sample] else []


rule mhcmatch_rerank:
    input:
        table=lambda w: SAMPLES_MAP[w.sample][f"table_{w.cls}"],
        # `--context` is not redundancy. A candidate table carries the MUTANT k-mer and nothing the
        # germline counterpart is recoverable from; the window FASTA carries the wild-type arm beside
        # it, which is where agretopicity comes from. Without it every row reads `wt_absent` --
        # correct, and a weaker model.
        context=lambda w: [SAMPLES_MAP[w.sample][f"fasta_{w.cls}"]]
        if f"fasta_{w.cls}" in SAMPLES_MAP[w.sample] else [],
        alleles=_alleles_in,
    output:
        f"{OUT}/rerank/{{sample}}.{{cls}}.epitopes.mhcmatch.tsv",
    params:
        prefix=lambda w: config["rerank_prefix"],
        ctx=lambda w, input: opt("--context", input.context[0] if input.context else None),
        tier=lambda w: opt("--tier", config["tier"]),
        tumor=lambda w: opt("--tumor", config.get("tumor")),
        prev=lambda w: opt("--prevalence", config.get("prevalence")),
        epitope=lambda w: opt("--epitope", config["rank_epitope"]),
        extra=lambda w: " ".join(filter(None, [
            flag("--extended", config["rank_extended"]),
            flag("--annotate", config["rank_annotate"]),
        ])),
    threads: 8
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch rank pairs {input.table} --cls {wildcards.cls} "
        "--passthrough --prefix {params.prefix} "
        "{params.ctx} {params.tier} {params.tumor} {params.prev} {params.epitope} {params.extra} "
        "--out {output}"
