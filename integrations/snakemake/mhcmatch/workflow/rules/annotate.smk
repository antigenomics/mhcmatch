"""Prior evidence and mimicry risk, both OFF by default and both reported beside the score.

Neither changes an ordering. A near-exact match to an already-tested neoantigen is prior evidence,
not a prediction, and is only meaningful for a cohort that did not contribute to the reference;
mimicry's two channel families carry opposite signs and are read separately or not at all.
"""


rule mhcmatch_neoag:
    input:
        lambda w: pool_for(w.arm, w.sample),
    output:
        f"{OUT}/{{arm}}/{{sample}}.mhc1.mhcmatch.neoag.tsv",
    params:
        subs=lambda w: opt("--max-subs", config["neoag_max_subs"]),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=240,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch neoag --peptides {input} --cls mhc1 {params.subs} --out {output}"


rule mhcmatch_mimicry:
    input:
        lambda w: pool_for(w.arm, w.sample),
    output:
        f"{OUT}/{{arm}}/{{sample}}.mhc1.mhcmatch.mimicry.tsv",
    params:
        ann=lambda w: flag("--annotate", config["mimicry_annotate"]),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=240,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch mimicry --peptides {input} --cls mhc1 {params.ann} --out {output}"
