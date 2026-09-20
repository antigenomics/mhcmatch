"""HLA typing file -> the allele list every scoring rule takes.

**The step whose absence is silent.** Every HLA caller writes the G-group form (`A*01:01:01G`) and
the pseudosequence tables are keyed at two fields, so an untrimmed name resolves to NOTHING --
and `Store._allele_set` drops what it cannot find without a word, so the run scores against an empty
panel and exits 0. `mhcmatch alleles` also splits the classes and joins the DP/DQ alpha-beta pairs,
neither of which a `cut -f2` does.
"""


rule mhcmatch_alleles:
    input:
        typing=lambda w: SAMPLES_MAP[w.sample]["typing"],
    output:
        f"{OUT}/alleles/{{sample}}.{{cls}}.txt",
    params:
        cls=lambda w: w.cls,
    threads: 1
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch alleles {input.typing} --cls {params.cls} --out {output}"
