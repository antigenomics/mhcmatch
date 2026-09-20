// The SELECTION seam: the one place a host's construct changes because of this model.
//
// **A reranked table published beside the host's own changes nothing.** A host pipeline picks the
// peptides that go into its construct in its own selector, and a selector sorts on a column whose
// name it fixes -- typically a literal `score`. `rank pairs --passthrough --prefix mm_` adds
// `mm_score` and leaves that column exactly as it found it, which is the right default and is also
// why `annotate` is a genuine no-op on the deliverable. Something has to write it.
//
// This process is that something, and it is the whole of `rerank` mode.
//
// **It preserves every multiplicative correction the host applied, because those are clinical and
// this model does not contain them.** The host's own scorer builds
//
//     score = base(affinity, agretopicity, expression) x driver x recurrence x hla_freq x hla_loss
//
// and only `base` is a question about peptides. A driver boost of 100.0 and an HLA-loss weight of
// 0.0 for an allele the tumour has lost are both ordinary settings -- the first is a deliberate two
// orders of magnitude, the second a veto. Replacing the product would put an epitope restricted to a
// deleted allele back into a construct. So the base is replaced and the factors are carried:
//
//     score := W x minmax(mm_score) x driver(driver_class) x hla_weight x recurrence_weight
//              x hla_frequency_weight
//
// W is the host's own total base weight, so the numbers stay in the range its operators read. Each
// factor is used only if its column is present -- a host that ran no recurrence library has no
// `recurrence_weight` and gets 1.0, which is what its own scorer would have done.
//
// The host's value is kept as `score_host`. Nothing is destroyed, and the two are one `paste` apart
// for anyone who wants to see what moved.
//
// **Every host-specific name here is a parameter with a default**, so the process carries no
// knowledge of any particular pipeline: `mhcmatch_overlay_score_column` is the column written,
// `mhcmatch_overlay_rescore_column` the one read, `mhcmatch_overlay_keep_weights` the factors
// carried, and the four base weights and two driver boosts are read from the host's own params when
// they happen to be loaded in the same session.

process MHCMATCH_RESCORE {
    tag "${meta.id}:${cls}"
    label 'process_single'

    conda "${moduleDir}/../mhcmatch/environment.yml"

    input:
    tuple val(meta), path(reranked), val(cls)

    output:
    tuple val(meta), path("*.epitopes.mhcmatch.scored.csv"), val(cls), emit: csv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pre    = params.mhcmatch_rerank_prefix ?: 'mm_'
    def src    = params.mhcmatch_overlay_rescore_column ?: "${pre}score"
    // The host's four base weights. Read from ITS params when they are loaded in the same session
    // -- `params.foo` on a key nobody set is null, so the fallbacks are what a standalone run gets.
    def total  = ((params.weight_affinity ?: 4.0) as double) +
                 ((params.weight_agretopicity ?: 1.0) as double) +
                 ((params.weight_expr_gene ?: 2.0) as double) +
                 ((params.weight_expr_local ?: 0.0) as double)
    def bd     = params.driver_boost_driver ?: 1.0
    def bl     = params.driver_boost_likely_driver ?: 1.0
    def keep   = params.mhcmatch_overlay_keep_weights ?:
                     'hla_weight,recurrence_weight,hla_frequency_weight'
    def dst    = params.mhcmatch_overlay_score_column ?: 'score'
    """
    python - ${reranked} ${prefix}.${cls}.epitopes.mhcmatch.scored.csv \\
            '${src}' '${total}' '${bd}' '${bl}' '${keep}' '${dst}' <<'PY'
import csv, sys

src_path, dst_path, score_col, total, bd, bl, keep, out_col = sys.argv[1:9]
total = float(total)
boost = {"driver": float(bd), "likely_driver": float(bl)}
factors = [c for c in keep.split(",") if c.strip()]

with open(src_path, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    cols = list(reader.fieldnames or [])
    rows = list(reader)

if score_col not in cols:
    sys.exit("MHCMATCH_RESCORE: %s has no column %r -- rerank mode needs the reranked table, "
             "not the host's own." % (src_path, score_col))

kept_col = out_col + "_host"
out_cols = cols + [c for c in (out_col, kept_col) if c not in cols]

def num(cell):
    try:
        return float(cell)
    except (TypeError, ValueError):
        return None

vals = [v for v in (num(r.get(score_col)) for r in rows) if v is not None]
lo, hi = (min(vals), max(vals)) if vals else (0.0, 0.0)
span = hi - lo

unscored = 0
for r in rows:
    v = num(r.get(score_col))
    if v is None:
        unscored += 1
        base = 0.0
    elif span == 0.0:
        # One distinct value: the host's own scaler gives every row 1.0 here, not 0.0.
        base = total
    else:
        base = total * (v - lo) / span
    factor = boost.get((r.get("driver_class") or "").strip(), 1.0)
    for c in factors:
        w = num(r.get(c))
        if w is not None:
            factor *= w
    r[kept_col] = r.get(out_col, "")
    r[out_col] = "%.6f" % (base * factor)

rows.sort(key=lambda r: float(r[out_col]), reverse=True)

with open(dst_path, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=out_cols, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print("MHCMATCH_RESCORE: %d row(s), %s := %.1f x minmax(%s) x %s; %d row(s) unscored -> 0 "
      "(previous value kept as %s)"
      % (len(rows), out_col, total, score_col,
         " x ".join(["driver(driver_class)"] + factors), unscored, kept_col),
      file=sys.stderr)
PY
    """

    stub:
    """
    head -1 ${reranked} | tr '\\t' ',' > ${task.ext.prefix ?: meta.id}.${cls}.epitopes.mhcmatch.scored.csv
    """
}
