# Proteomics validation preparation

`prepare_proteomics_queries.py` converts the supplied two-column protein
spectral-count files into canonical pCTSEA gene queries. It is a validation
tool and is not required to run pCTSEA.

The tool writes two query variants:

- `raw` contains gene-level spectrum counts.
- `nsaf` contains length-corrected values. For each protein, SAF is spectrum
  count divided by protein sequence length. Each gene value is then divided by
  the sample SAF total and scaled to one million.

Multiple accessions mapped to the same gene use the maximum value by default.
This conservative rule avoids summing repeated evidence from isoforms. Use
`--aggregation sum` only when the input protein inference supports summation.

Example:

```bash
python3 tools/prepare_proteomics_queries.py \
  --input-dir /path/to/proteomics/files \
  --fasta /path/to/pinned-human-uniprot.fasta \
  --output-dir validation/proteomics_queries
```

Use `--top-genes 100` to create a deliberately capped smoke-test corpus. Every
run writes `summary.tsv`, `unresolved.tsv`, and `provenance.tsv` so filtering and
mapping decisions remain inspectable.

This process does not create intensities. NSAF is a length correction for
spectral counts and is only an approximation when peptide-level assignments are
unavailable.
