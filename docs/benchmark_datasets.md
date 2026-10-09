# Benchmark datasets

## Location

`/projects/VitekLab/Data/MS/Benchmarking` on Explorer. 26 folders, roughly 37 GB.
Work against these paths rather than downloading.

## Scope

Label-free only. The four TMT folders and
`SingleCell-Leduc-MaxQ(TMT)` are out of scope. Single cell is an experimental
design rather than a format, so the remaining single cell folders are treated as
whatever tool produced them.

## Available formats

| Format | Datasets | Downstream ready |
|---|---|---|
| MaxQuant | 8 DDA folders | Only Puyvelde2022, which is the one folder with an MSstats annotation. Others have FragPipe's `experiment_annotation.tsv`, which would need reshaping |
| FragPipe | 8 DDA folders | Yes. `MSstats.csv` is FragPipe's own MSstats export and `FragPipetoMSstatsFormat` takes no annotation argument |
| DIA-NN | Gotti2021, Lou2022, Puyvelde2022, Dowell2021, Kalxdorf2021, Koopmans2023 | Gotti2021 and Puyvelde2022 only. The last three have no annotation file |
| Spectronaut | Navarro2016, Puyvelde2022, Dowell2021, Kalxdorf2021, Koopmans2023 | Navarro2016 and Puyvelde2022 only |
| MetaMorpheus | Solivais2024 | Yes |
| Proteome Discoverer | Payne (custom export) | No converter, see below |

Dowell2021, Kalxdorf2021, and Koopmans2023 each hold both DIA-NN and Spectronaut
output for the same raw data, which is useful for mapping accuracy even without an
annotation.

**Coverage gap:** MSstatsConvert supports roughly twelve label-free formats and there
is data here for five. The results table will either need additional data from
elsewhere or should report five formats and state the limitation.

## Reference pipelines

Eight R scripts on the server show how the lab ran these datasets. These are the
reference implementations to compare against.

| Converter | Arguments | Script location |
|---|---|---|
| `MaxQtoMSstatsFormat` | `(evidence, annot, pg)` | `DDA-Puyvelde2022/.../Maxquant/code-tony.R` |
| `FragPipetoMSstatsFormat` | `(input)`, no annotation | `DDA-Puyvelde2022/instrument_analysis.R` |
| `MetamorpheusToMSstatsFormat` | `(input, annot)` | `DDA-Solivais2024_Metamorpheus/Current/code-tony.R` |
| `DIANNtoMSstatsFormat` | `(input, annotation)` | `DIA-Gotti2021/example_diann_script.R` |
| `SpectronauttoMSstatsFormat` | `(input, annot)` | `DIA-Puyvelde2022-HYEtims735_DIA/.../code-tony.R` |

These are working files, not runnable code. One references an undefined object,
another contains local absolute paths. Read for intent.

`TOP0` and `TOP3` in FragPipe folders are quantification parameter variants of the
same dataset, not separate datasets.

## FDR baselines

The mixture scripts end with the same calculation. Human is constant background so
any significant human protein is a false positive, while yeast and E. coli are spiked
at known ratios and should come out significant:

```r
FDR = nrow(human) / (nrow(yeast) + nrow(human) + nrow(ecoli))
```

Recorded in script comments, computed with the hand written converters:

| Dataset | Path | FDR |
|---|---|---|
| Puyvelde2022 DDA | MaxQuant | 0.08 |
| Puyvelde2022 DIA | DIA-NN | 0.21 |
| Puyvelde2022 DIA | Spectronaut | 0.21 |

This is a stronger validation target than column matching: whether the science comes
out right, against a number already computed. `DDA-Solivais2024_Metamorpheus` goes
further with a contrast matrix encoding the designed spike ratios (1.5x, 2x, 2.5x,
3x), so recovered fold changes can be checked against designed ones.

## Three things that affect other tasks

**`SingleCell-Payne-Custom` is both the unseen format case and a reference for the
transformation executor.** It is the only dataset with no working converter, and
`msstatsconverter.R` does exactly what the executor needs to: `pivot_longer` over the
`Abundances:` columns for a wide to long reshape, then rename `Annotated.Sequence` to
`PeptideSequence`, `Master.Protein.Accessions` to `ProteinName`, and `Abundance` to
`Intensity`, with `ProductCharge`, `FragmentIon`, and `PrecursorCharge` set to NA and
`IsotopeLabelType` set to `L`.

**Pre-processing is an open scope question.** The Solivais script does substantial
work before calling the converter: dropping rows where `Protein Group` contains a
semicolon, dropping decoys, joining `QuantifiedProteins.tsv` to recover the organism,
then appending `|ECOLI` or `|HUMAN` to protein names. The species labels the FDR
calculation depends on do not exist in the raw file. Any comparison must apply
identical pre-processing to both paths or the numbers are not comparable. Proposed
position: this is curation rather than format mapping and stays outside the LLM's
remit. To be confirmed.

**The reference scripts give the pipeline, not the column mapping.** Ground truth for
mapping accuracy lives in the `clean_<Tool>` functions in MSstatsConvert's source,
which still need reading.