# Tal SINE — Talpidae Analysis

Analysis of **Tal SINEs** (tRNA-derived Short Interspersed Nuclear Elements) across Talpidae (moles, desmans, shrew-moles).

## Contents

| Path | Description |
|------|-------------|
| [WORKFLOW.md](WORKFLOW.md) | Full analysis workflow (Steps 1–8) |
| [db/](db/) | Tal SINE consensus databases used as search queries |
| [saq/](saq/) | *Scalopus aquaticus* results |
| [ccr/](ccr/) | *Condylura cristata* results |
| [toc/](toc/) | *Talpa occidentalis* subfamily alignments |
| [gpy/](gpy/) | *Galemys pyrenaicus* alignments |

## Species

| Code | Species | Common name | Assembly |
|------|---------|-------------|----------|
| saq  | *Scalopus aquaticus* | Eastern mole | [GCA_004024925.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_004024925.1/) |
| ccr  | *Condylura cristata* | Star-nosed mole | [GCF_000260355.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000260355.1/) |
| toc  | *Talpa occidentalis* | Western mole | TBD |
| gpy  | *Galemys pyrenaicus* | Pyrenean desman | [GCA_019455555.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_019455555.1/) |
| dmo  | *Desmana moschata* | Russian desman | TBD |

## Reports

- [saq — Scalopus aquaticus SINEderella report](saq/report.html) *(7 MB, open locally)*

## Analysis Status

| Species | Search | Subfamilies | Stats | Report | Tribes | PolyA | Term. | Orthologs |
|---------|:------:|:-----------:|:-----:|:------:|:------:|:-----:|:-----:|:---------:|
| saq | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| ccr | ✓ | ✓ | ✓ | ✓ | — | — | — | — |
| toc | ✓ | partial | — | — | — | — | — | — |
| gpy | partial | — | — | — | — | — | — | — |
| dmo | — | — | — | — | — | — | — | — |

## Tools

- [SINEderella](https://github.com/Toki-bio/SINEderella) — SINE search + subfamily assignment
- [SINEderella-dev](https://github.com/Toki-bio/SINEderella-dev) — development branch
