# Code availability

The following statement can be pasted into the manuscript. Replace
`sayaneshome` with your GitHub account and, after the first tagged
release, optionally mint a DOI via [Zenodo](https://zenodo.org/) and add it.

---

## Manuscript-ready statement (full workflow)

> **Code availability.** All code used for the bioinformatics analyses
> (differential expression with DESeq2, transcription-factor motif-enrichment
> analysis, gene-ontology/GSEA, and figure generation) and for the machine-learning
> analysis (random-forest and multitask models of transcription-factor–target gene
> expression) is openly available at
> https://github.com/sayaneshome/Sp3_manuscript. The repository
> includes the processed input data, a pinned computational environment
> (Python and R), and step-by-step instructions to reproduce all figures.

## Shorter variant (machine-learning focus)

> **Code availability.** The code for the machine-learning analysis, together with
> the full bioinformatics workflow and processed input data, is available at
> https://github.com/sayaneshome/Sp3_manuscript.

## With a Zenodo DOI (after release)

> **Code availability.** All analysis code is available at
> https://github.com/sayaneshome/Sp3_manuscript and archived at
> Zenodo (DOI: 10.5281/zenodo.XXXXXXX).

---

### How to add a Zenodo DOI
1. Log in to https://zenodo.org with your GitHub account.
2. Under **Settings → GitHub**, toggle the `Sp3_manuscript` repository **on**.
3. On GitHub, create a release (e.g. `v1.0.0`). Zenodo automatically archives it and
   issues a DOI; copy the DOI badge into `README.md` and the statement above.
