# _Klebsiella pneumoniae_ MAGs

Repository with custom code used for investigating the genomic and functional diversity of gut-associated _Klebsiella pneumoniae_ from MAGs and isolate genomes.

## Workflows

[ML-Microbiome - Machine learning classification of microbiome data](https://github.com/alexmsalmeida/ml-microbiome)

[GenoFan - Genome functional annotation pipeline](https://github.com/alexmsalmeida/genofan)


## Post-processing scripts

These scripts were used for downstream analysis and plotting:

#### <ins>Phylogenetic analysis</ins>

* [itol_gwas.R](scripts/itol_gwas.R)
* [itol_species.R](scripts/itol_species.R)
* [itol_novel.R](scripts/itol_novel.R)
* [metadata_adonis.R](scripts/metadata_adonis.R)
* [metadata_sankey.R](scripts/metadata_sankey.R)
* [tanglegram_panaroo-vs-snippy.R](scripts/tanglegram_panaroo-vs-snippy.R)
* [tanglegram_pan-vs-core.R](scripts/tanglegram_pan-vs-core.R)
* [pd_fold-change.R](scripts/pd_fold-change.R)
* [pd_fold-change_country.R](scripts/pd_fold-change_country.R)

#### <ins>Pan-genome analysis</ins>

* [panaroo_cog-plot.R](scripts/panaroo_cog-plot.R)
* [panaroo_combined.R](scripts/panaroo_combined.R)
* [panaroo_params.R](scripts/panaroo_params.R)
* [mags_genes.R](scripts/mags_genes.R)

#### <ins>Machine Learning modelling</ins>

* [ml_boxplots.R](scripts/ml_boxplots.R)
* [ml_mags-vs-isolates.R](scripts/ml_mags-vs-isolates.R)
* [ml_roc-curves.R](scripts/ml_roc-curves.R)

#### <ins>Microbial Genome-Wide Association Study (mGWAS)</ins>

* [pyseer-to-cog.R](scripts/pyseer-to-cog.R)
* [pyseer_volcano-cog.R](scripts/pyseer_volcano-cog.R)

#### <ins>Other</ins>

* [STs_gtype.R](scripts/STs.R)
* [mash_refseq.R](scripts/mash_refseq.R)
