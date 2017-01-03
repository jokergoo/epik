## epik

The **epik** package provides tools for integrative analysis and visualization of epigenomic sequencing data such as whole genome bisulfite sequencing, ChIP-Seq and RNA-Seq.

### Functions that deal with general genomic regions

- `basic_genomic_regions_stat()`
- `annotate_to_gene_models()`
- `annotate_to_genomic_features()`
- `genomic_regions_correlation()`
- `subgroup_specific_genomic_regions()`

### Functions that analyze methylation data

- `wgbs_qcplot()`
- `global_methylation_distribution()`
- `plot_multiple_samples_methylation_on_genome()`
- `heatmap_diff_methylation_in_genomic_features()`
- `methylation_subtype_classfication()`

### Functions that analyze ChIP-Seq data

- `chromatin_states_transition_chord_diagram()`

### Functions that do integrative analysis

- `correlated_regions()`
- `cr_enriched_at_tss()`
- `cr_hilbert()`
- `cr_gviz()`
- `enriched_heatmap_list_on_gene()`
- `enriched_heatmap_list_on_tss_cgi()`
- `enriched_heatmap_list_on_genomic_features()`
