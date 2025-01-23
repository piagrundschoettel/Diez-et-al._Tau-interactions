# Diez-et-al._Tau-interactions
_Github repository for the bulk TE analysis, TF analysis and H3K9ac Cut&Tag analysis of Diez et al. 2025 by Pia Grundschoettel, Koki Sakurai, and Andrew Newman._

In this repository, you can find all the code that has been written to analyse and visualize the bulk RNA-seq data and H3K9ac Cut&Tag data in Diez et al. 2025:
__Tau interactions with the nuclear envelope modulate neuronal chromatin and nuclear stability__
For this project the data is available under the following GEO accession number: GSE283514. Hippocampal mouse neurons were genetically modified with adeno-associated viruses to express and overexpress human Tau and knockdown mouse Tau to investigate the consequences of Tau overexpression or knockdown. 

## bulkRNA-seq_analysis

The alignment for the bulkRNA-seq data was done with the settings --outFilterMultimapNmax 100, --winAchnorMultimapNmax 100 to maximize repeat mapping. The TETranscripts package was used to generate the count table with the docker container mhammelllab/tetranscripts, which can be found [here](https://hubgw.docker.com/r/mhammelllab/tetranscripts). For the quanitfication of both gene and TE transcripts the M25 annotation was used and the TE annotation from hammel labs [mm10_rmsk_TE.gtf](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/).
```
singularity exec --bind [necessary folders]  tetranscripts.sif TEcount --sortByPos --format BAM --mode multi -b [sample.bam] --GTF M25_annotation.gtf --TE mm10_rmsk_TE.gtf --project sample_name --outdir /output
```
The following bulk RNA-seq analysis have been performed in the docker container jsschrepping/r_docker:jss_R430_bioc317_slim, which can be found [here](https://hub.docker.com/layers/jsschrepping/r_docker/jss_R430_bioc317_slim/images/sha256-43ffc9fe6d3951590c2b0158d764de9e48ace35a4d3756eff0acf3ce6454baa2)

## H3K9ac Cut&Tag

