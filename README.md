# APPLd_paper
Repository for files and scripts associated with APPLd manuscript. All possible files archived here (Github), files too large reside on Zenodo (https://doi.org/10.5281/zenodo.17751181).

Contents: 

folder:,File name:,Description:,File #:,GitHub?:,Zenodo?:,Associated Figure:
folder 1,SampleQC_precleanup.pdf,Sample QC summary of all findings on the initial submitted sample before AMPure cleanup,supp file 1,Yes,Yes,Figure 1
,SampleQC_data_precleanup.pdf,Sample QC data supporting data summary,supp file 2,Yes,Yes,
,AZENTA_NGS_PacBio_Data_Report.pdf,Sequencing run report data ,supp file 3,Yes,Yes,
,APPLdM.hifireads.fastq.gz,HiFi raw reads from APPLd PacBio Revio,supp file 4,No,Yes,
folder 2,APPLd_assembly.merfin_racon2.fa,Final Merfin-polished genome assembly,supp file 5,No,Yes,
,appld_reads.histo.txt,Jellyfish histogram input for GenomeScope,supp file 6,Yes,Yes,
,quast_output/ directory,Full QUAST stats and misassembly info,supp file 7,No,Yes,
,appld_busco/ directory,Raw BUSCO analysis output files,supp file 8,No,Yes,
,coverage_bins.txt,Raw per-bin depth data for heatmap,supp file 9,Yes,Yes,
,coverage.tsv,Samtools per-region coverage summary,supp file 10,No,Yes,
,primary_assembly_QUAST.pdf,Quality assesment pre polishing,supp file 11,Yes,Yes,
,polished_assembly_QUAST.pdf,Quality assesment post polishing,supp file 12,Yes,Yes,
folder 3,appl_vs_appld.blast,BLAST alignment of Appl CDS against the APPLd genome (for panel D),supp file 13,Yes,Yes,Figure 2
,dm6_r6.62.fa,"Reference genome used for mapping, variant calling, and visualizations",supp file 14,No,Yes,
,dmel-all-r6.62.gtf.gz,Reference genome annotation file used for gene track visualization,supp file 15,Yes,Yes,
,appld_freebayes.pacbio.vcf.gz,Variant calls (FreeBayes) mapped to the reference for VCF visualization,supp file 16,No,Yes,
,appld_longshot.vcf.gz,"Variant calls (Longshot) for comparison, including panel A Venn counts",supp file 17,Yes,Yes,
,mapped_dm6.sorted.bam,Reads from APPLd flies mapped to dm6 for coverage and deletion analysis,supp file 18,No,Yes,
,mapped_dm6.sorted.bam.bai,BAM index required for IGV/jBrowse viewing,supp file 19,Yes,Yes,
folder 4,pseudoX_singleline.fa,Fasta of PseudoX reconstruction from ptg000067l + 20nt + ptg000024l RC,supp file 20,Yes,Yes,Figure 3
,final_genesblast_bed.bed,"BED file of filtered BLAST hits for elav, Appl, and vnd to PseudoX",supp file 21,Yes,Yes,
,high_coverage_only.gtf,Liftoff annotation (coverage ≥0.95) for genes mapped to PseudoX,supp file 22,Yes,Yes,
,combined_dgenies_alignments/,PDF/PNG alignments of PseudoX and full genome to dm6 chrX (from D-GENIES),supp file 23,No,Yes,
,blastn_appl_vs_pseudoX.tsv,Raw BLAST output of Appl CDS to PseudoX,supp file 24,Yes,Yes,
,contig_orient_mapping/ directory,"Folder containing all eight pseudo-X orientation assemblies and their mapped read BAMs, edge read extraction, and testing of potential bridging contigs (004, 081). Demonstrates all tested configurations and mapping evidence used to rule out alternate orientations of the Appl locus in APPLd.",supp file 25,No,Yes,
,appl_vs_appld_on_appl.bed,File of contig hits mapped on the Appl locus for visualization on Flybase Appl CDS,supp file 26,Yes,Yes,
folder 5,T2A_stainfree_gel.jpg,Image of stainfree gel for T2A blot normalization ,supp file 27,Yes,Yes,Figure 4
,T2A_AP2_blot.jpg,Image of membrane for T2A blot using chicken anti-APPL (Philip Copenhaver),supp file 28,Yes,Yes,
,Appl-GFP_blot_quant.png,Image of Appl-GFP blot and quantification in figure format (individual files available on request),supp file 29,Yes,Yes,
folder 6,timsTOF_data_allgroups.xlsx,timsTOF data from all comparisons,supp file 32,Yes,Yes,Figure 5
,APPLdvW1118cytoscape.cys,APPLd vs W1118 cytoscape PPI from Metascape,supp file 33,Yes,Yes,
,APPLdvW1118_cytoscape_clusters.xlsx,APPLd vs W1118 pathways from Metascape,supp file 34,Yes,Yes,
,AppldVsW1118_volcano,Prism file containing information from APPLd to W1118 comparison of timsTOF,supp file 35,Yes,Yes,
,T2AvCantoncytoscape.cys,Appl-T2A-GAL4 v CantonS cytoscape PPI from Metascape,supp file 36,Yes,Yes,
,T2AvCantonS_cytoscape_clusters.xlsx,Appl-T2A-GAL4 v CantonS pathways from Metascape,supp file 37,Yes,Yes,
,T2AVsCanton_volcano,Prism file containing information from Appl-T2A-GAL4 to CantonS comparison of timsTOF,supp file 38,Yes,Yes,
,APPLdvT2Acytoscape.cys,APPLd v Appl-T2A-GAL4 cytoscape PPI from Metascape,supp file 39,Yes,Yes,
,APPLdvT2A_cytoscape_clusters.xlsx,APPLd v Appl-T2A-GAL4 pathways from Metascape,supp file 40,Yes,Yes,
,AppldVsT2A_volcano,Prism file containing information from APPLd to Appl-T2A-GAL4 comparison of timsTOF,supp file 41,Yes,Yes,
folder 7,Appl_signature_table.xlsx,List of proteins found to be significant in T2A/CantonS and APPLd/W1118 but not in APPLd/T2A,supp file 42,Yes,Yes,Figure 5 and Figure 6
,APPLd_signature_table.xlsx,List of proteins found to be significant in APPLd/T2A and APPLd/W1118 but not in CantonS/W1118,supp file 43,Yes,Yes,
,CG3156 CDS.prot,Normal predicted protein sequence of CG3156,supp file 44,Yes,Yes,Figure 6
,CG3156 altered CDS.prot,Altered predicted protein sequence of CG3156 in APPLd,supp file 45,Yes,Yes,
,Pdzd8 CDS.prot,Normal predicted protein sequence of Pdzd8,supp file 46,Yes,Yes,
,Pdzd8 altered CDS.prot,Altered predicted protein sequence of Pdzd8 in APPLd,supp file 47,Yes,Yes,
folder 8,orbitrap APPLd Volcano Prism Data.xlsx,"Processed protein-level differential expression results used for plotting the volcano plot in Figure 4B. Includes fold change, adjusted P-values, and filtering.",supp file 48,Yes,Yes,Supplemental Figure 3
,orbitrap GO_MCODE.csv,Enrichment summary for MCODE clusters generated in Metascape. Contains functional GO term labels shown in Figure 4A.,supp file 49,Yes,Yes,
,orbitrap metascape_scaling_values.xlsx,Node metadata for Cytoscape PPI map — includes matched fold change and significance values for scaling node size and color.,supp file 50,Yes,Yes,
,orbitrap metascapefinalAPPLd_BACKGROUND.xlsx,"List of all background proteins detected and used in the Metascape analysis, used for enrichment statistics.",supp file 51,Yes,Yes,
,orbitrap metascapefinalAPPLd_SIG.xlsx,Filtered list of significantly differentially expressed proteins (adj. P < 0.05) input to Metascape for PPI and clustering.,supp file 52,Yes,Yes,
folder 9,MetaComparisonProteomics.xlsx,"Raw and processed comparison between the current 3–5 day APPLd proteomics dataset and the significant gene list from Nithianandam et al., 2023 (~10d APPLd flies). Includes gene symbols, log2FCs from both studies, difference values, and classification into expression categories.",supp file 53,Yes,Yes,Supplemental Figure 7
folder 10,APPLd_signature_primers.xlsx,Primers uses for verification of APPLd predicted mutations from proteomic signature and genome assembly for use in sanger sequencing.,supp file 54,Yes,Yes,Supplemental Figure 6
folder 11,SleepAnticipationScript.R,Custom R pipeline for calculating Morning (MAI) and Evening (EAI) Anticipation Indices,supp file 55,Yes,Yes,Figure 7 and Supplemental Figure 8
,ClimbingLmmScript.R,R script for Linear Mixed-Model (LMM) analysis of negative geotaxis velocity,supp file 56,Yes,Yes,
