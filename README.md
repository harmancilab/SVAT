# SVAT: Secure Variant Annotation and Aggregation Tool

## Contents 

* [Introduction](#introduction)
* [Build](#build)
* [Usage](#usage)
    * [Target Region Generation](#target-region-generation)
    * [Variant Annotation](#variant-annotation)    
    * [Genotype Matrix Aggregation](#genotype-matrix-aggregation)      


This repository contains the source code and documentation for SVAT, a computer program that can be used for secure outsourcing of 2 basic tasks: 

1. **Variant Annotation:** Given a list of small variant loci (single nucleotide variants and small insertion/deletions -- Indels) with alternate alleles, SVAT can be used to annotate the variants with respect to their impact on the protein-coding sequences.
2. **Allele Frequency Aggregation:** This is basically secure pooling and counting of variant allele frequencies. SVAT can perform counting at the allele count (i.e., allele frequency) or  variant existence level (similar to beacons).

SVAT provides a general framework for storing the sensitive variant and genotype information securely in an encrypted format. SVAT makes use of homomorphic encryption (HE) security framework. HE has numerous advantages that make it very appealing for protecting data:
1. The data stays encrypted at-rest, in-transit, and even while it is being analyzed. In other words, data does not ever need to be decrypted.
2. The encrypted data is secure (including quantum attacks) even if it is stolen/hacked by unauthorized entities.
3. Encrypted data can be used on a cloud with no repercussions. Thus, HE frameworks mix very well with the cloud's massive computational and storage resources.

A major limitation of HE is that it can be notoriously slow for certain computations. **Thanks to recent advancements in HE literature**, machine learning models can be securely run at speeds comparable to computing on the data in plain format. 

SVAT provides a framework that combines an efficient data representation with ultra-fast secure HE evaluation schemes to provide high-speed aggregation/annotation operations.

## Introduction
SVAT is developed as 2 modules that perform parts of the annotation and aggregation tasks. 

**First module (SVAT)** implements the operations for preprocessing/postprocessing:
1. Preprocessing of the variant/genotype data for annotation and aggregation.
2. Postprocessing of the data after it is received.
2. Variant vectorization procedure that is necessary for the efficient and secure analysis of variant positions.

**Second module (HESVAT)** is the security module. This module performs:
1. Public/Private Key generation [Data Owner]
2. Encryption/Re-encryption/Decryption of variant and genotype data [Data Owner/Server]
3. Secure annotation/aggregation operations [Server]

## Build

You can download the main project with all contained submodule to your local computer by using the clone command:

```
git clone --recurse-submodules https://github.com/harmancilab/SVAT.git
```

To build SVAT, run following commands on the base directory:

```
make clean
make
cd HESVAT
make clean
make
```

These commands build the SVAT and HESVAT executables.

## Usage
We provide step-by-step examples of running SVAT for secure annotation and aggregation.

### Target Region Generation
The first step in annotation and aggregation tasks is the generation of the target regions.

```
# Download GFF formatted gene annotations. We use gencode v31.
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz

# Extract the gene id's.
gzip -cd gencode.v31.annotation.gff3.gz | grep "gene_type=protein_coding" | awk {'if($3=="exon" && $1=="chr1")print $0'} | awk 'BEGIN{FS="\t"}{split($9, arr, ";");for(i=1;i<=length(arr);i++){if(substr(arr[i], 1, 7)=="Parent="){print substr(arr[i], 8);;break;}}}' | sort -u > transcript_ids.list

# Generate the sorted target regions.
gzip -cd gencode.v31.annotation.gff3.gz | grep "gene_type=protein_coding" | awk -v start_pos=1 -v end_pos=250000000 {'if($1=="chr1" && $4>start_pos && $4<end_pos){print $0}'} | SVAT -get_EOIs_BED_per_exonic_regions transcript_ids.list stdin Transcript ${l_exon_extension} EOI_gencode31_transcript_exons_merged.bed
```

### Variant Annotation

The variant annotation takes a VCF file as input. SVAT can generate simulated VCF files that can be used for testing and comparison. A VCF file can be simulated using:

```
# Generate simulated variants.
HG38_DIR=hg38_genome

# First, download and process the genome sequence.
mkdir $HG38_DIR
cd $HG38_DIR
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
SVAT -binarize_fasta_per_file hg38.fa.gz .
cd ..

# Now we can simulate.
EOI_exons_fp=EOI_gencode31_transcript_exons_merged.bed
SVAT -generate_simulated_InDels_on_EOI $EOI_exons_fp 10 $mutation_prob 1.999 ${HG38_DIR} simulated_dels.vcf.gz
```

Above command simulates deletion variants.

1. **VCF Conversion:** First, SVAT processes VCF file and separates into chromosomes.

```
# Separate VCF into chromosomes.
rm -f -r per_chrom_VCF/
mkdir per_chrom_VCF/
SVAT -separate_VCF_per_chromosome simulated_dels.vcf.gz hg38.list per_chrom_VCF
```

2. **Vectorize the VCF file on target regions:**
Finally, the deletion variants are vectorized.
```
EOI_exons_fp=EOI_gencode31_transcript_exons_merged.bed
SVAT -signalize_Deletion_genotypes_per_VCF $EOI_exons_fp per_chrom_VCF VCF_signal
```

3. **Encryption and Secure Annotation and Decryption of the Vectorized Variants (HESVAT)**
The variants are encrypted and secured using the commands under HESVAT's annotation task (task1): [Secure Annotation](https://github.com/K-miran/HESVAT/tree/808fc23#Task-1-Secure-Annotation)

4. **Translation of the Annotated Variants Vector**
Given the decrypted signal is stored under 'Annotated_Signal', the deletions can be translated using the following command:
```
SVAT -translate_annotated_Deletions_from_annotated_signals Annotated_Signal per_chrom_VCF annotated_variants.vcf
```

### Genotype Matrix Aggregation
The aggregation process starts with a genotype matrix. It is necessary to separate the genotype matrix into chromosomes. It should be noted that in many cases the genotype matrices are by-default chromosome separated and this step can be skipped.

1. **Genotype encoding and coordinate vectorization:** Following separation of VCF per chromosomes, SVAT re-encodes and generates a genotype matrix in vectorized coordinates.

```
encoding=0 # Set to 1 for counting alleles.
SVAT -encode_SNV_genotypes_2_matrix_per_VCF per_chrom_VCF sample_ids.list EOI_gencode19_gene_exons.bed_merged.bed ${encoding} Encoded_Genotypes
```

2. **Encryption (Re-encryption), Secure Aggregation and Decryption of the Genotype Matrix**: HESVAT is used to encrypt the encoded genotype matrix and securely aggregate it, as described here [Secure Aggregation](https://github.com/K-miran/HESVAT/tree/808fc23#Task-2-Secure-Aggregation)
HESVAT can also perform re-encryption and aggregation of the matrices. The re-encryption is described here [Genotype re-encryption](https://github.com/K-miran/HESVAT/tree/808fc23#Task-3-Secure-Aggregation-by-proxy-encryption)
