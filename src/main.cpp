#include <stdio.h>
#include <stdlib.h>
#include "svat_utils.h"
#include "svat_aggregation_utils.h"
#include "svat_ansi_string.h"

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
Variant/Genotype Encryption:\n\
	Variant Encryption:\n\
		-encrypt_vectorized_SNVs\n\
	Genotype Matrix (Re)Encryption\n\
		-encrypt_encoded_SNV_genotypes_2_matrix\n\
		-reencrypt_genotype_matrix\n\
Variant Simulation:\n\
	-generate_simulated_InDels_on_EOI\n\
	-generate_simulated_SNVs_on_EOI\n\
VEP Track Generation:\n\
	-get_EOIs_BED_per_exonic_regions\n\
	-get_frame_gene_unique_EOIs_BED_per_CDS_regions\n\
	-generate_EOI_enumerating_SNV_VCF\n\
	-generate_EOI_enumerating_Indel_VCF\n\
Vectorization:\n\
	VEP Annotated Variants:\n\
		-separate_VEP_output_per_chromosome\n\
		-signalize_VEP_annotated_SNVs\n\
		-signalize_VEP_annotated_Inserts\n\
		-signalize_VEP_annotated_Deletes\n\
		-signalize_VEP_annotated_Deletes_w_CDS_vicinity\n\
	VCF Variants:\n\
		-separate_VCF_per_chromosome\n\
		-signalize_SNV_genotypes_per_VCF\n\
		-signalize_Insertion_genotypes_per_VCF\n\
		-signalize_Deletion_genotypes_per_VCF\n\
Annotation Signal Processing:\n\
	Plaintext:\n\
		-multiply_signalized_SNVs_and_annotations\n\
		-multiply_signalized_Indels_and_annotations\n\
	Secure:\n\
		-secure_multiply_signalized_SNVs_and_annotations\n\
		-secure_multiply_signalized_Indels_and_annotations\n\
Variant Genotype Matrix Aggregation:\n\
	Plaintext Variant Encoding:\n\
		-encode_SNV_genotypes_2_matrix_per_VCF\n\
		-encode_DEL_genotypes_2_matrix_per_VCF\n\
	Plaintext Aggregation:\n\
		-plain_aggregate_encoded_Del_genotype_matrix\n\
		-plain_aggregate_encoded_SNV_genotype_matrix\n\
	Secure Aggregation:\n\
		-secure_aggregate_encrpyted_SNV_genotype_matrix\n\
Plain text Impact Analysis:\n\
	-compare_VEP_annotations_with_SVAT_annotations\n\
	-extract_annotated_SNVs_from_signal\n\
	-translate_annotated_Deletions_from_annotated_signals\n\
	-translate_annotated_Insertions_from_annotated_signals\n\
	-translate_annotated_SNVs_from_annotated_signals\n\
----------------------------------------------------------------\n\n", argv[0]);

		exit(0);
	}
	
	if (t_string::compare_strings(argv[1], "-secure_multiply_signalized_SNVs_and_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Plaintext Annotation Vector directory] \
[Encrypted vectorized SNVs directory] \
[Encrypted Annotated Variants Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* plaintext_annotation_signal_dir = argv[2];
		char* encrypted_SNV_vector_dir = argv[3];
		char* op_dir = argv[4];

		secure_multiply_SNV_variant_and_annotation_signals(plaintext_annotation_signal_dir, encrypted_SNV_vector_dir, op_dir);
	} // -secure_multiply_signalized_SNVs_and_annotations
	else if (t_string::compare_strings(argv[1], "-secure_aggregate_encrypted_SNV_genotype_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encrypted SNV genotypes directory] \
[EOI regions BED path] \
[Variant regions BED path] \
[Sample id's list file path] \
[Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encrypted_vectorized_snv_genotypes_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* VOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		secure_aggregate_encrypted_SNV_genotype_matrix(encrypted_vectorized_snv_genotypes_dir,
			EOI_regs_BED_fp,
			VOI_regs_BED_fp,
			sample_ids_list_fp,
			genotype_encoding_type);
	} // secure_aggregate_encrpyted_SNV_genotype_matrix option.
	else if (t_string::compare_strings(argv[1], "-encrypt_vectorized_SNVs"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Vectorized variants directory] \
[EOI regions BED path] \
[Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vectorized_SNVs_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* encrypted_SNVs_vector_op_dir = argv[4];

		encrypt_vectorized_SNVs(vectorized_SNVs_dir, EOI_regs_BED_fp, encrypted_SNVs_vector_op_dir);
	} // -encrypt_encoded_SNV_genotypes_2_matrix option end.
	else if (t_string::compare_strings(argv[1], "-encrypt_encoded_SNV_genotypes_2_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encoded SNV genotypes directory] \
[Encrypted SNV genotypes output directory] \
[EOI regions BED path] \
[Sample id's list file path] \
[Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encoded_vectorized_snvs_dir = argv[2];
		char* encrypted_snv_genotype_matrix_dir = argv[3]; // This is the output directory.
		char* EOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		encrypt_encoded_SNV_genotypes_2_matrix(encoded_vectorized_snvs_dir,
			encrypted_snv_genotype_matrix_dir,
			EOI_regs_BED_fp,
			sample_ids_list_fp, // We need this only to know the number of inidividuals in the dataset, nothing else.
			genotype_encoding_type);
	} // -encrypt_encoded_SNV_genotypes_2_matrix option end.
	else if (t_string::compare_strings(argv[1], "-reencrypt_SNV_genotypes_matrix"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Encrypted SNV genotypes directory] \
[Re-Encrypted SNV genotypes output directory] \
[EOI regions BED path] \
[Sample id's list file path] \
[Genotype Encoding (0:Existence/1:Allele Count)] \
[Public keys file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encrypted_snv_genotype_matrix_dir = argv[2];
		char* reencrypted_snv_genotype_matrix_dir = argv[3];
		char* EOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);
		char* public_keys_fp = argv[7];

		re_encrypt_genotype_matrix(encrypted_snv_genotype_matrix_dir, // This is the output directory and the reencrypted data.
			reencrypted_snv_genotype_matrix_dir, // This is the output directory and the reencrypted data.
			EOI_regs_BED_fp, // This is necessary only to know the length of vectorized representation.
			sample_ids_list_fp, 
			genotype_encoding_type,
			public_keys_fp);
	} // -reencrypt_SNV_genotypes_matrix option.
	if (t_string::compare_strings(argv[1], "-plain_aggregate_encoded_SNV_genotype_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encoded genotypes directory] [EOI regions BED path] [VOI regions BED path] [Sample id's list file path] [Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encoded_vectorized_snv_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* VOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		plain_aggregate_encoded_SNV_genotype_Full_Matrix(encoded_vectorized_snv_dir,
			EOI_regs_BED_fp,
			VOI_regs_BED_fp,
			sample_ids_list_fp,
			genotype_encoding_type);
	} // -plain_aggregate_encoded_SNV_genotype_matrix option.
	else if (t_string::compare_strings(argv[1], "-plain_aggregate_encoded_Del_genotype_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encoded genotypes directory] [EOI regions BED path] [VOI regions BED path] [Sample id's list file path] [Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encoded_vectorized_del_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* VOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		plain_aggregate_encoded_Del_genotype_matrix(encoded_vectorized_del_dir,
													EOI_regs_BED_fp,
													VOI_regs_BED_fp,
													sample_ids_list_fp,
													genotype_encoding_type);
	} // -plain_aggregate_encoded_Del_genotype_matrix option.
	if (t_string::compare_strings(argv[1], "-encode_DEL_genotypes_2_matrix_per_VCF"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Per chrom VCF dir.] [Sample ids list path] [EOI regions BED path] [Genotype Encoding (0:Existence/1:Allele Count)] [Genotypes output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* per_chrom_VCF_dir = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* EOI_regs_BED_fp = argv[4];
		int genotype_encoding_type = atoi(argv[5]);
		char* genotypes_op_dir = argv[6];

		encode_DEL_genotypes_2_matrix_per_VCF(per_chrom_VCF_dir, sample_ids_list_fp,
			EOI_regs_BED_fp,
			genotype_encoding_type, genotypes_op_dir);
	} // -compare_VEP_annotations_with_SVAT_annotations option.
	else if (t_string::compare_strings(argv[1], "-encode_SNV_genotypes_2_matrix_per_VCF"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Per chrom VCF dir.] [Sample ids list path] [EOI regions BED path] [Genotype Encoding (0:Existence/1:Allele Count)] [Genotypes output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* per_chrom_VCF_dir = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* EOI_regs_BED_fp = argv[4];
		int genotype_encoding_type = atoi(argv[5]);
		char* genotypes_op_dir = argv[6];

		encode_SNV_genotypes_2_matrix_per_VCF(per_chrom_VCF_dir, sample_ids_list_fp,
			EOI_regs_BED_fp,
			genotype_encoding_type, genotypes_op_dir);
	} // -compare_VEP_annotations_with_SVAT_annotations option.
	else if (t_string::compare_strings(argv[1], "-compare_VEP_annotations_with_SVAT_annotations"))
	{ 
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [VEP annotated output file path] [VEP impact strings list (full)] \
[VEP impacts 2 focus list file path] [High impacts list file path] [SVAT annotated BED file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vep_annotated_op_fp = argv[2];
		char* VEP_impact_string_ctx_fp = argv[3];
		char* vep_impacts_2_focus_list_fp = argv[4];
		char* hi_impacts_list_fp = argv[5];
		char* svat_annotated_vcf_fp = argv[6];

		compare_VEP_annotations_with_SVAT_annotations(vep_annotated_op_fp,
														VEP_impact_string_ctx_fp,
														vep_impacts_2_focus_list_fp,
														hi_impacts_list_fp,
														svat_annotated_vcf_fp);
	} // -compare_signal_annotated_vs_direct_VEP_annotations option.
	else if (t_string::compare_strings(argv[1], "-separate_VCF_per_chromosome"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [VCF file path (stdin ok)] [chromosome ids/lengths list file path] [Per chromosome output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vcf_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* per_chrom_op_dir = argv[4];

		separate_VCF_per_chromosome(vcf_fp, chr_ids_lengths_list_fp, per_chrom_op_dir);
	} // -separate_VEP_output_per_chromosome option.
	else if (t_string::compare_strings(argv[1], "-separate_VEP_output_per_chromosome"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [VEP output file path (stdin ok)] [chromosome ids/lengths list file path] [Per chromosome output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vep_op_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* per_chrom_op_dir = argv[4];

		separate_VEP_output_per_chromosome(vep_op_fp, chr_ids_lengths_list_fp, per_chrom_op_dir);
	} // -separate_VEP_output_per_chromosome option.
	else if (t_string::compare_strings(argv[1], "-generate_simulated_InDels_on_EOI"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Maximum indel length] [Per position variant Probability] [Deletion Probability (<1.0)] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		int l_max_indel = atoi(argv[3]);
		double per_posn_var_prob = atof(argv[4]);
		double del_prob = atof(argv[5]);
		char* bin_seq_dir = argv[6];
		char* op_vcf_fp = argv[7];

		generate_EOI_randomized_Indel_VCF(EOI_bed_fp, l_max_indel, per_posn_var_prob, del_prob, bin_seq_dir, op_vcf_fp);
	} // -generate_simulated_InDels_on_EOI option.
	else if (t_string::compare_strings(argv[1], "-generate_simulated_SNVs_on_EOI"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Per position variant Probability] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		double per_posn_var_prob = atof(argv[3]);
		char* bin_seq_dir = argv[4];
		char* op_vcf_fp = argv[5];

		generate_EOI_randomized_SNV_VCF(EOI_bed_fp, per_posn_var_prob, bin_seq_dir, op_vcf_fp);
	} // -generate_simulated_InDels_on_EOI option.
	else if (t_string::compare_strings(argv[1], "-extract_annotated_SNVs_from_signal"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s %s [Annotated SNV signals directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* op_fp = argv[3];

		extract_annotated_SNVs_from_signal(annotation_signal_dir, op_fp);
	} // -extract_annotated_SNVs_from_signal option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_Deletions_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_Deletions_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Deletions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_Insertions_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_Insertions_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Insertions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_SNVs_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_SNVs_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Insertions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-multiply_signalized_SNVs_and_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotation signal directory] [Variant signal directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* variant_signal_dir = argv[3];
		char* op_dir = argv[4];

		multiply_SNV_variant_and_annotation_signals(annotation_signal_dir, variant_signal_dir, op_dir);
	} // -multiply_signalized_SNVs_and_annotations option.
	else if (t_string::compare_strings(argv[1], "-multiply_signalized_Indels_and_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotation signal directory] [Variant signal directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* variant_signal_dir = argv[3];
		char* op_dir = argv[4];

		multiply_Indel_variant_and_annotation_signals(annotation_signal_dir, variant_signal_dir, op_dir);
	} // -multiply_signalized_Indels_and_annotations option.
	else if (t_string::compare_strings(argv[1], "-get_frame_gene_unique_EOIs_BED_per_CDS_regions"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Element IDs list file path] [Sub-element selector: \"CDS\"/\"exon\"] [GFF file path] [Exon extension length (5 bps)] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* element_ids_list_fp = argv[2];
		char* CDS_exon_selector = argv[3];
		char* gff_fp = argv[4];
		int l_exon_ext = atoi(argv[5]);
		char* op_fp = argv[6];

		get_frame_gene_unique_EOIs_BED_per_CDS_regions(element_ids_list_fp, CDS_exon_selector, gff_fp, l_exon_ext, op_fp);
	}
	else if (t_string::compare_strings(argv[1], "-get_EOIs_BED_per_exonic_regions"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Element IDs list file path] [GFF file path] [CDS/Exon selection: \"CDS\"/\"exon\"] [Element type (\"Gene\"/\"Transcript\")] [Exon extension length (5 bps)] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* element_ids_list_fp = argv[2];
		char* gff_fp = argv[3];
		char* CDS_exon_selector = argv[4];
		char* EOI_element_type_str = argv[5];
		int l_exon_ext = atoi(argv[6]);
		char* op_fp = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "EOIs with genes.\n");
			EOI_element_type = EOI_GENE_ELEMENTS;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "EOIs with transcripts.\n");
			EOI_element_type = EOI_TRANSCRIPT_ELEMENTS;
		}

		get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF(element_ids_list_fp, gff_fp, CDS_exon_selector, EOI_element_type, l_exon_ext, op_fp);
	} // -get_EOIs_BED_per_exonic_regions
	else if (t_string::compare_strings(argv[1], "-generate_EOI_enumerating_SNV_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		char* bin_seq_dir = argv[3];
		char* op_fp = argv[4];

		generate_EOI_enumerating_SNV_VCF(EOI_bed_fp, bin_seq_dir, op_fp);
	} // -generate_Exome_enumerating_VCF option.
	else if (t_string::compare_strings(argv[1], "-generate_EOI_enumerating_Indel_VCF"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Maximum indel length] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		int l_max_indel = atoi(argv[3]);
		char* bin_seq_dir = argv[4];
		char* op_fp = argv[5];

		generate_EOI_enumerating_Indel_VCF(EOI_bed_fp, l_max_indel, bin_seq_dir, op_fp);
	} // -generate_Exome_enumerating_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_Insertion_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [VCF file path (stdin not ok)] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* VCF_fp = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_Insertion_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			VCF_fp,
			op_dir);
	} // -signalize_Indel_genotypes_per_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_Deletion_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Per chromosome VCF directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_Deletion_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			per_chrom_VCF_dir,
			op_dir);
	} // -signalize_Indel_genotypes_per_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_SNV_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [VCF file path (stdin not ok)] [Output directory]\n", argv[0], argv[1]);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* VCF_fp = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_SNV_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			VCF_fp,
			op_dir);
	} // -signalize_VCF_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_SNVs"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP Output file path (stdin not ok)] \
[Genome sequence directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_fp = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "DO NOT USE SUMMARIZATION OVER GENES, THIS OPTION IS NOT TESTED WELL.\n");exit(0);
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_SNVs_per_EOI_regs_element_summarization(EOI_regs_BED_fp, 
																		genome_dir, 
																		EOI_element_type, 
																		vep_op_fp,
																		sorted_impact_value_strings_list_fp,
																		op_dir);
	} // -signalize_VEP_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Inserts"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP Output file path (stdin not ok)] \
[Genome sequence directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_fp = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Inserts_per_EOI_regs_element_summarization(EOI_regs_BED_fp,
			genome_dir, 
			EOI_element_type,
			vep_op_fp,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.



	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Deletes"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP per chromosome output directory] \
[Per chromosome sequences directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_dir = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization(EOI_regs_BED_fp,
			genome_dir,
			EOI_element_type,
			vep_op_dir,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Deletes_w_CDS_vicinity"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP per chromosome output directory] \
[Per chromosome sequences directory] [CDS BED file path] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_dir = argv[4];
		char* genome_dir = argv[5];
		char* cds_bed_fp = argv[6];
		char* sorted_impact_value_strings_list_fp = argv[7];
		char* op_dir = argv[8];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization_with_CDS_vicinity(EOI_regs_BED_fp,
			genome_dir, cds_bed_fp,
			EOI_element_type,
			vep_op_dir,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.



}
