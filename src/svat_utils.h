#ifndef __CRYPTANNOT__
#define __CRYPTANNOT__

#include <vector>
using namespace std;

enum {
	CDS_ctx,
	splice_region_ctx,
	splice_donor_ctx,
	splice_acceptor_ctx,
	intron_ctx,
	fp_UTR_ctx,
	tp_UTR_ctx,
	start_codon_ctx,
	stop_codon_ctx,
	N_VEP_CTX
};

struct t_VEP_term_ctx
{
	char* id;
	bool* context_array;

	bool CDS_ctx;
	bool splice_region_ctx;
	bool splice_donor_ctx;
	bool splice_acceptor_ctx;
	bool intron_ctx;
	bool fp_UTR_ctx;
	bool tp_UTR_ctx;
	bool start_codon_ctx;
	bool stop_codon_ctx;
};

void unpack_impact_signal_values(unsigned long long packed_impact_signal_value,
	unsigned int& coding_frame,
	unsigned long long& impact_value,
	char* neigh_seq_nucs, int n_bits_per_nuc, int n_neigh_nucs);

void unpack_impact_signal_values_w_CDS_neighborhood(unsigned long long packed_impact_signal_value,
	unsigned int& coding_frame,
	unsigned long long& impact_value,
	char* geno_neigh_seq_nucs,
	char* CDS_neigh_seq_nucs, int n_bits_per_nuc, int n_neigh_nucs);

unsigned long long pack_set_impact_signal_values_w_CDS_neighborhood(unsigned int coding_frame,
	unsigned long long impact_value,
	unsigned int geno_neigh_seq_signal,
	unsigned int cds_neigh_seq_signal,
	int n_bits_per_nuc, int n_neigh_nucs);

unsigned long long pack_set_impact_signal_values(unsigned int coding_frame,
	unsigned long long impact_value,
	unsigned int neigh_seq_signal, int n_bits_per_nuc, int n_neigh_nucs);

void unpack_neigh_nucs(unsigned int neigh_seq_signal, char* buff, int n_bits_per_nuc, int n_nucs);

//void compare_signal_annotated_vs_direct_VEP_annotations(char* vcf_fp,
//	char* VEP_impact_string_ctx_fp,
//	char* signal_annotated_VCF_BED_fp,
//	char* vep_annotated_fp);

bool check_element_context_per_impact_bitmap(vector<t_VEP_term_ctx*>* vep_term_ctx, int context_id, unsigned long long impact_entry_bitmap);

void compare_VEP_annotations_with_SVAT_annotations(char* vep_annotated_op_fp,
	char* VEP_impact_string_ctx_fp,
	char* vep_impacts_2_focus_list_fp,
	char* high_impacts_list_fp,
	char* svat_annotated_vcf_fp);

vector<t_VEP_term_ctx*>* load_VEP_annotation_term_context(char* vep_annotation_term_ctx_fp, vector<char*>* sorted_vep_terms);

void separate_VEP_output_per_chromosome(char* vep_op_fp, char* chr_ids_lengths_list_fp, char* per_chrom_op_dir);
void separate_VCF_per_chromosome(char* vcf_fp, char* chr_ids_lengths_list_fp, char* per_chrom_op_dir);

enum { VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION, VEP_SIGNALIZE_GENE_SUMMARIZATION };

void signalize_VEP_annotated_SNVs_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
	char* genome_dir,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir);

void signalize_VCF_Deletion_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
	char* per_chrom_VCF_dir,
	char* op_dir);

void signalize_VCF_Insertion_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
	char* per_chrom_VCF_dir,
	char* op_dir);

void signalize_VCF_SNV_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
	char* per_chrom_VCF_dir,
	char* op_dir);

void encrypt_vectorized_SNVs(char* vectorized_SNVs_dir, char* EOI_regs_BED_fp, char* encrypted_SNVs_vector_op_dir);

void secure_multiply_SNV_variant_and_annotation_signals(char* plaintext_annotation_signal_dir, char* encrypted_SNV_vector_dir, char* op_dir);

void generate_EOI_randomized_Indel_VCF(char* EOI_bed_fp, int l_max_indel, double per_posn_var_prob, double del_prob, char* bin_seq_dir, char* op_vcf_fp);
void generate_EOI_randomized_SNV_VCF(char* EOI_bed_fp, double per_posn_var_prob, char* bin_seq_dir, char* op_vcf_fp);

void multiply_SNV_variant_and_annotation_signals(char* annotation_signal_dir, char* variant_signal_dir, char* op_dir);

void multiply_Indel_variant_and_annotation_signals(char* annotation_signal_dir, char* variant_signal_dir, char* op_dir);

void extract_annotated_SNVs_from_signal(char* annotation_signal_dir, char* op_fp);
//void extract_annotated_Indels_from_signal(char* annotated_variant_signal_dir, char* op_fp);
void translate_annotated_Deletions_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp);
void translate_annotated_Insertions_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp);
void translate_annotated_SNVs_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp);

enum { EOI_GENE_ELEMENTS, EOI_TRANSCRIPT_ELEMENTS };
void get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF(char* element_ids_list_fp, char* GENCODE_gff_fp, char* CDS_exon_selector, int eoi_element_type, int l_exon_ext, char* op_fp);
void get_frame_gene_unique_EOIs_BED_per_CDS_regions(char* gene_ids_list_fp, char* CDS_exon_selector, char* GENCODE_gff_fp, int l_exon_ext, char* op_fp);

void generate_EOI_enumerating_SNV_VCF(char* EOI_bed_fp, char* bin_seq_dir, char* op_fp);
void generate_EOI_enumerating_Indel_VCF(char* EOI_bed_fp, int l_max_indel, char* bin_seq_dir, char* op_vcf_fp);

void signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
	char* genome_dir,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir);

void signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization_with_CDS_vicinity(char* EOI_regs_BED_fp,
	char* genome_dir,
	char* CDS_regs_fp,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir);

void signalize_VEP_annotated_Inserts_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
	char* genome_dir,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir);

#endif // __CRYPTANNOT__