#ifndef __SVAT_ENCRYPTION_UTILS__
#define __SVAT_ENCRYPTION_UTILS__

enum {GENOTYPE_ENC_EXISTENCE, GENOTYPE_ENC_GENO, GENOTYPE_ENC_HAPLO};

void encode_DEL_genotypes_2_matrix_per_VCF(char* per_chrom_VCF_dir, char* sample_ids_list_fp,
	char* EOI_regs_BED_fp,
	int genotype_encoding_type, char* genotypes_op_dir);

void encode_SNV_genotypes_2_matrix_per_VCF(char* per_chrom_VCF_dir, char* sample_ids_list_fp,
	char* EOI_regs_BED_fp,
	int genotype_encoding_type, char* genotypes_op_dir);

bool sort_encoded_geno_per_track_posn(void** encoded_geno_entry1, void** encoded_geno_entry2);

void encrypt_encoded_SNV_genotypes_2_matrix(char* encoded_vectorized_snvs_dir,
	char* encrypted_snv_genotype_matrix_dir,
	char* EOI_regs_BED_fp,
	char* sample_ids_list_fp, 
	int genotype_encoding_type);

void re_encrypt_genotype_matrix(char* encrypted_genotype_matrix_directory_list,
	char* pooled_encrypted_genotype_matrix_dir,
	char* EOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type,
	char* public_keys_fp);

unsigned int* secure_aggregate_encrypted_SNV_genotype_matrix(char* encrypted_vectorized_snv_genotypes_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type);

unsigned int* plain_aggregate_encoded_Del_genotype_matrix(char* encoded_vectorized_del_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type);

unsigned int* plain_aggregate_encoded_SNV_genotype_Full_Matrix(char* encoded_vectorized_snvs_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type);

unsigned int* plain_aggregate_encoded_SNV_genotype_matrix(char* encoded_vectorized_del_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type);

#endif // __SVAT_ENCRYPTION_UTILS__
