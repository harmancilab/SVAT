#include <stdio.h>
#include <stdlib.h>
#include "file_utils.h"
#include "svat_ansi_string.h"
#include "svat_rng.h"
#include "svat_seed_manager.h"
#include "svat_annot_region_tools.h"
#include "svat_genomics_coords.h"
#include "svat_genome_sequence_tools.h"
#include "svat_nucleotide.h"
#include "svat_utils.h"
#include "svat_aggregation_utils.h"
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <algorithm>
#include <vector>
using namespace std;

bool __DUMP_CRYPTAGGREGATE_INDEL_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

///////////////////////////////////////////////////////////////////////////////////////////////
// Encoding functions: Converts the VCF file genotypes to 1 or 2 bit per individual matrices.
///////////////////////////////////////////////////////////////////////////////////////////////
void encode_DEL_genotypes_2_matrix_per_VCF(char* per_chrom_VCF_dir, char* sample_ids_list_fp,
	char* EOI_regs_BED_fp,
	int genotype_encoding_type, char* genotypes_op_dir)
{
	char* genotype_encoding_strs[3];
	genotype_encoding_strs[GENOTYPE_ENC_EXISTENCE] = t_string::copy_me_str("GENO_ENC_EXISTENCE");
	genotype_encoding_strs[GENOTYPE_ENC_GENO] = t_string::copy_me_str("GENO_ENC_GENO");
	genotype_encoding_strs[GENOTYPE_ENC_HAPLO] = t_string::copy_me_str("GENO_ENC_HAPLO");

	fprintf(stderr, "Encoding deletion genotypes using genotype encoding %s\n", genotype_encoding_strs[genotype_encoding_type]);

	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Encoding SNV genotype matrix on %d regions.\n", EOI_regs->size());

	vector<char*>* chr_ids = get_chr_ids(EOI_regs);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Using %d samples in VCF data under %s\n", sample_ids->size(), per_chrom_VCF_dir);

	int n_bits_per_sample = 2;
	if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
	{
		n_bits_per_sample = 1;
	}

	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * (sample_ids->size() + 1)) / 8));
	fprintf(stderr, "Using %d bits per genotype encoding in total %d bytes.\n", n_bits_per_sample, n_bytes_per_var);

	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Assign the starting index to each EOI region.
		vector<t_annot_region*>* cur_chr_EOI_regs = get_regions_per_chromosome(EOI_regs, chr_ids->at(i_chr));
		fprintf(stderr, "%d EOI regions on %s\n", (int)cur_chr_EOI_regs->size(), chr_ids->at(i_chr));
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		// Sanity check.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// This array holds a vector for each track position.
		vector<void**>** per_track_posn_encoded_geno_data = new vector<void**>*[cur_chr_EOI_covg + 2];
		for (int track_posn = 0; track_posn < cur_chr_EOI_covg; track_posn++)
		{
			per_track_posn_encoded_geno_data[track_posn] = new vector<void**>();
		} // track_posn loop.

		// 1       9374375 rs926250        G       A       100     PASS    AC=3207;AF=0.640375;AN=5008;NS=2504;DP=13999;EAS_AF=0.2698;AMR_AF=0.6023;AFR_AF=0.8759;EUR_AF=0.6889;SAS_AF=0.681;AA=A|||;VT=SNP        GT      1|0
		char cur_chrom_vcf_fp[1000];
		sprintf(cur_chrom_vcf_fp, "%s/%s.vcf.gz", per_chrom_VCF_dir, chr_ids->at(i_chr));
		if (check_file(cur_chrom_vcf_fp))
		{
			fprintf(stderr, "Reading variants in %s\n", cur_chrom_vcf_fp);
		}
		else
		{
			fprintf(stderr, "Could not find VCF file @ %s\n", cur_chrom_vcf_fp);
		}

		FILE* f_vcf = open_f(cur_chrom_vcf_fp, "r");
		int CHROM_col_i = 0;
		int POSITION_col_i = 1;
		int VAR_ID_col_i = 2;
		int REF_ALLELE_col_i = 3;
		int ALT_ALLELE_col_i = 4;
		int GENOTYPE_col_i = 9;

		char** per_col_buffers = new char*[10];
		for (int col_i = 0; col_i < 10; col_i++)
		{
			per_col_buffers[col_i] = new char[1000];
		} // col_i loop.
		char geno_str_buff[100];

		int n_processed_dels = 0;
		int n_skipped_dels = 0;
		while (1)
		{
			char* cur_line = getline(f_vcf);
			if (cur_line == NULL)
			{
				break;
			}

			if (cur_line[0] == '#')
			{
				delete[] cur_line;
				continue;
			}

			// Read the first 9 columns.
			int i_cur_char = 0;
			for (int col_i = 0; col_i < 9; col_i++)
			{
				if (t_string::get_next_token(cur_line, per_col_buffers[col_i], 1000, "\t", i_cur_char) == false)
				{
					fprintf(stderr, "Could not read required fields: %s\n", cur_line);
					exit(0);
				}
			} // col_i loop.

			// Parse the location: Insert location is used as is, unlike deletions.
			//char* var_chr = t_string::copy_me_str(cur_line_toks->at(CHROM_col_i)->str());
			char* var_chr = t_string::copy_me_str(per_col_buffers[CHROM_col_i]);
			//int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str());
			int var_posn = atoi(per_col_buffers[POSITION_col_i]) + 1;

			//int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_ref_allele = t_string::string_length(per_col_buffers[REF_ALLELE_col_i]);
			//int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(per_col_buffers[ALT_ALLELE_col_i]);

			// Check deletion criterion.
			if (l_alt_allele != 1 ||
				l_ref_allele < 2)
			{
				delete[] cur_line;
				n_skipped_dels++;
				continue;
			}

			// We take the first allele in every SNP.
			char* ref_allele = per_col_buffers[REF_ALLELE_col_i];
			char* alt_allele = per_col_buffers[ALT_ALLELE_col_i];

			// We take the first allele in every SNP.
			int l_deletion = -1;
			if (l_ref_allele > l_alt_allele)
			{
				l_deletion = l_ref_allele - l_alt_allele;
			}
			else
			{
				delete[] cur_line;
				continue;
			}

			n_processed_dels++;

			if (n_processed_dels % 1000 == 0)
			{
				fprintf(stderr, "@%d. Deletion (%d skipped).         \r", n_processed_dels, n_skipped_dels);
			}

			if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
			{
				fprintf(stderr, "l_indel: %d; %s\n", l_deletion, cur_line);
			}

			// Search the current variant.
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, (int)cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			// For every overlap, update the allelic signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < (int)cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				int overlap_start = MAX(cur_chr_EOI_regs->at(cur_EOI_reg_i)->start, var_posn);
				int overlap_end = MIN(cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn + l_deletion - 1);

				// Make sure that this deletion is engulfed in this EOI.
				if (overlap_start <= overlap_end &&
					overlap_start <= var_posn &&
					overlap_end >= var_posn + l_deletion - 1)
				{
					if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding Deletion Overlap %s: %s:%d-%d:%s ; [%d-%d]\n", cur_line,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							overlap_start, overlap_end);
					}

					// Parse the per sample genotypes.
					int* per_sample_genotypes = new int[(int)sample_ids->size() + 2];
					for (int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
					{
						int cur_sample_genotype = 0;
						if (t_string::get_next_token(cur_line, geno_str_buff, 100, "\t", i_cur_char) == false)
						{
							fprintf(stderr, "Could not parse %d. sample's genotype: %s\n", sample_i, cur_line);
							exit(0);
						}

						char* cur_sample_geno_str = geno_str_buff;
						int geno0 = (int)(cur_sample_geno_str[0] - '0');
						int geno1 = (int)(cur_sample_geno_str[2] - '0');
						if (genotype_encoding_type == GENOTYPE_ENC_GENO)
						{
							cur_sample_genotype = geno0 + geno1;
						}
						else if (genotype_encoding_type == GENOTYPE_ENC_HAPLO)
						{
							cur_sample_genotype = (geno0 << 1) + geno1;
						}

						// Check the assigned genotype value.
						if (genotype_encoding_type == GENOTYPE_ENC_GENO && 
							cur_sample_genotype > 2)
						{
							fprintf(stderr, "Sanity check failed; genotype is greater than 2.\n");
							exit(0);
						}

						if (genotype_encoding_type == GENOTYPE_ENC_HAPLO &&
							cur_sample_genotype > 3)
						{
							fprintf(stderr, "Sanity check failed; genotype is greater than 3.\n");
							exit(0);
						}

						per_sample_genotypes[sample_i] = cur_sample_genotype;

						if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
						{
							fprintf(stderr, "Sample_i: %d ; Genotype=%s (%d)\n", sample_i,
								cur_sample_geno_str, cur_sample_genotype);
						}
					} // sample_i loop.

					// Set all the nucleotides between start and end to 1.
					for (int del_nuc_posn = overlap_start;
						del_nuc_posn <= overlap_end;
						del_nuc_posn++)
					{
						// For each overlap, update ref and alternate alleles.
						int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (del_nuc_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);						

						unsigned char* cur_var_per_sample_encoded_genotypes = new unsigned char[n_bytes_per_var + 2];
						memset(cur_var_per_sample_encoded_genotypes, 0, sizeof(unsigned char) * n_bytes_per_var);

						if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)

						{
							fprintf(stderr, "Adding encoded genotypes (%d) @ %d(%d)\n", genotype_encoding_type, del_nuc_posn, track_posn);
						}

						for (int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
						{
							int byte_i = (sample_i * n_bits_per_sample) / 8;
							int bit_i = (sample_i * n_bits_per_sample) % 8;
							int cur_sample_genotype = per_sample_genotypes[sample_i];

							if (genotype_encoding_type == GENOTYPE_ENC_GENO &&
								cur_sample_genotype > 2)
							{
								fprintf(stderr, "Sanity check failed; genotype is greater than 2.\n");
								exit(0);
							}		

							if (genotype_encoding_type == GENOTYPE_ENC_HAPLO &&
								cur_sample_genotype > 3)
							{
								fprintf(stderr, "Sanity check failed; genotype is greater than 3.\n");
								exit(0);
							}

							if (byte_i >= n_bytes_per_var)
							{
								fprintf(stderr, "Sanity check failed byte_i out of range: %d > %d\n", byte_i, n_bytes_per_var);
								exit(0);
							}

							if (track_posn >= cur_chr_EOI_covg)
							{
								fprintf(stderr, "Track position is out of range: %d > %d\n", track_posn, cur_chr_EOI_covg);
								exit(0);
							}

							if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
							{
								fprintf(stderr, "Sample_i: %d ; Genotype=%d; byte_i=%d/%d, bit_i=%d\n", sample_i,
									cur_sample_genotype,
									byte_i, n_bytes_per_var, bit_i);
							}

							unsigned char cur_geno_added_byte = -1;
							if (n_bits_per_sample == 1)
							{
								unsigned char geno_existence = (0) ? (cur_sample_genotype == 0) : (1);
								cur_geno_added_byte = (geno_existence << bit_i);

								if (((1 << bit_i) & cur_var_per_sample_encoded_genotypes[byte_i]) > 0)
								{
									fprintf(stderr, "Sanity check failed: This value is already set for this individual.\n");
									exit(0);
								}

							}
							else if (n_bits_per_sample == 2)
							{
								cur_geno_added_byte = (cur_sample_genotype << bit_i);

								if(((3 << bit_i) & cur_var_per_sample_encoded_genotypes[byte_i]) > 0)
								{
									fprintf(stderr, "Sanity check failed: This value is already set for this individual.\n");
									exit(0);
								}
							}

							if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
							{
								fprintf(stderr, "%d: Sample %d (%d): %X; OR'ing with: %X\n",
									var_posn,
									sample_i, per_sample_genotypes[sample_i],
									cur_var_per_sample_encoded_genotypes[byte_i],
									cur_geno_added_byte);
							}

							cur_var_per_sample_encoded_genotypes[byte_i] = cur_var_per_sample_encoded_genotypes[byte_i] | cur_geno_added_byte;
						} // col_i loop.						

						// Add this variant.
						void** cur_var_encode_geno_data = new void*[2];
						int* cur_track_posn_allele = new int[2];
						cur_track_posn_allele[0] = 0;
						cur_track_posn_allele[1] = 0; // Allele is meaningless in this context.
						cur_var_encode_geno_data[0] = cur_track_posn_allele;
						cur_var_encode_geno_data[1] = cur_var_per_sample_encoded_genotypes;

						// Add this genotype to the current track position.
						per_track_posn_encoded_geno_data[track_posn]->push_back(cur_var_encode_geno_data);
					} // del_nuc_posn loop.

					delete[] per_sample_genotypes;
				} // overlap check.

				cur_EOI_reg_i++;
			} // looping over EOIs.

			delete[] cur_line;
		} // VCF file reading loop.

		// Write and clean memory.
		fprintf(stderr, "Saving encoded genotypes matrix.\n");
		char matrix_op_fp[1000];
		sprintf(matrix_op_fp, "%s/%s_matrix.bin.gz", genotypes_op_dir, chr_ids->at(i_chr));
		FILE* f_matrix_op = open_f(matrix_op_fp, "w");

		unsigned char* empty_geno_array = new unsigned char[n_bytes_per_var + 2];
		memset(empty_geno_array, 0, sizeof(unsigned char) * n_bytes_per_var);

		for (int track_posn = 0; track_posn < cur_chr_EOI_covg; track_posn++)
		{
			if (per_track_posn_encoded_geno_data[track_posn]->size() == 0)
			{
				// >>>>>>>>>Note that these must be included in the secure analysis case.<<<<<<<<<
				//// Dont write empty genotype arrays.
				//fwrite(&track_posn, sizeof(int), 1, f_matrix_op);
				//fwrite(empty_geno_array, sizeof(unsigned char), n_bytes_per_var, f_matrix_op);
			}
			else
			{
				// Write the position.
				fwrite(&track_posn, sizeof(int), 1, f_matrix_op);

				// First, build the OR'ed array of deletion states.
				unsigned char* cur_posn_genotypes = new unsigned char[n_bytes_per_var];
				memset(cur_posn_genotypes, 0, sizeof(unsigned char) * n_bytes_per_var);

				fprintf(stderr, "Pooling %d variants at track position %d            \r", (int)per_track_posn_encoded_geno_data[track_posn]->size(), track_posn);

				for (int var_i = 0; var_i < (int)per_track_posn_encoded_geno_data[track_posn]->size(); var_i++)
				{
					void** cur_genotype_data = per_track_posn_encoded_geno_data[track_posn]->at(var_i);
					unsigned char* cur_var_per_sample_encoded_genotypes = (unsigned char*)(cur_genotype_data[1]);

					for (int byte_i = 0; byte_i < n_bytes_per_var; byte_i++)
					{
						cur_posn_genotypes[byte_i] = cur_posn_genotypes[byte_i] | cur_var_per_sample_encoded_genotypes[byte_i];
					} // byte_i loop.

					delete[] cur_var_per_sample_encoded_genotypes;
				} // var_i loop

				// Now write the aggregated genotypes.
				fwrite(cur_posn_genotypes, sizeof(unsigned char), n_bytes_per_var, f_matrix_op);

				delete[] cur_posn_genotypes;
			}
		} // track posn loop.

		close_f(f_matrix_op, matrix_op_fp);
	} // chr ids loop.
}

// Genotype processing: This should be done SNV/Insert/Del specific manner.
void encode_SNV_genotypes_2_matrix_per_VCF(char* per_chrom_VCF_dir, char* sample_ids_list_fp, 
	char* EOI_regs_BED_fp,
	int genotype_encoding_type, char* genotypes_op_dir)
{
	char* genotype_encoding_strs[3];
	genotype_encoding_strs[GENOTYPE_ENC_EXISTENCE] = t_string::copy_me_str("GENO_ENC_EXISTENCE");
	genotype_encoding_strs[GENOTYPE_ENC_GENO] = t_string::copy_me_str("GENO_ENC_GENO");
	genotype_encoding_strs[GENOTYPE_ENC_HAPLO] = t_string::copy_me_str("GENO_ENC_HAPLO");

	fprintf(stderr, "Encoding SNV genotypes under %s using genotype encoding %s\n", per_chrom_VCF_dir, genotype_encoding_strs[genotype_encoding_type]);

	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Encoding SNV genotype matrix on %d regions.\n", (int)EOI_regs->size());

	vector<char*>* chr_ids = get_chr_ids(EOI_regs);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Using %d samples in VCF data under %s\n", (int)sample_ids->size(), per_chrom_VCF_dir);

	int n_bits_per_sample = 2;
	if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
	{
		n_bits_per_sample = 1;
	}

	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * ((int)sample_ids->size() + 1)) / 8));
	fprintf(stderr, "Using %d bits per genotype encoding in total %d bytes.\n", n_bits_per_sample, n_bytes_per_var);

	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Assign the starting index to each EOI region.
		vector<t_annot_region*>* cur_chr_EOI_regs = get_regions_per_chromosome(EOI_regs, chr_ids->at(i_chr));
		fprintf(stderr, "%d EOI regions on %s\n", (int)cur_chr_EOI_regs->size(), chr_ids->at(i_chr));
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		  // Sanity check.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}
	
		// 1       9374375 rs926250        G       A       100     PASS    AC=3207;AF=0.640375;AN=5008;NS=2504;DP=13999;EAS_AF=0.2698;AMR_AF=0.6023;AFR_AF=0.8759;EUR_AF=0.6889;SAS_AF=0.681;AA=A|||;VT=SNP        GT      1|0
		char cur_chrom_vcf_fp[1000];
		sprintf(cur_chrom_vcf_fp, "%s/%s.vcf.gz", per_chrom_VCF_dir, chr_ids->at(i_chr));
		if (check_file(cur_chrom_vcf_fp))
		{
			fprintf(stderr, "Reading variants in %s\n", cur_chrom_vcf_fp);
		}
		else
		{
			fprintf(stderr, "Could not find VCF file @ %s\n", cur_chrom_vcf_fp);
		}

		// This array holds a vector for each track position.
		vector<void**>** per_track_posn_encoded_geno_data = new vector<void**>*[cur_chr_EOI_covg + 2];
		for (int track_posn = 0; track_posn < cur_chr_EOI_covg; track_posn++)
		{
			per_track_posn_encoded_geno_data[track_posn] = new vector<void**>();
		} // track_posn loop.

		FILE* f_vcf = open_f(cur_chrom_vcf_fp, "r");
		int CHROM_col_i = 0;
		int POSITION_col_i = 1;
		int VAR_ID_col_i = 2;
		int REF_ALLELE_col_i = 3;
		int ALT_ALLELE_col_i = 4;
		int GENOTYPE_col_i = 9;

		char** per_col_buffers = new char*[10];
		for (int col_i = 0; col_i < 10; col_i++)
		{
			per_col_buffers[col_i] = new char[1000];
		} // col_i loop.
		char geno_str_buff[100];

		fprintf(stderr, "Allocating the allelic signal.\n");

		int n_processed_vars = 0;

		while (1)
		{
			char* cur_line = getline(f_vcf);
			if (cur_line == NULL)
			{
				break;
			}

			if (cur_line[0] == '#')
			{
				delete[] cur_line;
				continue;
			}

			// Read the first 9 columns.
			int i_cur_char = 0;
			for (int col_i = 0; col_i < 9; col_i++)
			{
				if (t_string::get_next_token(cur_line, per_col_buffers[col_i], 1000, "\t", i_cur_char) == false)
				{
					fprintf(stderr, "Could not read required fields: %s\n", cur_line);
					exit(0);
				}
			} // col_i loop.

			// Parse the location: Insert location is used as is, unlike deletions.
			//char* var_chr = t_string::copy_me_str(cur_line_toks->at(CHROM_col_i)->str());
			char* var_chr = t_string::copy_me_str(per_col_buffers[CHROM_col_i]);
			//int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str());
			int var_posn = atoi(per_col_buffers[POSITION_col_i]);

			// 
			//int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_ref_allele = t_string::string_length(per_col_buffers[REF_ALLELE_col_i]);
			//int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(per_col_buffers[ALT_ALLELE_col_i]);

			if (l_alt_allele != 1 ||
				l_ref_allele != 1)
			{
				//t_string::clean_tokens(cur_line_toks);
				delete[] cur_line;
				continue;
			}

			n_processed_vars++;

			if (n_processed_vars % 10000 == 0)
			{
				fprintf(stderr, "@ %d. SNV.     \n", n_processed_vars);
			}

			// We take the first allele in every SNP.
			char ref_allele = per_col_buffers[REF_ALLELE_col_i][0];
			char alt_allele = per_col_buffers[ALT_ALLELE_col_i][0];

			char ref_alt_allele[10];
			ref_alt_allele[0] = ref_allele;
			ref_alt_allele[1] = alt_allele;

			// Search the current variant.
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, (int)cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			// For every overlap, update the allelic signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < (int)cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				bool overlap = (cur_chr_EOI_regs->at(cur_EOI_reg_i)->start <= var_posn &&
					cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn);

				if (overlap)
				{
					// For each overlap, update ref and alternate alleles.
					int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

					unsigned char* cur_var_per_sample_encoded_genotypes = new unsigned char[n_bytes_per_var + 2];
					memset(cur_var_per_sample_encoded_genotypes, 0, sizeof(unsigned char) * n_bytes_per_var);

					if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding encoded genotypes (%d) @ %d(%d)\n", genotype_encoding_type, var_posn, track_posn);
					}

					char* cur_ref_allele_str = per_col_buffers[REF_ALLELE_col_i];
					char* cur_alt_allele_str = per_col_buffers[ALT_ALLELE_col_i];

					for(int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
					{
						int byte_i = (sample_i * n_bits_per_sample) / 8;
						int bit_i = (sample_i * n_bits_per_sample) % 8;
						int cur_sample_genotype = 0;
						if (t_string::get_next_token(cur_line, geno_str_buff, 100, "\t", i_cur_char) == false)
						{
							fprintf(stderr, "Could not parse %d. sample's genotype: %s\n", sample_i, cur_line);
							exit(0);
						}

						char* cur_sample_geno_str = geno_str_buff;
						int geno0 = cur_sample_geno_str[0] - '0';
						int geno1 = cur_sample_geno_str[2] - '0';
						if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
						{
							if ((geno0 + geno1) > 0)
							{
								cur_sample_genotype = 1;
							}
						}
						else if (genotype_encoding_type == GENOTYPE_ENC_GENO)
						{
							cur_sample_genotype = geno0 + geno1;
						}
						else if (genotype_encoding_type == GENOTYPE_ENC_HAPLO)
						{
							cur_sample_genotype = (geno0 << 1) + geno1;
						}

						// Check the assigned genotype value.
						if (genotype_encoding_type == GENOTYPE_ENC_GENO &&
							cur_sample_genotype > 2)
						{
							fprintf(stderr, "Sanity check failed; genotype is greater than 2.\n");
							exit(0);
						}

						if (genotype_encoding_type == GENOTYPE_ENC_HAPLO &&
							cur_sample_genotype > 3)
						{
							fprintf(stderr, "Sanity check failed; genotype is greater than 2.\n");
							exit(0);
						}

						if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
						{
							fprintf(stderr, "Sample_i: %d Genotype=%s (%d); byte_i=%d/%d, bit_i=%d\n", sample_i,
								cur_sample_geno_str, cur_sample_genotype,
								byte_i, n_bytes_per_var, bit_i);
						}

						if (byte_i >= n_bytes_per_var)
						{
							fprintf(stderr, "Sanity check failed: byte_i out of range: %d > %d\n", byte_i, n_bytes_per_var);
							exit(0);
						}

						if (track_posn >= cur_chr_EOI_covg)
						{
							fprintf(stderr, "Track position is out of range: %d > %d\n", track_posn, cur_chr_EOI_covg);
							exit(0);
						}

						unsigned char cur_geno_added_byte = -1;
						if (n_bits_per_sample == 1)
						{
							unsigned char geno_existence = (cur_sample_genotype == 0) ? (0) : (1);
							cur_geno_added_byte = (geno_existence << bit_i);

							if (((1 << bit_i) & cur_var_per_sample_encoded_genotypes[byte_i]) > 0)
							{
								fprintf(stderr, "Sanity check failed: This value is already set for this individual.\n");
								exit(0);
							}
						}
						else if (n_bits_per_sample == 2)
						{
							cur_geno_added_byte = (cur_sample_genotype << bit_i);

							if (((3 << bit_i) & cur_var_per_sample_encoded_genotypes[byte_i]) > 0)
							{
								fprintf(stderr, "Sanity check failed: This value is already set for this individual.\n");
								exit(0);
							}
						}

						if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
						{
							fprintf(stderr, "%d: Sample %d (%s): %X; OR'ing with: %X\n", 
									var_posn,
									sample_i, cur_sample_geno_str, 
									cur_var_per_sample_encoded_genotypes[byte_i],
									cur_geno_added_byte);
						}

						cur_var_per_sample_encoded_genotypes[byte_i] = cur_var_per_sample_encoded_genotypes[byte_i] | cur_geno_added_byte;
					} // col_i loop.

					// Add this variant.
					void** cur_var_encode_geno_data = new void*[2];
					int* cur_track_posn_allele = new int[2];
					cur_track_posn_allele[0] = track_posn;
					cur_track_posn_allele[1] = nuc_2_num(alt_allele);
					cur_var_encode_geno_data[0] = cur_track_posn_allele;
					cur_var_encode_geno_data[1] = cur_var_per_sample_encoded_genotypes;

					// Add this genotype to the current track position.
					per_track_posn_encoded_geno_data[track_posn]->push_back(cur_var_encode_geno_data);
				} // overlap check.

				cur_EOI_reg_i++;
			} // looping over EOIs.

			delete[] cur_line;
		} // VCF file reading loop.

		// Write and clean memory.
		fprintf(stderr, "Saving encoded genotypes matrix.\n");
		char matrix_op_fp[1000];
		sprintf(matrix_op_fp, "%s/%s_matrix.bin.gz", genotypes_op_dir, chr_ids->at(i_chr));
		FILE* f_matrix_op = open_f(matrix_op_fp, "w");

		unsigned char* empty_geno_array = new unsigned char[n_bytes_per_var + 2];
		memset(empty_geno_array, 0, sizeof(unsigned char) * n_bytes_per_var);

		for (int track_posn = 0; track_posn < cur_chr_EOI_covg; track_posn++)
		{
			if (per_track_posn_encoded_geno_data[track_posn]->size() == 0)
			{
				// >>>>>>>>>Note that these must be included in the secure analysis case.<<<<<<<<<
				// >>>>>>>>>Note that these must be included in the secure analysis case.<<<<<<<<<
				// >>>>>>>>>Note that these must be included in the secure analysis case.<<<<<<<<<
				//fwrite(&track_posn, sizeof(int), 1, f_matrix_op);

				//// Dont write empty genotype arrays.
				//for (int allele_i = 0; allele_i < 5; allele_i++)
				//{
				//	fwrite(empty_geno_array, sizeof(unsigned char), n_bytes_per_var, f_matrix_op);
				//} // allele_i loop.
			}
			else
			{
				// Write the position.
				fwrite(&track_posn, sizeof(int), 1, f_matrix_op);

				// First, build the OR'ed array of deletion states.
				for (int allele_i = 0; allele_i < 5; allele_i++)
				{
					unsigned char* cur_posn_allele_genotypes = new unsigned char[n_bytes_per_var];
					memset(cur_posn_allele_genotypes, 0, sizeof(unsigned char) * n_bytes_per_var);

					if (track_posn % 100000 == 0 ||
						per_track_posn_encoded_geno_data[track_posn]->size() > 0)
					{
						fprintf(stderr, "Pooling %d variants at track position %d            \r", per_track_posn_encoded_geno_data[track_posn]->size(), track_posn);
					}

					if (per_track_posn_encoded_geno_data[track_posn]->size() > 1)
					{
						fprintf(stderr, "Multiple variants @ track_posn: %d            \n", track_posn);
					}

					for (int var_i = 0; var_i < (int)per_track_posn_encoded_geno_data[track_posn]->size(); var_i++)
					{
						// Must check the allele.
						void** cur_genotype_data = per_track_posn_encoded_geno_data[track_posn]->at(var_i);
						int* cur_var_track_posn_allele_i = (int*)(cur_genotype_data[0]);						

						// Make sure that this allele matches.
						if (allele_i == cur_var_track_posn_allele_i[1])
						{
							unsigned char* cur_var_per_sample_encoded_genotypes = (unsigned char*)(cur_genotype_data[1]);

							for (int byte_i = 0; byte_i < n_bytes_per_var; byte_i++)
							{
								cur_posn_allele_genotypes[byte_i] = cur_posn_allele_genotypes[byte_i] | cur_var_per_sample_encoded_genotypes[byte_i];
							} // byte_i loop.

							delete[] cur_var_per_sample_encoded_genotypes;
						} // allele_i check.
					} // var_i loop

					// Now write the aggregated genotypes.
					fwrite(cur_posn_allele_genotypes, sizeof(unsigned char), n_bytes_per_var, f_matrix_op);

					delete[] cur_posn_allele_genotypes;
				} // allele_i loop.
			} // data check.
		} // track posn loop.

		close_f(f_matrix_op, matrix_op_fp);
	} // chr ids loop.
	fprintf(stderr, "\nDone Encoding SNVs.\n");
}

bool sort_encoded_geno_per_track_posn(void** encoded_geno_entry1, void** encoded_geno_entry2)
{
	int* geno1_posn_allele = (int*)(encoded_geno_entry1[0]);
	int* geno2_posn_allele = (int*)(encoded_geno_entry2[0]);

	return(geno1_posn_allele[0] < geno2_posn_allele[0]);
}

// This count and returns results for testing.
unsigned int* plain_aggregate_encoded_Del_genotype_matrix(char* encoded_vectorized_del_dir, 
	char* EOI_regs_BED_fp, 
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type)
{
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", (int)EOI_regs->size());

	if (genotype_encoding_type == GENOTYPE_ENC_HAPLO)
	{
		fprintf(stderr, "***Haplotype encoding does not work, yet.***\n");
		exit(0);
	}

	int n_bits_per_sample = 2;
	if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
	{
		n_bits_per_sample = 1;
	}

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * ((int)sample_ids->size() + 1)) / 8));
	fprintf(stderr, "Using %d bits per genotype encoding in total %d bytes.\n", n_bits_per_sample, n_bytes_per_var);

	t_restr_annot_region_list* restr_EOI_regs = restructure_annot_regions(EOI_regs);

	vector<t_annot_region*>* VOI_regs = load_BED(VOI_regs_BED_fp);
	fprintf(stderr, "Aggregating genotype signal for %d variants.\n", (int)VOI_regs->size());

	vector<char*>* chr_ids = restr_EOI_regs->chr_ids;
	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Assign the starting index to each EOI region.
		vector<t_annot_region*>* cur_chr_EOI_regs = restr_EOI_regs->regions_per_chrom[i_chr];
		fprintf(stderr, "%d EOI regions on %s\n", (int)cur_chr_EOI_regs->size(), chr_ids->at(i_chr));
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		vector<t_annot_region*>* cur_chr_VOI_regs = get_regions_per_chromosome(VOI_regs, chr_ids->at(i_chr));

		// Setup the target regions.
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		int* per_track_posn_genomic_posn = new int[cur_chr_EOI_covg + 1];
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;

			for (int cur_track_posn = cumul_signal_index;
				cur_track_posn < cumul_signal_index + cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start;
				cur_track_posn++)
			{
				per_track_posn_genomic_posn[cur_track_posn] = cur_chr_EOI_regs->at(i_reg)->start + (cur_track_posn - cumul_signal_index);
			} // cur_track_posn loop.

			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		// Sanity check on the coverage.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		char cur_chr_encoded_geno_matrix_fp[1000];
		sprintf(cur_chr_encoded_geno_matrix_fp, "%s/%s_matrix.bin.gz", encoded_vectorized_del_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_encoded_geno_matrix_fp))
		{
			fprintf(stderr, "Could not find %s\n", cur_chr_encoded_geno_matrix_fp);
			exit(0);
		}
		else
		{
			fprintf(stderr, "Loading encoded genotypes from %s\n", cur_chr_encoded_geno_matrix_fp);
		}

		vector<void**>* encoded_posn_genotype_data = new vector<void**>();
		FILE* f_geno_matrix = open_f(cur_chr_encoded_geno_matrix_fp, "rb");
		while(1)
		{
			int cur_track_posn = 0;
			if (fread(&cur_track_posn, sizeof(int), 1, f_geno_matrix) != 1)
			{
				fprintf(stderr, "Finished reading file.\n");
				break;
			}			

			unsigned char* cur_posn_del_state = new unsigned char[n_bytes_per_var+1];
			if (fread(cur_posn_del_state, sizeof(unsigned char), n_bytes_per_var, f_geno_matrix) != n_bytes_per_var)
			{
				fprintf(stderr, "Could not read %d encoded genotypes.\n", n_bytes_per_var);
				exit(0);
			}

			if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
			{
				fprintf(stderr, "Read the encoded genotypes @ %d\n", cur_track_posn);
			}

			void** cur_var_encoded_genotype_data = new void*[2];
			int* cur_posn_data = new int[2];
			cur_posn_data[0] = cur_track_posn;
			cur_var_encoded_genotype_data[0] = cur_posn_data;
			cur_var_encoded_genotype_data[1] = cur_posn_del_state;
			encoded_posn_genotype_data->push_back(cur_var_encoded_genotype_data);
		} // file reading loop.
		close_f(f_geno_matrix, cur_chr_encoded_geno_matrix_fp);

		fprintf(stderr, "Loaded encoded genotypes for %d variants.\n", (int)encoded_posn_genotype_data->size());
		for (int i_voi = 0; i_voi < (int)cur_chr_VOI_regs->size(); i_voi++)
		{
			fprintf(stderr, "--------------------------\n%s:%d\n", cur_chr_VOI_regs->at(i_voi)->chrom, cur_chr_VOI_regs->at(i_voi)->start);
			int cur_var_count = 0;

			// Build the encoded genotype matrix for the current variant.
			unsigned char** per_posn_encoded_genotypes = new unsigned char*[cur_chr_VOI_regs->at(i_voi)->end - cur_chr_VOI_regs->at(i_voi)->start + 3];
			memset(per_posn_encoded_genotypes, 0, sizeof(unsigned char*) * (cur_chr_VOI_regs->at(i_voi)->end - cur_chr_VOI_regs->at(i_voi)->start + 3));

			// Following copies the encoded genotypes arrays to a list for the positions on the current deletion variant.
			int n_copied_posns = 0;
			for (int track_pos_i = 0; track_pos_i < (int)encoded_posn_genotype_data->size(); track_pos_i++)
			{
				void** cur_pos_data = encoded_posn_genotype_data->at(track_pos_i);
				int* pos_ptr = (int*)(cur_pos_data[0]);
				unsigned char* encoded_genotypes = (unsigned char*)(cur_pos_data[1]);
				int cur_track_posn = pos_ptr[0];

				// Make sure we have the one position to the left and right of the variant.
				if (per_track_posn_genomic_posn[cur_track_posn] >= (cur_chr_VOI_regs->at(i_voi)->start - 1) &&
					per_track_posn_genomic_posn[cur_track_posn] <= (cur_chr_VOI_regs->at(i_voi)->end + 1))
				{
					fprintf(stderr, "Matched track posn. %d :: %d (%d-%d)\n", cur_track_posn, per_track_posn_genomic_posn[cur_track_posn], 
						cur_chr_VOI_regs->at(i_voi)->start-1, cur_chr_VOI_regs->at(i_voi)->end+1);

					per_posn_encoded_genotypes[per_track_posn_genomic_posn[cur_track_posn] - (cur_chr_VOI_regs->at(i_voi)->start - 1)] = encoded_genotypes;
					n_copied_posns++;
				}
			} // track_pos_i loop.

			if (n_copied_posns == 0)
			{
				continue;
			}

			// Copy the empty genotypes for the empty ones.
			unsigned char* empty_geno = new unsigned char[n_bytes_per_var + 2];
			memset(empty_geno, 0, n_bytes_per_var);

			for (int genomic_pos = (cur_chr_VOI_regs->at(i_voi)->start - 1);
				genomic_pos <= (cur_chr_VOI_regs->at(i_voi)->end + 1);
				genomic_pos++)
			{
				unsigned char* cur_var_genotypes = per_posn_encoded_genotypes[genomic_pos - (cur_chr_VOI_regs->at(i_voi)->start - 1)];
				if (cur_var_genotypes == NULL)
				{
					per_posn_encoded_genotypes[genomic_pos - (cur_chr_VOI_regs->at(i_voi)->start - 1)] = empty_geno;
				}
			} // genomic_pos loop.

			fprintf(stderr, "%d-%d: Copied %d/%d posns\n", cur_chr_VOI_regs->at(i_voi)->start, cur_chr_VOI_regs->at(i_voi)->end, n_copied_posns, cur_chr_VOI_regs->at(i_voi)->end + 1 - cur_chr_VOI_regs->at(i_voi)->start + 2);

			// Go over all the samples and count the consistent ones.
			for (int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
			{
				int cur_sample_del_geno = -1;
				int cur_sample_left_geno = -1;
				int cur_sample_right_geno = -1;

				for(int genomic_pos = cur_chr_VOI_regs->at(i_voi)->start - 1; 
					genomic_pos <= cur_chr_VOI_regs->at(i_voi)->end + 1;
					genomic_pos++)
				{
					unsigned char* cur_var_genotypes = per_posn_encoded_genotypes[genomic_pos - (cur_chr_VOI_regs->at(i_voi)->start - 1)];

					int byte_i = (sample_i * n_bits_per_sample) / 8;
					int bit_i = (sample_i * n_bits_per_sample) % 8;

					int cur_sample_geno = 0;
					if (n_bits_per_sample == 1)
					{
						cur_sample_geno = (cur_var_genotypes[byte_i] & (1 << bit_i)) >> bit_i;
					}

					if (n_bits_per_sample == 2)
					{
						cur_sample_geno = (cur_var_genotypes[byte_i] & (3 << bit_i)) >> bit_i;
					}

					if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
					{
						if (cur_sample_geno > 2)
						{
							fprintf(stderr, "The allele count is not as expected while aggregating.\n");
							exit(0);
						}
					}

					// Based on position, set the left and right genotypes.
					if (genomic_pos == (cur_chr_VOI_regs->at(i_voi)->start - 1))
					{
						cur_sample_left_geno = cur_sample_geno;
					}
					else if (genomic_pos == (cur_chr_VOI_regs->at(i_voi)->end + 1))
					{
						cur_sample_right_geno = cur_sample_geno;
					}
					else
					{
						if (cur_sample_geno == 0)
						{
							cur_sample_del_geno = 0;
							break;
						}

						if (cur_sample_del_geno == -1)
						{
							cur_sample_del_geno = cur_sample_geno;
						}

						// If the genotype state changed, set the genotype to 0.
						if (cur_sample_del_geno > 0 &&
							cur_sample_del_geno != cur_sample_geno)
						{
							fprintf(stderr, "The genotype state changed for sample %d @ %d: %d;%d\n", sample_i, genomic_pos, cur_sample_del_geno, cur_sample_geno);
							cur_sample_del_geno = 0;
							break;
						}
					}
				} // genomic_pos loop.

				// Based on the deletion states on the deletion and at the edges, update the allele count.
				if (cur_sample_del_geno > 0 &&
					cur_sample_del_geno != cur_sample_left_geno &&
					cur_sample_del_geno != cur_sample_right_geno)
				{
					if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Sample %d: %d\n", sample_i, cur_sample_del_geno);
					}

					int allele_count = cur_sample_del_geno;

					cur_var_count += allele_count;
				}
			} // sample_i loop.

			fprintf(stderr, "%s:%d (%s) -- %d/%d over positions %d/%d variant matches\n", cur_chr_VOI_regs->at(i_voi)->chrom, cur_chr_VOI_regs->at(i_voi)->start, cur_chr_VOI_regs->at(i_voi)->name, cur_var_count, (int)sample_ids->size(),
					n_copied_posns, cur_chr_VOI_regs->at(i_voi)->end + 1 - cur_chr_VOI_regs->at(i_voi)->start + 2);
		} // i_voi loop.
	} // i_chr loop.

	return(NULL);
}

// This count and returns results for testing.
unsigned int* plain_aggregate_encoded_SNV_genotype_matrix(char* encoded_vectorized_snvs_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type)
{
	fprintf(stderr, "Aggregating SNV's.\n");

	// EOI stands for elements of interest. These are the regions
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", (int)EOI_regs->size());

	// We will process each chromosome separately. Get the chromosome id's.
	vector<char*>* chr_ids = get_chr_ids(EOI_regs);

	//t_restr_annot_region_list* restr_EOI_regs = restructure_annot_regions(EOI_regs);

	vector<t_annot_region*>* VOI_regs = load_BED(VOI_regs_BED_fp);
	fprintf(stderr, "Aggregating genotype signal for %d variants.\n", VOI_regs->size());

	// Check the encoding of the genotypes.
	if (genotype_encoding_type == GENOTYPE_ENC_HAPLO)
	{
		fprintf(stderr, "***Haplotype encoding does not work, yet.***\n");
		exit(0);
	}

	int n_bits_per_sample = 2;
	if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
	{
		n_bits_per_sample = 1;
	}

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	// This is the # of bytes for each variant that necessary to store the genotypes of all individuals; We need this to allocate the genotype vectors later.
	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * ((int)sample_ids->size() + 1)) / 8));
	fprintf(stderr, "Using %d bits per genotype encoding in total %d bytes.\n", n_bits_per_sample, n_bytes_per_var);
	
	// Process each chromosome separately.
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Assign the starting index to each EOI region.
		vector<t_annot_region*>* cur_chr_EOI_regs = get_regions_per_chromosome(EOI_regs, chr_ids->at(i_chr));
		fprintf(stderr, "%d EOI regions on %s\n", (int)cur_chr_EOI_regs->size(), chr_ids->at(i_chr));

		// Following is a special sorting that is used throughout to set the "vectorized coordinate system". It first sorts with respect to start position, then wrt name of the element.
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		// We need to look at the variants of interest on this chromosome.
		vector<t_annot_region*>* cur_chr_VOI_regs = get_regions_per_chromosome(VOI_regs, chr_ids->at(i_chr));

		// Setup the target regions on this chromosome: These are the 
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		int* per_track_posn_genomic_posn = new int[cur_chr_EOI_covg + 1];
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;

			for (int cur_track_posn = cumul_signal_index;
				cur_track_posn < cumul_signal_index + cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start;
				cur_track_posn++)
			{
				per_track_posn_genomic_posn[cur_track_posn] = cur_chr_EOI_regs->at(i_reg)->start + (cur_track_posn - cumul_signal_index);
			} // cur_track_posn loop.

			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		// Sanity check on the coverage.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// Set the encoded genotype matrix name for this chromosome.
		char cur_chr_encoded_geno_matrix_fp[1000];
		sprintf(cur_chr_encoded_geno_matrix_fp, "%s/%s_matrix.bin.gz", encoded_vectorized_snvs_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_encoded_geno_matrix_fp))
		{
			fprintf(stderr, "Could not find %s\n", cur_chr_encoded_geno_matrix_fp);
			exit(0);
		}
		else
		{
			fprintf(stderr, "Loading encoded genotypes from %s\n", cur_chr_encoded_geno_matrix_fp);
		}

		// Following loads the plaintext genotype matrix in a vector. 
		// Each entry in the vector corresponds to a position on the vectorized coordinates and holds the encoded genotype vector for the individuals at that position.
		// Note that for each position on the we need 4 different vectors. This information is stored with the genotypes.
		vector<void**>* encoded_posn_genotype_data = new vector<void**>();
		FILE* f_geno_matrix = open_f(cur_chr_encoded_geno_matrix_fp, "rb");
		while (1)
		{
			// "Track" refers the "vector". 
			// First, read the position on the vectorized coordinates for this SNV (or SNP).
			int cur_track_posn = 0;
			if (fread(&cur_track_posn, sizeof(int), 1, f_geno_matrix) != 1)
			{
				fprintf(stderr, "Finished reading file.\n");
				break;
			}

			// We are storing 5 alleles (A, C, G, T, N).
			for (int allele_i = 0; allele_i < 5; allele_i++)
			{
				// Allocate the genotypes for the current position: n_bytes_per_var are necessary to store the genotypes for all individuals.
				unsigned char* cur_posn_snv_genotypes = new unsigned char[n_bytes_per_var + 1];
				if (fread(cur_posn_snv_genotypes, sizeof(unsigned char), n_bytes_per_var, f_geno_matrix) != n_bytes_per_var)
				{
					fprintf(stderr, "Could not read %d encoded genotypes.\n", n_bytes_per_var);
					exit(0);
				}

				if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
				{
					fprintf(stderr, "Read the encoded genotypes @ %d\n", cur_track_posn);
				}

				// This is the weird data we use to store the genotypes: One 2-integer array that contains the position on the vectorized coordinates. Second contains
				// the genotypes vector for this variant.
				void** cur_var_encoded_genotype_data = new void*[2];
				int* cur_posn_data = new int[2];
				cur_posn_data[0] = cur_track_posn;
				cur_posn_data[1] = allele_i;
				cur_var_encoded_genotype_data[0] = cur_posn_data;
				cur_var_encoded_genotype_data[1] = cur_posn_snv_genotypes;
				encoded_posn_genotype_data->push_back(cur_var_encoded_genotype_data);
			}
		} // file reading loop.
		close_f(f_geno_matrix, cur_chr_encoded_geno_matrix_fp);

		fprintf(stderr, "Loaded encoded genotypes for %d vectorized positions.\n", (int)encoded_posn_genotype_data->size());

		// Go over all the Variants-Of-Interest (VOI).
		for (int i_voi = 0; i_voi < (int)cur_chr_VOI_regs->size(); i_voi++)
		{
			fprintf(stderr, "--------------------------\n%s:%d\n", cur_chr_VOI_regs->at(i_voi)->chrom, cur_chr_VOI_regs->at(i_voi)->start);

			if (cur_chr_VOI_regs->at(i_voi)->start != cur_chr_VOI_regs->at(i_voi)->end)
			{
				fprintf(stderr, "Sanity check failed: The SNP start and end positions are different.\n");
				exit(0);
			}

			// This is the position of the SNV.
			int cur_SNV_posn = cur_chr_VOI_regs->at(i_voi)->start;

			// This is the final count we would like to have.
			int cur_var_count = 0;

			// We assume that the name of the variant contains reference and alternate nucleotides in the first 2 columns separated by space. 
			// For example: "1       10504   10505   A T AC=1        .       +"
			char ref_allele[100];
			char alt_allele[100];
			sscanf(cur_chr_VOI_regs->at(i_voi)->name, "%s %s", ref_allele, alt_allele);

			// This is the vector that contains the encoded genotypes for this variant.
			unsigned char* cur_snv_encoded_genotypes = NULL;

			// Following copies the encoded genotypes arrays to a list for the positions on the current deletion variant.
			// We go over all the track positions and select the one with 
			int n_copied_vector_posns = 0;
			for (int track_pos_i = 0; track_pos_i < (int)encoded_posn_genotype_data->size(); track_pos_i++)
			{
				void** cur_pos_data = encoded_posn_genotype_data->at(track_pos_i);
				int* pos_ptr = (int*)(cur_pos_data[0]);
				unsigned char* encoded_genotypes = (unsigned char*)(cur_pos_data[1]);
				int cur_track_posn = pos_ptr[0];
				int allele_i = pos_ptr[1];

				// Make sure we have the one position to the left and right of the variant.
				if (per_track_posn_genomic_posn[cur_track_posn] == cur_chr_VOI_regs->at(i_voi)->start &&
					allele_i == nuc_2_num(alt_allele[0]))
				{
					fprintf(stderr, "Matched track posn. %d :: %d (%d-%d)\n", cur_track_posn, per_track_posn_genomic_posn[cur_track_posn],
						cur_chr_VOI_regs->at(i_voi)->start - 1, cur_chr_VOI_regs->at(i_voi)->end + 1);

					cur_snv_encoded_genotypes = encoded_genotypes;
					n_copied_vector_posns++;
					break;
				}
			} // track_pos_i loop.

			// If we could not find the current variant, do not process any more.
			if (n_copied_vector_posns == 0)
			{
				continue;
			}

			// Copy the empty genotypes for the empty ones.
			unsigned char* empty_geno = new unsigned char[n_bytes_per_var + 2];
			memset(empty_geno, 0, n_bytes_per_var);

			if (cur_snv_encoded_genotypes == NULL)
			{
				cur_snv_encoded_genotypes = empty_geno;
			}

			fprintf(stderr, "%d-%d: Copied %d/%d posns\n", cur_chr_VOI_regs->at(i_voi)->start, cur_chr_VOI_regs->at(i_voi)->end, n_copied_vector_posns, cur_chr_VOI_regs->at(i_voi)->end - cur_chr_VOI_regs->at(i_voi)->start + 1);

			// Go over all the samples and count the consistent ones.
			for (int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
			{
				int byte_i = (sample_i * n_bits_per_sample) / 8;
				int bit_i = (sample_i * n_bits_per_sample) % 8;

				int cur_sample_snv_geno = 0;
				if (n_bits_per_sample == 1)
				{
					cur_sample_snv_geno = (cur_snv_encoded_genotypes[byte_i] & (1 << bit_i)) >> bit_i;
				}

				if (n_bits_per_sample == 2)
				{
					cur_sample_snv_geno = (cur_snv_encoded_genotypes[byte_i] & (3 << bit_i)) >> bit_i;
				}

				if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
				{
					if (cur_sample_snv_geno > 1)
					{
						fprintf(stderr, "The allele count is not as expected while aggregating.\n");
						exit(0);
					}
				}

				// Based on the deletion states on the deletion and at the edges, update the allele count.
				if (cur_sample_snv_geno > 0)
				{
					if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Sample %d: %d\n", sample_i, cur_sample_snv_geno);
					}

					int allele_count = cur_sample_snv_geno;
					cur_var_count += allele_count;
				}
			} // sample_i loop.

			fprintf(stderr, "%s:%d (%s) -- %d/%d over positions %d/%d variant matches\n", cur_chr_VOI_regs->at(i_voi)->chrom, cur_chr_VOI_regs->at(i_voi)->start, cur_chr_VOI_regs->at(i_voi)->name, cur_var_count, (int)sample_ids->size(),
						n_copied_vector_posns, cur_chr_VOI_regs->at(i_voi)->start);
		} // i_voi loop.
	} // i_chr loop.

	return(NULL);
}

// This count and returns results for testing.
unsigned int* plain_aggregate_encoded_SNV_genotype_Full_Matrix(char* encoded_vectorized_snvs_dir,
	char* EOI_regs_BED_fp,
	char* VOI_regs_BED_fp,
	char* sample_ids_list_fp,
	int genotype_encoding_type)
{
	fprintf(stderr, "Aggregating SNV's.\n");

	// EOI stands for elements of interest. These are the regions
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", (int)EOI_regs->size());

	// We will process each chromosome separately. Get the chromosome id's.
	vector<char*>* chr_ids = get_chr_ids(EOI_regs);

	// This is the Variants-of-Interest, that contains the list of SNVs that we would like to aggregate the frequency.
	vector<t_annot_region*>* VOI_regs = load_BED(VOI_regs_BED_fp);
	fprintf(stderr, "Aggregating genotype signal for %d variants.\n", (int)VOI_regs->size());

	// Check the encoding of the genotypes.
	int n_bits_per_sample = 2;
	if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
	{
		n_bits_per_sample = 1;
	}
	else if (genotype_encoding_type == GENOTYPE_ENC_GENO)
	{
		n_bits_per_sample = 2;
	}
	else
	{
		fprintf(stderr, "***Haplotype encoding does not work, yet.***\n");
		exit(0);
	}

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", sample_ids->size());

	// This is the # of bytes for each variant that necessary to store the genotypes of all individuals; We need this to allocate/load the genotype matrix later.
	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * (sample_ids->size() + 1)) / 8));
	fprintf(stderr, "Using %d bits per genotype encoding in total %d bytes.\n", n_bits_per_sample, n_bytes_per_var);

	// Process each chromosome separately.
	for (int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		// Assign the starting index to each EOI region.
		vector<t_annot_region*>* cur_chr_EOI_regs = get_regions_per_chromosome(EOI_regs, chr_ids->at(i_chr));
		fprintf(stderr, "%d EOI regions on %s\n", cur_chr_EOI_regs->size(), chr_ids->at(i_chr));

		// Following is a special sorting that is used throughout to set the "vectorized coordinate system". It first sorts with respect to start position, then wrt name of the element.
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		// We need to look at the variants of interest on this chromosome.
		vector<t_annot_region*>* cur_chr_VOI_regs = get_regions_per_chromosome(VOI_regs, chr_ids->at(i_chr));

		// Setup the target regions on this chromosome.
		// This is necessary to setup the mapping between "genomic coordinates" and the "vectorized coordinates".
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		int* per_track_posn_genomic_posn = new int[cur_chr_EOI_covg + 1];
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;

			// Following sets mapping between vector-coordinates to the genomic-coordinates. Note that I write "track", which means "vector".
			for (int cur_track_posn = cumul_signal_index;
				cur_track_posn < cumul_signal_index + cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start;
				cur_track_posn++)
			{
				per_track_posn_genomic_posn[cur_track_posn] = cur_chr_EOI_regs->at(i_reg)->start + (cur_track_posn - cumul_signal_index);
			} // cur_track_posn loop.

			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		// Sanity check on the coverage.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// Set the encoded genotype matrix name for this chromosome.
		char cur_chr_encoded_geno_matrix_fp[1000];
		sprintf(cur_chr_encoded_geno_matrix_fp, "%s/%s_matrix.bin.gz", encoded_vectorized_snvs_dir, chr_ids->at(i_chr));
		if (!check_file(cur_chr_encoded_geno_matrix_fp))
		{
			fprintf(stderr, "Could not find %s\n", cur_chr_encoded_geno_matrix_fp);
			exit(0);
		}
		else
		{
			fprintf(stderr, "Loading encoded genotypes from %s\n", cur_chr_encoded_geno_matrix_fp);
		}

		// Following allocates and then loads the plaintext genotype matrix into a 3-dimensional array.
		// First dimension holds the alternate allele (the nucleotide A,C,G,T)
		// Second dimension holds the vectorized position.
		// Third dimension holds the byte index for the encoded genotypes.
		unsigned char*** per_allele_encoded_posn_genotype_data = new unsigned char**[5];
		for (int allele = 0; allele < 5; allele++)
		{
			per_allele_encoded_posn_genotype_data[allele] = new unsigned char*[cur_chr_EOI_covg + 3];
			for (int track_posn = 0; track_posn < cur_chr_EOI_covg; track_posn++)
			{
				// Load the and initialize to 0.
				per_allele_encoded_posn_genotype_data[allele][track_posn] = new unsigned char[n_bytes_per_var];
				memset(per_allele_encoded_posn_genotype_data[allele][track_posn], 0, sizeof(unsigned char) * n_bytes_per_var);
			} // posn loop.
		} // allele loop.

		FILE* f_geno_matrix = open_f(cur_chr_encoded_geno_matrix_fp, "rb");
		while (1)
		{
			// "Track" refers the "vector". 
			// First, read the position on the vectorized coordinates for this SNV (or SNP).
			int cur_track_posn = 0;
			if (fread(&cur_track_posn, sizeof(int), 1, f_geno_matrix) != 1)
			{
				fprintf(stderr, "Finished reading file.\n");
				break;
			}

			// We are storing 5 alleles (A, C, G, T, N): Load these in turn.
			for (int allele_i = 0; allele_i < 5; allele_i++)
			{
				// Allocate the genotypes for the current position: n_bytes_per_var are necessary to store the genotypes for all individuals.
				unsigned char* cur_posn_snv_genotypes = new unsigned char[n_bytes_per_var + 1];
				if (fread(cur_posn_snv_genotypes, sizeof(unsigned char), n_bytes_per_var, f_geno_matrix) != n_bytes_per_var)
				{
					fprintf(stderr, "Could not read %d encoded genotypes.\n", n_bytes_per_var);
					exit(0);
				}

				if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
				{
					fprintf(stderr, "Read the encoded genotypes @ %d\n", cur_track_posn);
				}

				per_allele_encoded_posn_genotype_data[allele_i][cur_track_posn] = cur_posn_snv_genotypes;
			}
		} // file reading loop.
		close_f(f_geno_matrix, cur_chr_encoded_geno_matrix_fp);

		// Go over all the Variants-Of-Interest (VOI).
		for (int i_voi = 0; i_voi < (int)cur_chr_VOI_regs->size(); i_voi++)
		{
			fprintf(stderr, "--------------------------\n%s:%d\n", cur_chr_VOI_regs->at(i_voi)->chrom, cur_chr_VOI_regs->at(i_voi)->start);

			if (cur_chr_VOI_regs->at(i_voi)->start != cur_chr_VOI_regs->at(i_voi)->end)
			{
				fprintf(stderr, "Sanity check failed: The SNP start and end positions are different.\n");
				exit(0);
			}

			// This is the position of the SNV.
			int cur_SNV_posn = cur_chr_VOI_regs->at(i_voi)->start;

			// This is the final count we would like to have.
			int cur_var_count = 0;

			// We assume that the name of the variant contains reference and alternate nucleotides in the first 2 columns separated by space. 
			// For example: "1       10504   10505   A T AC=1        .       +"
			char ref_allele[100];
			char alt_allele[100];
			sscanf(cur_chr_VOI_regs->at(i_voi)->name, "%s %s", ref_allele, alt_allele);

			// Find the current SNV's vectorized coordinate from the genomic-2-vector coordinates mapping that we set above.
			unsigned char* cur_snv_encoded_genotypes = NULL;
			for (int track_pos = 0; track_pos < cur_chr_EOI_covg; track_pos++)
			{
				// Check the mapping between the SNV's genomic position and the vectorized coordinate.
				if (per_track_posn_genomic_posn[track_pos] == cur_SNV_posn)
				{
					cur_snv_encoded_genotypes = per_allele_encoded_posn_genotype_data[nuc_2_num(alt_allele[0])][track_pos];
					break;
				}
			} // track_pos loop.

			// If we could not find this SNV, skip it and move to the next SNV.
			if (cur_snv_encoded_genotypes == NULL)
			{
				continue;
			}

			fprintf(stderr, "Found genotypes for SNV @ %d, counting the alleles over %d samples\n", cur_SNV_posn, sample_ids->size());

			// Go over all the samples and count the consistent ones.
			for (int sample_i = 0; sample_i < (int)sample_ids->size(); sample_i++)
			{
				// Find the byte index that is corresponding to this sample, then select the bit index.
				int byte_i = (sample_i * n_bits_per_sample) / 8;
				int bit_i = (sample_i * n_bits_per_sample) % 8;

				// Depending on the # of bits per sample, extract the genotype.
				int cur_sample_snv_geno = 0;
				if (n_bits_per_sample == 1)
				{
					cur_sample_snv_geno = (cur_snv_encoded_genotypes[byte_i] & (1 << bit_i)) >> bit_i;
				}

				if (n_bits_per_sample == 2)
				{
					cur_sample_snv_geno = (cur_snv_encoded_genotypes[byte_i] & (3 << bit_i)) >> bit_i;
				}

				// This is a sanity check on the extracted value, can skip.
				if (genotype_encoding_type == GENOTYPE_ENC_EXISTENCE)
				{
					if (cur_sample_snv_geno > 2)
					{
						fprintf(stderr, "The allele count is not as expected while aggregating.\n");
						exit(0);
					}
				}

				if (cur_sample_snv_geno > 0)
				{
					if (__DUMP_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Sample %d: %d\n", sample_i, cur_sample_snv_geno);
					}

					int allele_count = cur_sample_snv_geno;
					cur_var_count += allele_count;
				}

			} // sample_i loop.

			fprintf(stderr, "%s:%d (%s) -- %d/%d over positions.\n", 
					cur_chr_VOI_regs->at(i_voi)->chrom, cur_SNV_posn, cur_chr_VOI_regs->at(i_voi)->name, 
					cur_var_count, sample_ids->size());
		} // i_voi loop.
	} // i_chr loop.

	return(NULL);
}