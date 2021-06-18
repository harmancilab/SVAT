#include <stdio.h>
#include <stdlib.h>
#include "svat_rng.h"
#include "svat_seed_manager.h"
#include "file_utils.h"
#include "svat_ansi_string.h"
#include "svat_annot_region_tools.h"
#include "svat_genomics_coords.h"
#include "svat_genome_sequence_tools.h"
#include "svat_nucleotide.h"
#include "svat_utils.h"
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <algorithm>
#include <vector>
using namespace std;

bool __DUMP_CRYPTANNOT_SECURE_SNV_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

void encrypt_vectorized_SNVs(char* vectorized_SNVs_dir, char* EOI_regs_BED_fp, char* encrypted_SNVs_vector_op_dir)
{

}

// Following is the plaintext version of the code. 
void secure_multiply_SNV_variant_and_annotation_signals(char* plaintext_annotation_signal_dir, char* encrypted_SNV_vector_dir, char* op_dir)
{
	// Each chromosome is processed separately. We first load the list of chromosome id's of the annotation.
	char anno_sig_chr_ids_fp[1000];
	sprintf(anno_sig_chr_ids_fp, "%s/chr_ids.txt", plaintext_annotation_signal_dir);
	vector<char*>* anno_chr_ids = buffer_file(anno_sig_chr_ids_fp);
	if (anno_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotation chromosome id's from %s\n", anno_sig_chr_ids_fp);
		exit(0);
	}

	// Load the variant vector's chromosome id's.
	char var_sig_chr_ids_fp[1000];
	sprintf(var_sig_chr_ids_fp, "%s/chr_ids.txt", encrypted_SNV_vector_dir);
	vector<char*>* var_chr_ids = buffer_file(var_sig_chr_ids_fp);
	if (var_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotation chromosome id's from %s\n", var_sig_chr_ids_fp);
		exit(0);
	}

	for (int i_chr = 0; i_chr < var_chr_ids->size(); i_chr++)
	{
		// First load the EOI regions.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", encrypted_SNV_vector_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_sorted_BED_fp))
		{
			fprintf(stderr, "Skipping %s, no impact signals for this chromosome.\n", var_chr_ids->at(i_chr));
		}

		// Load and sort the EOI regions.
		vector<t_annot_region*>* sorted_cur_chr_EOI_regs = load_BED(cur_sorted_BED_fp);
		sort_set_sorting_info(sorted_cur_chr_EOI_regs, sort_regions_per_start_end_name);

		int cur_chr_EOI_covg = coverage(sorted_cur_chr_EOI_regs);
		fprintf(stderr, "Loaded %d EOI regions with %d nucleotide coverage on %s.\n",
			sorted_cur_chr_EOI_regs->size(),
			cur_chr_EOI_covg,
			var_chr_ids->at(i_chr));

		// Check signal files.
		char cur_chr_annotation_signal_fp[1000];
		sprintf(cur_chr_annotation_signal_fp, "%s/impact_signal_%s.bin.gz", plaintext_annotation_signal_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_annotation_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_annotation_signal_fp, var_chr_ids->at(i_chr));
			exit(0);
		}

		char cur_chr_variant_signal_fp[1000];
		sprintf(cur_chr_variant_signal_fp, "%s/variant_signal_%s.bin.gz", encrypted_SNV_vector_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_variant_signal_fp, var_chr_ids->at(i_chr));
			exit(0);
		}

		// Load the annotation signals for each of the 4 alleles in a 64bit unsigned integer array. 
		fprintf(stderr, "Loading per allele annotation signals.\n");
		unsigned long long** per_allele_annotation_signal = new unsigned long long*[4];
		FILE* f_annotation_signal = open_f(cur_chr_annotation_signal_fp, "rb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Reading impacts per allele %d\n", i_allele);
			per_allele_annotation_signal[i_allele] = new unsigned long long[cur_chr_EOI_covg + 2];
			if (fread(per_allele_annotation_signal[i_allele], sizeof(unsigned long long), cur_chr_EOI_covg, f_annotation_signal) != cur_chr_EOI_covg)
			{
				fprintf(stderr, "Could not read allele_i=%d from %s\n", i_allele, cur_chr_annotation_signal_fp);
				exit(0);
			}
		} // i_allele loop.
		close_f(f_annotation_signal, cur_chr_annotation_signal_fp);

		// Load the encrypted variant vector.
		fprintf(stderr, "Loading per allele variant signals.\n");
		unsigned int** per_allele_variant_signal = new unsigned int*[4];
		FILE* f_variant_signal = open_f(cur_chr_variant_signal_fp, "rb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Reading impacts per allele %d\n", i_allele);
			per_allele_variant_signal[i_allele] = new unsigned int[cur_chr_EOI_covg + 2];
			if (fread(per_allele_variant_signal[i_allele], sizeof(unsigned int), cur_chr_EOI_covg, f_variant_signal) != cur_chr_EOI_covg)
			{
				fprintf(stderr, "Could not read allele_i=%d from %s\n", i_allele, cur_chr_variant_signal_fp);
				exit(0);
			}
		} // i_allele loop.
		close_f(f_variant_signal, cur_chr_variant_signal_fp);

		// Allocate the multiplied signal.
		fprintf(stderr, "Allocating per allele annotated variants signal.\n");
		unsigned long long** per_allele_annotation_cross_variant_signal = new unsigned long long*[4];
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			per_allele_annotation_cross_variant_signal[i_allele] = new unsigned long long[cur_chr_EOI_covg + 2];
		} // i_allele loop.

		  // Simply multiply every integer on variant and annotation signals.
		for (int i = 0; i < cur_chr_EOI_covg; i++)
		{
			for (int i_allele = 0; i_allele < 4; i_allele++)
			{
				per_allele_annotation_cross_variant_signal[i_allele][i] = per_allele_variant_signal[i_allele][i] * per_allele_annotation_signal[i_allele][i];
			} // i_allele loop.
		} // i loop.

		  // Write the crossed signal to output directory.
		fprintf(stderr, "Saving annotated variant signal.\n");
		char cur_chr_crossed_signal_fp[1000];
		sprintf(cur_chr_crossed_signal_fp, "%s/annotated_variant_signal_%s.bin.gz", op_dir, var_chr_ids->at(i_chr));
		FILE* f_annotated_variant_signal = open_f(cur_chr_crossed_signal_fp, "wb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Saving impacts per allele %d\n", i_allele);
			fwrite(per_allele_annotation_cross_variant_signal[i_allele], sizeof(unsigned long long), cur_chr_EOI_covg, f_annotated_variant_signal);
		} // i_allele loop.
		close_f(f_annotated_variant_signal, cur_chr_crossed_signal_fp);

		// Save the target regions.
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, var_chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, sorted_cur_chr_EOI_regs);
	} // i_chr loop.

	char impact_string_fp[1000];
	sprintf(impact_string_fp, "%s/impact_value_strings.list", plaintext_annotation_signal_dir);
	char impact_string_op_fp[1000];
	sprintf(impact_string_op_fp, "%s/impact_value_strings.list", op_dir);
	fprintf(stderr, "Copying the impact strings.\n");
	copy_file(impact_string_fp, impact_string_op_fp);

	fprintf(stderr, "Saving the chromosome id's.\n");
	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", op_dir);
	FILE* f_chr_op = open_f(chr_ids_fp, "w");
	for (int i_chr = 0; i_chr < var_chr_ids->size(); i_chr++)
	{
		fprintf(f_chr_op, "%s\n", var_chr_ids->at(i_chr));
	} // i_chr loop
	fclose(f_chr_op);
}
