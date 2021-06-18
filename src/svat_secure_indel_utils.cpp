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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <algorithm>
#include <vector>
using namespace std;

bool __DUMP_CRYPTANNOT_INDEL_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

void encrypt_vectorized_Indels(char* vectorized_Indels_dir, char* EOI_regs_BED_fp, char* encrypted_Indel_vector_op_dir)
{

}

// Following is the plaintext version of the code. 
void secure_multiply_Indel_variant_and_annotation_signals(char* plaintext_annotation_signal_dir, char* encrypted_Indel_vector_dir, char* op_dir)
{
	// Load the annotation signal.
	char anno_sig_chr_ids_fp[1000];
	sprintf(anno_sig_chr_ids_fp, "%s/chr_ids.txt", plaintext_annotation_signal_dir);
	vector<char*>* anno_chr_ids = buffer_file(anno_sig_chr_ids_fp);
	if (anno_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotation chromosome id's from %s\n", anno_sig_chr_ids_fp);
		exit(0);
	}

	// Load the variant signal.
	char var_sig_chr_ids_fp[1000];
	sprintf(var_sig_chr_ids_fp, "%s/chr_ids.txt", encrypted_Indel_vector_dir);
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
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", encrypted_Indel_vector_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_sorted_BED_fp))
		{
			fprintf(stderr, "Skipping %s, no impact signals for this chromosome.\n", var_chr_ids->at(i_chr));
			continue;
		}

		// Load the EOI regions.
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
			continue;
		}

		char cur_chr_variant_signal_fp[1000];
		sprintf(cur_chr_variant_signal_fp, "%s/variant_signal_%s.bin.gz", encrypted_Indel_vector_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_variant_signal_fp, var_chr_ids->at(i_chr));
			exit(0);
		}

		// Load the annotation signals.
		fprintf(stderr, "Loading indel annotation signals.\n");
		
		// Allocate the per indel per l_indel selector.
		unsigned long long* indel_annotation_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		memset(indel_annotation_signal, 0, sizeof(unsigned long long) * cur_chr_EOI_covg);

		// Load the annotation signals.
		FILE* f_annotation_signal = open_f(cur_chr_annotation_signal_fp, "rb");
		unsigned int n_read_annotations = fread(indel_annotation_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_annotation_signal);
		if(n_read_annotations != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Could not read annotation signal from %s\n", cur_chr_annotation_signal_fp);
			exit(0);
		}

		fprintf(stderr, "Loaded %d long annotation signal.\n", n_read_annotations);
		close_f(f_annotation_signal, cur_chr_annotation_signal_fp);

		// Load the variant signals.
		fprintf(stderr, "Loading per length indel variant signals.\n");
		unsigned int* indel_genotype_signal = new unsigned int[cur_chr_EOI_covg + 2];
		memset(indel_genotype_signal, 0, sizeof(unsigned int) * cur_chr_EOI_covg);

		FILE* f_variant_signal = open_f(cur_chr_variant_signal_fp, "rb");
		unsigned int n_read_genotypes = fread(indel_genotype_signal, sizeof(unsigned int), cur_chr_EOI_covg, f_variant_signal);
		fprintf(stderr, "Loaded %d long genotype signal.\n", n_read_genotypes);
		if (n_read_genotypes != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Could not read genotype signal from %s\n", cur_chr_variant_signal_fp);
			exit(0);
		}
		close_f(f_variant_signal, cur_chr_variant_signal_fp);

		// Allocate the multiplied signal.
		fprintf(stderr, "Allocating annotated variants signal.\n");
		unsigned long long* cross_variant_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		memset(cross_variant_signal, 0, sizeof(unsigned long long) * cur_chr_EOI_covg);

		// Simply multiply every integer on variant and annotation signals.
		for (int i = 0; i < cur_chr_EOI_covg; i++)
		{
			cross_variant_signal[i] = indel_genotype_signal[i] * indel_annotation_signal[i];

			/*if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				if (cross_variant_signal[i] > 0)
				{
					char cross_seq[100];
					unpack_neigh_nucs(cross_variant_signal[i], cross_seq, 3, 6);

					char annot_seq[100];
					unpack_neigh_nucs(indel_annotation_signal[i], annot_seq, 3, 6);

					unsigned int impact_score = cross_variant_signal[i] >> 18;
					fprintf(stderr, "%s:%d [geno: %u] (Neighseq: %s - %s ;; Impact Score: %u)\n", var_chr_ids->at(i_chr), i,
							indel_genotype_signal[i],
							cross_seq, annot_seq, 
							impact_score);
				}
			}*/
		} // i loop.

		// Save the files.
		char cur_chr_anno_var_signal_op_fp[1000];
		sprintf(cur_chr_anno_var_signal_op_fp, "%s/annotated_variant_signal_%s.bin.gz", op_dir, var_chr_ids->at(i_chr));
		fprintf(stderr, "Saving annotated variant signal to %s.\n", cur_chr_anno_var_signal_op_fp);

		FILE* f_cur_chr_anno_var_signal_op = open_f(cur_chr_anno_var_signal_op_fp, "wb");
		fwrite(cross_variant_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_cur_chr_anno_var_signal_op);
		close_f(f_cur_chr_anno_var_signal_op, cur_chr_anno_var_signal_op_fp);

		// Save the target regions.
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, var_chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, sorted_cur_chr_EOI_regs);

		delete[] cross_variant_signal;
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
} // multiply_Indel_variant_and_annotation_signals



