#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svat_variation_tools.h"
#include "svat_genomics_coords.h"
#include "svat_annot_region_tools.h"
#include "svat_genome_sequence_tools.h"
#include "svat_gff_utils.h"
#include "file_utils.h"
#include "svat_genomics_coords.h"
#include "svat_ansi_string.h"
#include "svat_nucleotide.h"
#include <string.h>
#include <ctype.h>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svat_nomenclature.h"
#include "svat_rng.h"
#include "svat_seed_manager.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

bool __DUMP_VARIATION_TOOLS_MSGS__ = false;

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

double get_recomb_rate_difference_per_regions(t_restr_annot_region_list* restr_recomb_rate_regs, 
											t_annot_region* reg1, t_annot_region* reg2)
{
	if (!t_string::compare_strings(reg1->chrom, reg2->chrom))
	{
		return(1);
	}

	int i_chr = t_string::get_i_str(restr_recomb_rate_regs->chr_ids, reg1->chrom);

	t_annot_region* left_reg = NULL;
	t_annot_region* right_reg = NULL;
	if (reg1->start < reg2->start)
	{
		left_reg = reg1;
		right_reg = reg2;
	}
	else
	{
		left_reg = reg2;
		right_reg = reg1;
	}
	
	// Overlap the left region:
	t_annot_region* left_oling_recomb_reg = NULL;
	if (left_reg->start <= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->start)
	{
		left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0);
	}
	else if(left_reg->start >= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->end)
	{
		left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->back();
	}
	else
	{
		int reg_i = locate_posn_region_per_region_starts(left_reg->start, restr_recomb_rate_regs->regions_per_chrom[i_chr], 0, restr_recomb_rate_regs->regions_per_chrom[i_chr]->size());
		while (reg_i > 0 &&
			reg_i < restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_end >= left_reg->start)
		{
			reg_i--;
		} // reg_i loop.

		while (reg_i >= 0 &&
			reg_i < restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_start <= left_reg->end)
		{
			int ol_start = MAX(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->start, left_reg->start);
			int ol_end = MIN(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->end, left_reg->end);
			if (ol_end >= ol_start)
			{
				// Found overlap, set it.
				left_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Found right reg @ %s:%d:: %s:%d-%d\n",
						left_reg->chrom, left_reg->start,
						left_oling_recomb_reg->chrom, left_oling_recomb_reg->start, left_oling_recomb_reg->end);
				}

				break;
			}
			reg_i++;
		} // reg_i loop.
	}

	// Overlap the right region:
	t_annot_region* right_oling_recomb_reg = NULL;
	if (right_reg->start <= restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0)->start)
	{
		right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(0);
	}
	else if (right_reg->start >= restr_recomb_rate_regs->regions_per_chrom[i_chr]->back()->end)
	{
		right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->back();
	}
	else
	{
		int reg_i = locate_posn_region_per_region_starts(right_reg->start, restr_recomb_rate_regs->regions_per_chrom[i_chr], 0, restr_recomb_rate_regs->regions_per_chrom[i_chr]->size());
		while (reg_i > 0 &&
			reg_i < restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_end >= right_reg->start)
		{
			reg_i--;
		} // reg_i loop.

		while (reg_i >= 0 &&
			reg_i < restr_recomb_rate_regs->regions_per_chrom[i_chr]->size() &&
			restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->sort_info->cumulative_sorted_start <= right_reg->end)
		{
			int ol_start = MAX(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->start, right_reg->start);
			int ol_end = MIN(restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i)->end, right_reg->end);
			if (ol_end >= ol_start)
			{
				// Found overlap, set it.
				right_oling_recomb_reg = restr_recomb_rate_regs->regions_per_chrom[i_chr]->at(reg_i);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Found right reg @ %s:%d:: %s:%d-%d\n",
						right_reg->chrom, right_reg->start,
						right_oling_recomb_reg->chrom, right_oling_recomb_reg->start, right_oling_recomb_reg->end);
				}

				break;
			}
			reg_i++;
		} // reg_i loop.
	}

	double recomb_diff = -1;
	if (left_oling_recomb_reg != NULL &&
		right_oling_recomb_reg != NULL)
	{		
		double* first_reg_recomb_info = (double*)(left_oling_recomb_reg->data);
		double* last_reg_recomb_info = (double*)(right_oling_recomb_reg->data);
		recomb_diff = last_reg_recomb_info[0] - first_reg_recomb_info[1];

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "delta recomb(%s:%d, %s:%d)::[%s:%d-%d (%.4f) ---> %s:%d-%d (%.4f)]\n",
				reg1->chrom, reg1->start, reg2->chrom, reg2->start,
				left_oling_recomb_reg->chrom, left_oling_recomb_reg->start, left_oling_recomb_reg->end, first_reg_recomb_info[1],
				right_oling_recomb_reg->chrom, right_oling_recomb_reg->start, right_oling_recomb_reg->end, last_reg_recomb_info[0]);
		}
	}
	else
	{
		fprintf(stderr, "Could not match one of the regions for delta recomb(%s:%d, %s:%d)\n",
			reg1->chrom, reg1->start, reg2->chrom, reg2->start);
		exit(0);
	}

	return(recomb_diff);
}

// https://en.wikipedia.org/wiki/Centimorgan
double get_recomb_prob_per_delta_cM(double delta_cM)
{
	double recomb_prob = 0.5 * (1 - exp(-2 * delta_cM / 100));
	return(recomb_prob);
}

vector<t_annot_region*>* load_recombination_rate_regions_per_map_file(char* recombination_rate_dir)
{
	char recomb_rate_chr_ids_fp[1000];
	sprintf(recomb_rate_chr_ids_fp, "%s/chr_ids.txt", recombination_rate_dir);
	vector<char*>* recomb_chr_ids = buffer_file(recomb_rate_chr_ids_fp);
	if (recomb_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load chromosome id's from recombination directory @ %s\n", recomb_rate_chr_ids_fp);
		exit(0);
	}

	fprintf(stderr, "Loaded %d chr ids for the recombination rate.\n", recomb_chr_ids->size());
	vector<t_annot_region*>* recomb_rate_regs = new vector<t_annot_region*>();
	for (int i_chr = 0; i_chr < recomb_chr_ids->size(); i_chr++)
	{
		char cur_chr_recomb_rate_fp[1000];
		sprintf(cur_chr_recomb_rate_fp, "%s/%s.map", recombination_rate_dir, recomb_chr_ids->at(i_chr));

		fprintf(stderr, "Loading the recombination rates on %s from %s\n", recomb_chr_ids->at(i_chr), cur_chr_recomb_rate_fp);
		FILE* f_recomb_map = open_f(cur_chr_recomb_rate_fp, "r");
		vector<t_annot_region*>* cur_chr_recomb_regs = new vector<t_annot_region*>();
		while (1)
		{
			char* cur_line = getline(f_recomb_map);
			if (cur_line == NULL)
			{
				break;
			}

			// 22      .       0.000000        16051347
			char read_chr_id[100];
			double cM = 0;
			int cur_posn = 0;
			if (sscanf(cur_line, "%s %*s %lf %d", read_chr_id, &cM, &cur_posn) != 3)
			{
				fprintf(stderr, "Could not parse: %s\n", cur_line);
				exit(0);
			}

			t_annot_region* cur_recomb_reg = get_empty_region();
			cur_recomb_reg->chrom = t_string::copy_me_str(read_chr_id);
			cur_recomb_reg->start = cur_posn;
			cur_recomb_reg->end = cur_posn;
			cur_recomb_reg->dbl_score = cM;
			cur_recomb_reg->strand = '+';

			recomb_rate_regs->push_back(cur_recomb_reg);
			cur_chr_recomb_regs->push_back(cur_recomb_reg);
			delete[] cur_line;
		} // file reading loop.
		fclose(f_recomb_map);

		fprintf(stderr, "Resetting the coordinates and the recombination info for the %d recombination regions.\n", cur_chr_recomb_regs->size());

		// Fix the coordinates by neighboring regions positions.
		sort(cur_chr_recomb_regs->begin(), cur_chr_recomb_regs->end(), sort_regions);

		for (int i_reg = 1; i_reg < cur_chr_recomb_regs->size(); i_reg++)
		{
			cur_chr_recomb_regs->at(i_reg)->start = cur_chr_recomb_regs->at(i_reg - 1)->end+1;

			// Setup the recombination rates for this region.
			double* cur_reg_recomb_info = new double[2];
			cur_chr_recomb_regs->at(i_reg)->data = cur_reg_recomb_info;
			cur_reg_recomb_info[0] = cur_chr_recomb_regs->at(i_reg)->dbl_score;
			cur_reg_recomb_info[1] = cur_chr_recomb_regs->at(i_reg-1)->dbl_score;
		} // i_reg loop.
	} // i_chr loop.

	fprintf(stderr, "Loaded %d recombination regions in total.\n", recomb_rate_regs->size());

	return(recomb_rate_regs);
}

vector<t_annot_region*>* load_variant_signal_regions_wrapper(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp)
{
	vector<t_annot_region*>* var_sig_regs = NULL;

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	if (sample_ids == NULL)
	{
		fprintf(stderr, "Could not load the sample ids from %s", sample_ids_list_fp);
		exit(0);
	}

	if (t_string::compare_strings(geno_sig_regs_BED_fp, "stdin") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed") || 
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".txt.gz") ||
		t_string::ends_with(geno_sig_regs_BED_fp, ".bed.gz"))
	{
		var_sig_regs = load_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else if (t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".matbed.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedmat.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".bedsig.gz") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed") ||
			t_string::ends_with(geno_sig_regs_BED_fp, ".sigbed.gz"))
	{
		var_sig_regs = load_binarized_variant_genotype_signal_regions(geno_sig_regs_BED_fp, sample_ids);
	}
	else
	{
		fprintf(stderr, "%s(%d): Could not identify variant signal regions file type: %s\n", __FILE__, __LINE__, geno_sig_regs_BED_fp);
		exit(0);
	}

	t_string::clean_string_list(sample_ids);

	return(var_sig_regs);
}

void merge_samples_per_matching_genotype_signal_regions(char* sample1_regs_fp, char* sample1_list_fp, 
														char* sample2_regs_fp, char* sample2_list_fp, bool match_region_names, char* op_fp)
{
	vector<char*>* sample1_ids = buffer_file(sample1_list_fp);
	vector<t_annot_region*>* geno1_sig_regs = load_variant_signal_regions_wrapper(sample1_regs_fp, sample1_list_fp);

	vector<char*>* sample2_ids = buffer_file(sample2_list_fp);
	vector<t_annot_region*>* geno2_sig_regs = load_variant_signal_regions_wrapper(sample2_regs_fp, sample2_list_fp);

	fprintf(stderr, "Merging %d regions of %d samples in %s and %d regions of %d samples and saving to %s\n",
		geno1_sig_regs->size(), sample1_ids->size(), sample1_regs_fp, 
		geno2_sig_regs->size(), sample2_ids->size(), sample2_regs_fp, 
		op_fp);

	vector<t_annot_region*>* intersects = intersect_annot_regions(geno1_sig_regs, geno2_sig_regs, true);
	fprintf(stderr, "Processing %d matching regions.\n", intersects->size());

	vector<t_annot_region*>* merged_sample_genotype_regions = new vector<t_annot_region*>();
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* geno1_reg = int_info->src_reg;
		t_annot_region* geno2_reg = int_info->dest_reg;

		if (geno1_reg->start == geno2_reg->start &&
			geno1_reg->end == geno2_reg->end)
		{
			if (!match_region_names ||
				t_string::compare_strings(geno1_reg->name, geno2_reg->name))
			{
				t_annot_region* new_reg = get_empty_region();
				new_reg->chrom = t_string::copy_me_str(geno1_reg->chrom);
				new_reg->start = geno1_reg->start;
				new_reg->end = geno1_reg->end;
				new_reg->strand = geno1_reg->strand;

				if (match_region_names)
				{
					new_reg->name = t_string::copy_me_str(geno1_reg->name);
				}

				void** reg1_info = (void**)(geno1_reg->data);
				char* reg1_geno_sig = (char*)(reg1_info[0]);
				void** reg2_info = (void**)(geno2_reg->data);
				char* reg2_geno_sig = (char*)(reg2_info[0]);

				char* pooled_geno_sig = new char[sample1_ids->size() + sample2_ids->size() + 2];
				for (int i_s = 0; i_s < sample1_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s] = reg1_geno_sig[i_s];
				}

				for (int i_s = 0; i_s < sample2_ids->size(); i_s++)
				{
					pooled_geno_sig[i_s + sample1_ids->size()] = reg2_geno_sig[i_s];
				}

				// Copy the region info.
				void** new_reg_info = new void*[2];
				new_reg_info[0] = pooled_geno_sig;
				new_reg->data = new_reg_info;

				merged_sample_genotype_regions->push_back(new_reg);
			} // region name check.
		} // coordinate check.
	} // i_int loop.

	vector<char*>* merged_sample_ids = new vector<char*>();
	merged_sample_ids->insert(merged_sample_ids->end(), sample1_ids->begin(), sample1_ids->end());
	merged_sample_ids->insert(merged_sample_ids->end(), sample2_ids->begin(), sample2_ids->end());

	fprintf(stderr, "Saving the %d regions with %d samples to %s\n", merged_sample_genotype_regions->size(), merged_sample_ids->size(), op_fp);

	// Write the sample ids.
	char merged_sample_ids_fp[1000];
	sprintf(merged_sample_ids_fp, "%s_samples.list", op_fp);
	FILE* f_merged_sample_ids = open_f(merged_sample_ids_fp, "w");
	for (int i_s = 0; i_s < merged_sample_ids->size(); i_s++)
	{
		fprintf(f_merged_sample_ids, "%s\n", merged_sample_ids->at(i_s));
	} // i_s loop.
	fclose(f_merged_sample_ids);

	binarize_variant_genotype_signal_regions(merged_sample_genotype_regions, NULL, merged_sample_ids, op_fp);
}

void extract_genotype_signals_per_region_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* regions_BED_fp, char* op_fp)
{
	vector<t_annot_region*>* roi_regs = load_BED(regions_BED_fp);
	for (int i_reg = 0; i_reg < roi_regs->size(); i_reg++)
	{
		roi_regs->at(i_reg)->data = NULL;
	} // i_reg loop.
	fprintf(stderr, "Extracting genotype signals for %d ROIs.\n", roi_regs->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", geno_sig_regs->size(), sample_ids->size());

	fprintf(stderr, "Intersecting ROI's with genotype signal regions.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(roi_regs, geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* cur_roi_reg = int_info->src_reg;
		t_annot_region* cur_geno_sig_reg = int_info->dest_reg;

		if (t_string::compare_strings(cur_roi_reg->chrom, cur_geno_sig_reg->chrom) &&
			cur_roi_reg->start == cur_geno_sig_reg->start &&
			cur_roi_reg->end == cur_geno_sig_reg->end)
		{
			if (cur_roi_reg->data != NULL)
			{
				fprintf(stderr, "Sanity check failed: Already assigned: %s:%d-%d\n", cur_geno_sig_reg->chrom, cur_geno_sig_reg->start, cur_geno_sig_reg->end);
				exit(0);
			}

			cur_roi_reg->data = cur_geno_sig_reg->data;
		}
	} // i_int loop.

	// Count the number of ROI with signal regions.
	vector<t_annot_region*>* roi_regs_w_signals = new vector<t_annot_region*>();
	int n_roi_w_signal = 0;
	for (int i_reg = 0; i_reg < roi_regs->size(); i_reg++)
	{
		if (roi_regs->at(i_reg)->data != NULL)
		{
			roi_regs_w_signals->push_back(roi_regs->at(i_reg));
		}
	} // i_reg loop.
	fprintf(stderr, "Matched genotype signals to %d regions.\n", roi_regs_w_signals->size());

	// Save the regions with signals on them.
	binarize_variant_genotype_signal_regions(roi_regs_w_signals, NULL, sample_ids, op_fp);
}

// Following function extracts the genotype signals for samples from an input VCF file. It focuses on a chromosome and a subset of regions to decrease memory usage.
void extract_genotype_signals_per_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
										char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
										char* var_regions_BED_fp,			// Regions to focus on while extracting.
										char* chr_id_2_process,				// Chromosome to process.
										char* bin_seq_dir,					// Sequences directory.
										bool match_ref_alleles_flag,		// This is the flag that tells to sanity-check the ref allele.
										bool match_region_names_flag,		// This is a flag to enforce matching of region names.
										bool haplotype_specific_encoding,	// Do we want to use phasing information in encoding? (i.e., haplotype specific information)
										char* op_fp)						// Output file path.
{
	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(0);
	}

	vector<t_annot_region*>* all_var_regions = load_BED(var_regions_BED_fp);

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", vcf_sample_ids->size());

	vector<char*>* chr_ids = get_chr_ids(all_var_regions);

	vector<t_annot_region*>* var_regions = NULL;
	int chrom_i_2_process = t_string::get_i_str(chr_ids, chr_id_2_process);
	if (chrom_i_2_process == chr_ids->size())
	{
		var_regions = all_var_regions;
	}
	else
	{
		var_regions = get_regions_per_chromosome(all_var_regions, chr_ids->at(chrom_i_2_process));
	}

	chr_ids = get_chr_ids(var_regions);

	fprintf(stderr, "Extracting signals on %d regions.\n", var_regions->size());

	if (haplotype_specific_encoding)
	{
		fprintf(stderr, "Using haplotype specific encoding.\n");
	}
	else
	{
	}

	char** per_chr_seq = new char*[chr_ids->size() + 2];

	if (match_ref_alleles_flag)
	{
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			int l_chrom = 0;
			fprintf(stderr, "Loading %s.\n", chr_ids->at(i_chr));
			char cur_chr_seq_fp[1000];
			sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));

			if (check_file(cur_chr_seq_fp))
			{
				per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
			}
			else
			{
				sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, chr_ids->at(i_chr));

				if (check_file(cur_chr_seq_fp))
				{
					per_chr_seq[i_chr] = load_binary_sequence_file(cur_chr_seq_fp, l_chrom);
				}
				else
				{
					fprintf(stderr, "Could not load the sequence for %s\n", chr_ids->at(i_chr));
					exit(0);
				}
			}
		} // i_chr loop.
	}
	else
	{
		fprintf(stderr, "Skipping ref genome loading.\n");
	}
	  // Set the genotype signal array on all regions.
	for (int i_v = 0; i_v < var_regions->size(); i_v++)
	{
		char* pooled_var_alleles = var_regions->at(i_v)->name;

		void** cur_reg_info = new void*[2];
		cur_reg_info[0] = NULL;
		cur_reg_info[1] = NULL;
		var_regions->at(i_v)->data = cur_reg_info;
	} // i_v loop.

	FILE* f_vcf = open_f(vcf_fp, "r");

	char* buff = new char[100 * 1000];
	vector<t_annot_region*>* vcf_regs = new vector<t_annot_region*>();
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

		if (vcf_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region.           \r", vcf_regs->size());
		}

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
		char chrom[100];
		char posn_str[100];
		char id[100];
		char ref[100];
		char alt[100];
		char qual[100];
		char filter[100];
		char info[10000];
		char format[100];

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(format, buff);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int GT_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];
		while (t_string::get_next_token(format, cur_tok, l_format, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, "GT"))
			{
				GT_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "Found GT entry @ %d\n", GT_entry_i);
		}

		if (GT_entry_i == -1)
		{
			fprintf(stderr, "Could not find GT entry in %s, skipping\n", cur_line);
			delete[] cur_line;
			continue;
		}

		if (GT_entry_i != 0)
		{
			fprintf(stderr, "Format string is not as expected: %s\n", format);
		}
		
		// Make sure we have the chromosome id match and we are using alleles.
		char* cur_var_geno_sig = new char[vcf_sample_ids->size() + 2];

		// Initialize everything to not set: This contains the variant signals for each individual.
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			cur_var_geno_sig[i_s] = -1;
		} // i_s loop.

		bool correctly_parsed_all_genotypes = true;
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);

			//if (buff[3] != 0)
			//{
			//	fprintf(stderr, "Failed to read the %d. genotype correctly for line:%s\n", i_s, cur_line);
			//	exit(0);
			//}

			// Check if all the variants are bi-allelic, only.
			if ((buff[0] != '.' && buff[0] != '.' && buff[0] != '0' && buff[0] != '1') ||
				(buff[1] != '|' && buff[1] != '/') ||
				(buff[2] != '0' && buff[2] != '1' && buff[2] != '.' && buff[2] != '.'))
			{
				correctly_parsed_all_genotypes = false;
				fprintf(stderr, "Failed to read the %d. genotype correctly for entry %s in line: %s\n", i_s, buff, cur_line);
				//exit(0);
			}

			if (haplotype_specific_encoding)
			{
				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[1] == '|')
				{
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
				else
				{
					// Randomly assign the haplotypes?
					char geno0_val = (char)(buff[0] - '0');
					char geno1_val = (char)(buff[2] - '0');
					char cur_geno = 2 * geno0_val + geno1_val;

					cur_var_geno_sig[i_s] = cur_geno;
				}
			}
			else
			{
				cur_var_geno_sig[i_s] = 0;

				if (buff[0] == '.' || buff[2] == '.')
				{
					cur_var_geno_sig[i_s] = -1;
				}
				else if (buff[0] == '1' && buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 2;
				}
				else if (buff[0] == '1' || buff[2] == '1')
				{
					cur_var_geno_sig[i_s] = 1;
				}
			}

			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(0);
		}

		if (correctly_parsed_all_genotypes)
		{
			// Intersect with the variant regions.
			t_annot_region* cur_vcf_reg = get_empty_region();
			cur_vcf_reg->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(cur_vcf_reg->chrom);
			cur_vcf_reg->start = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);

			int l_ref_allele = t_string::string_length(ref);
			cur_vcf_reg->end = cur_vcf_reg->start + l_ref_allele - 1;
			cur_vcf_reg->name = t_string::copy_me_str(id);
			cur_vcf_reg->strand = '+';

			// Set the signal to data.
			void** vcf_reg_info = new void*[2];

			// 1st data: signals.
			vcf_reg_info[0] = cur_var_geno_sig;

			// 2nd data: alleles.
			char** ref_alt_alleles = new char*[2];
			ref_alt_alleles[0] = t_string::copy_me_str(ref);
			ref_alt_alleles[1] = t_string::copy_me_str(alt);
			vcf_reg_info[1] = ref_alt_alleles;

			// If matching ref alleles is requested, check the alleles.
			if (l_ref_allele < 100 &&
				match_ref_alleles_flag)
			{
				int i_chr_cur_line = t_string::get_i_str(chr_ids, cur_vcf_reg->chrom);
				for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
				{
					if (toupper(per_chr_seq[i_chr_cur_line][i]) != toupper(ref[i - cur_vcf_reg->start]))
					{
						fprintf(stderr, "Could not match the reference allele to chromosome sequence @ %s:%d: %c, %c\n",
							cur_vcf_reg->chrom, cur_vcf_reg->start,
							toupper(per_chr_seq[i_chr_cur_line][i]), toupper(ref[i - cur_vcf_reg->start]));

						exit(0);
					}
				} // i loop.

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "Match @ %s:%d-%d (%s): \n%s\n", cur_vcf_reg->chrom, cur_vcf_reg->start, cur_vcf_reg->end, cur_vcf_reg->name, ref);
					for (int i = cur_vcf_reg->start; i <= cur_vcf_reg->end; i++)
					{
						fprintf(stderr, "%c", toupper(per_chr_seq[i_chr_cur_line][i]));
					} // i loop.
					fprintf(stderr, "\n");
				}
			} // allele check is done for variants smaller than 100 bp.

			cur_vcf_reg->data = vcf_reg_info;

			vcf_regs->push_back(cur_vcf_reg);
		} // genotype parsing check.

		delete[] cur_line;
	} // vcf file reading loop.

	vector<t_annot_region*>* intersects = intersect_annot_regions(var_regions, vcf_regs, true);
	fprintf(stderr, "Processing %d intersections using coordinate matching.               \n", intersects->size());

	if (match_region_names_flag)
	{
		fprintf(stderr, "Enforcing substring based name matching.\n");
	}

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* pooled_var_reg = cur_int_info->src_reg;
		t_annot_region* vcf_reg = cur_int_info->dest_reg;

		bool coord_check_pass = true;
		if (pooled_var_reg->start != vcf_reg->start ||
			pooled_var_reg->end != vcf_reg->end)
		{
			coord_check_pass = false;
		}

		// Do we need to have name matching?
		bool name_check_pass = true;
		if (match_region_names_flag &&
			pooled_var_reg->name != NULL &&
			vcf_reg->name != NULL)
		{
			int i_char = 0;

			// If there is a dot for the name, skip.
			if (!t_string::compare_strings(pooled_var_reg->name, ".") &&
				!t_string::compare_strings(vcf_reg->name, ".") &&
				!t_string::compare_substrings_ci(pooled_var_reg->name, vcf_reg->name, i_char) &&
				!t_string::compare_substrings_ci(vcf_reg->name, pooled_var_reg->name, i_char))
			{
				name_check_pass = false;
			}
		}

		if (name_check_pass &&
			coord_check_pass)
		{
			int i_chr = t_string::get_i_str(chr_ids, pooled_var_reg->chrom);
			char* chr_seq = per_chr_seq[i_chr];

			// Get the genotype signal for the vcf region.
			void** vcf_reg_info = (void**)(vcf_reg->data);
			char* cur_vcf_reg_geno_sig = (char*)(vcf_reg_info[0]);

			// Set the genotype signal.
			void** cur_var_reg_info = (void**)(pooled_var_reg->data);
			cur_var_reg_info[0] = cur_vcf_reg_geno_sig;
		}

		delete cur_int_info;
	} // i_int loop.
	delete_annot_regions(vcf_regs);

	// Copy the current haplotype alleles.
	for (int i_reg = 0; i_reg < var_regions->size(); i_reg++)
	{
		void** cur_var_reg_info = (void**)(var_regions->at(i_reg)->data);

		if (cur_var_reg_info[0] == NULL)
		{
			char* cur_var_geno_sig = new char[vcf_sample_ids->size()];
			for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
			{
				cur_var_geno_sig[i_s] = -1;
			} // i_s loop.

			cur_var_reg_info[0] = cur_var_geno_sig;
		}
	} // i_reg loop.

	// Save the binarized signals.
	binarize_variant_genotype_signal_regions(var_regions, NULL, vcf_sample_ids, op_fp);

	//// Dump the plain matrix, too.
	//dump_geno_sig_regs_plain(var_regions, vcf_sample_ids, "op.bed");
	
	//// Now dump the variant signal profiles.
	//FILE* f_geno_sig = open_f("op.bed", "w");
	//for (int i_v = 0; i_v < var_regions->size(); i_v++)
	//{
	//	fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
	//		var_regions->at(i_v)->chrom, 
	//		translate_coord(var_regions->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
	//		translate_coord(var_regions->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		var_regions->at(i_v)->name);

	//	void** cur_var_reg_info = (void**)(var_regions->at(i_v)->data);
	//	int* cur_var_signal = (int*)(cur_var_reg_info[0]);
	//	for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_geno_sig, "\t%d", cur_var_signal[i_s]);
	//	} // i_s loop.

	//	fprintf(f_geno_sig, "\n");
	//} // i_v loop.
	//fclose(f_geno_sig);
}

void convert_genotype_signal_regs_2_VCF(char* geno_sigs_reg_fp, char* sample_ids_list_fp, bool phase_op, char* op_VCF_fp)
{
	fprintf(stderr, "Writing the VCF from %s using phased=%d genotypes.\n", geno_sigs_reg_fp, phase_op);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples.\n", sample_ids->size());

	vector<t_annot_region*>* geno_sig_regs = load_variant_genotype_signal_regions(geno_sigs_reg_fp, sample_ids);

	fprintf(stderr, "Loaded %d regions.\n", geno_sig_regs->size());
	t_restr_annot_region_list* restr_geno_sig_regs = restructure_annot_regions(geno_sig_regs);

	FILE* f_VCF = open_f(op_VCF_fp, "w");

	// Write the header.
	fprintf(f_VCF, "##fileformat=VCFv4.2\n");
	fprintf(f_VCF, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < sample_ids->size(); i_s++)
	{
		fprintf(f_VCF, "\t%s", sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_VCF, "\n");

	// Write the genotypes.
	for (int i_chr = 0; i_chr < restr_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			if (i_reg % 1000 == 0)
			{
				fprintf(stderr, "Converting %d/%d variant.             \r", i_reg, cur_chr_input_geno_sig_regs->size());
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			if (toks->size() < 3)
			{
				fprintf(stderr, "Could not parse alleles from %s; make sure the region names are formatted as: [Variant ID]_[Ref allele]_[Alt allele]\n", cur_chr_input_geno_sig_regs->at(i_reg)->name);
				exit(0);
			}

			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_VCF, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < sample_ids->size(); i_s++)
			{
				if (!phase_op)
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					int geno = cur_hap0 + cur_hap1;

					if (geno == 0)
					{
						fprintf(f_VCF, "\t0/0");
					}
					else if (geno == 1)
					{
						fprintf(f_VCF, "\t0/1");
					}
					else if (geno == 2)
					{
						fprintf(f_VCF, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(0);
					}
				}
				else
				{
					int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
					int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

					char geno_str[10];
					strcpy(geno_str, "0|0");

					if (cur_hap0 == 1)
					{
						geno_str[0] = '1';
					}

					if (cur_hap1 == 1)
					{
						geno_str[2] = '1';
					}

					fprintf(f_VCF, "\t%s", geno_str);

					int geno = cur_hap0 + cur_hap1;
					if (geno > 2)
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(0);
					}
				}
			} // input genotype existence check.

			fprintf(f_VCF, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_VCF, op_VCF_fp);
}


int get_genotype_per_haplocoded_genotype(char haplocoded_geno)
{
	int all1 = get_allele_per_haplotype(haplocoded_geno, 0);
	int all2 = get_allele_per_haplotype(haplocoded_geno, 1);
	return(all1 + all2);
}

int get_allele_per_haplotype(char geno, int hap_i)
{
	int allele = ((geno & (1 << hap_i)) >> hap_i);
	return(allele);
}

void dump_haplo_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, int haplo_2_extract, const char* op_fp)
{
	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);


		void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
		char* cur_var_signal = (char*)(cur_var_reg_info[0]);
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			char cur_sample_haplo_str[10];
			memset(cur_sample_haplo_str, 0, 10);

			if (haplo_2_extract == 2)
			{
				cur_sample_haplo_str[1] = '|';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 3)
			{
				cur_sample_haplo_str[1] = '\t';
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
				cur_sample_haplo_str[2] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}
			else if (haplo_2_extract == 0)
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 0) + '0';
			}
			else
			{
				cur_sample_haplo_str[0] = (char)get_allele_per_haplotype(cur_var_signal[i_s], 1) + '0';
			}

			fprintf(f_geno_sig, "\t%s", cur_sample_haplo_str);
		} // i_s loop.

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void dump_geno_sig_regs_plain(vector<t_annot_region*>* geno_sig_regs, vector<char*>* vcf_sample_ids, bool dump_geno_sig_regs_plain, const char* op_fp)
{
	if (dump_geno_sig_regs_plain)
	{
		fprintf(stderr, "Saving regions only.\n");
	}

	// Now dump the variant signal profiles.
	FILE* f_geno_sig = open_f(op_fp, "w");
	for (int i_v = 0; i_v < geno_sig_regs->size(); i_v++)
	{
		char cur_reg_name[1000];
		if (geno_sig_regs->at(i_v)->name == NULL)
		{
			strcpy(cur_reg_name, "NONAME");
		}
		else
		{
			strcpy(cur_reg_name, geno_sig_regs->at(i_v)->name);
		}

		fprintf(f_geno_sig, "%s\t%d\t%d\t%s",
			geno_sig_regs->at(i_v)->chrom,
			translate_coord(geno_sig_regs->at(i_v)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(geno_sig_regs->at(i_v)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			cur_reg_name);

		if (!dump_geno_sig_regs_plain)
		{
			void** cur_var_reg_info = (void**)(geno_sig_regs->at(i_v)->data);
			char* cur_var_signal = (char*)(cur_var_reg_info[0]);
			for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
			{
				fprintf(f_geno_sig, "\t%d", (int)cur_var_signal[i_s]);
			} // i_s loop.
		}

		fprintf(f_geno_sig, "\n");
	} // i_v loop.
	fclose(f_geno_sig);
}

void extract_genotype_signals_per_subsample_list(char* geno_sig_regs_BED_fp, char* sample_ids_list_fp, char* subsample_ids_list_fp, char* op_fp)
{
	vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", geno_sig_regs->size(), sample_ids->size());

	vector<char*>* geno_subsample_ids = buffer_file(subsample_ids_list_fp);
	if (geno_subsample_ids == NULL)
	{
		fprintf(stderr, "Could not load subsample id's list from %s\n", subsample_ids_list_fp);
		exit(0);
	}
	else
	{
		fprintf(stderr, "Extracting genotype signals for %d individuals.\n", geno_subsample_ids->size());
	}

	vector<int>* per_subsample_geno_sample_i = new vector<int>();
	for (int i_ss = 0; i_ss < geno_subsample_ids->size(); i_ss++)
	{
		int cur_ss_sample_i = t_string::get_i_str(sample_ids, geno_subsample_ids->at(i_ss));
		per_subsample_geno_sample_i->push_back(cur_ss_sample_i);
	} // i_sub_s loop.

	vector<t_annot_region*>* geno_sig_regs_w_subsample_signals = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Extracting %d. region            \r", i_reg);
		}

		void** cur_reg_info = (void**)(geno_sig_regs->at(i_reg)->data);
		char* cur_reg_geno_sigs = (char*)(cur_reg_info[0]);

		// Copy the signals.
		char* cur_copy_reg_geno_sigs = new char[geno_subsample_ids->size() + 2];
		for (int i_ss = 0; i_ss < geno_subsample_ids->size(); i_ss++)
		{			
			if (per_subsample_geno_sample_i->at(i_ss) < sample_ids->size())
			{
				cur_copy_reg_geno_sigs[i_ss] = cur_reg_geno_sigs[per_subsample_geno_sample_i->at(i_ss)];
			}
			else
			{
				cur_copy_reg_geno_sigs[i_ss] = -1;
			}
		} // i_ss loop.

		void** cur_copy_reg_info = new void*[2];		
		cur_copy_reg_info[0] = cur_copy_reg_geno_sigs;

		t_annot_region* copy_reg = duplicate_region(geno_sig_regs->at(i_reg));
		copy_reg->data = cur_copy_reg_info;
		geno_sig_regs_w_subsample_signals->push_back(copy_reg);
	} // i_reg loop.

	// Save the regions.
	binarize_variant_genotype_signal_regions(geno_sig_regs_w_subsample_signals, NULL, geno_subsample_ids, op_fp);
}

// Following can dump the supplied regions or load them then dump.
void binarize_variant_genotype_signal_regions(vector<t_annot_region*>* genotype_signal_regions, char* variant_geno_sig_regs_BED_fp, vector<char*>* geno_sample_ids, const char* op_fp)
{
	if (genotype_signal_regions == NULL)
	{
		if (t_string::compare_strings(variant_geno_sig_regs_BED_fp, "stdin") || check_file(variant_geno_sig_regs_BED_fp))
		{
			genotype_signal_regions = load_variant_genotype_signal_regions(variant_geno_sig_regs_BED_fp, geno_sample_ids);
		}		
		else
		{
			fprintf(stderr, "Signal regions are not supplied and could not load them using %s.\n", variant_geno_sig_regs_BED_fp);
			exit(0);
		}
	}

	vector<char*>* chr_ids = get_chr_ids(genotype_signal_regions);

	FILE* f_op = open_f(op_fp, "wb");
	int n_chrs = chr_ids->size();
	fwrite(&n_chrs, sizeof(int), 1, f_op);
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chr[1000];
		strcpy(cur_chr, chr_ids->at(i_chr));
		fwrite(cur_chr, sizeof(char), 1000, f_op);
	} // i_chr loop.

	int sample_size = geno_sample_ids->size();
	int n_regs = genotype_signal_regions->size();

	fwrite(&sample_size, sizeof(int), 1, f_op);
	fwrite(&n_regs, sizeof(int), 1, f_op);
	for (int i_reg = 0; i_reg < genotype_signal_regions->size(); i_reg++)
	{
		// Write the chromosome index.
		int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		fwrite(&i_chr, sizeof(int), 1, f_op);
		fwrite(&(reg_BED_start), sizeof(int), 1, f_op);
		fwrite(&(reg_BED_end), sizeof(int), 1, f_op);

		// Write the regions's name.
		if (genotype_signal_regions->at(i_reg)->name == NULL)
		{
			genotype_signal_regions->at(i_reg)->name = t_string::copy_me_str(".");
		}

		int l_reg_name_str = t_string::string_length(genotype_signal_regions->at(i_reg)->name);
		fwrite(&l_reg_name_str, sizeof(int), 1, f_op);
		fwrite(genotype_signal_regions->at(i_reg)->name, sizeof(char), l_reg_name_str, f_op);

		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
		fwrite(cur_reg_geno_sig, sizeof(char), geno_sample_ids->size(), f_op);
	} // i_reg loop.

	fprintf(stderr, "\nClosing the file.\n");
	close_f(f_op, op_fp);
	fprintf(stderr, "Finished writing. Doing a sanity check.\n");

	bool perform_check = true;

	if (!perform_check)
	{
		fprintf(stderr, "Skipping sanity check.\n");
		return;
	}

	clock_t cur_time = clock();
	vector<t_annot_region*>* loaded_sig_regs = load_binarized_variant_genotype_signal_regions(op_fp, geno_sample_ids);
	fprintf(stderr, "Loaded in %.4f seconds.\n", ((double)(clock() - cur_time)) / CLOCKS_PER_SEC);

	if (loaded_sig_regs->size() != genotype_signal_regions->size())
	{
		fprintf(stderr, "Dumping/Loading failed.\n");
		exit(0);
	}

	for (int i_reg = 0; i_reg < loaded_sig_regs->size(); i_reg++)
	{
		if (loaded_sig_regs->at(i_reg)->start != genotype_signal_regions->at(i_reg)->start ||
			loaded_sig_regs->at(i_reg)->end != genotype_signal_regions->at(i_reg)->end)
		{
			fprintf(stderr, "Non-match: %s:%d-%d, %s:%d-%d\n", 
					loaded_sig_regs->at(i_reg)->chrom, loaded_sig_regs->at(i_reg)->start, loaded_sig_regs->at(i_reg)->end, 
					genotype_signal_regions->at(i_reg)->chrom, genotype_signal_regions->at(i_reg)->start, genotype_signal_regions->at(i_reg)->end);
			exit(0);
		}

		// Go over the genotype signals.
		void** cur_reg_info = (void**)(genotype_signal_regions->at(i_reg)->data);
		char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);

		void** cur_loaded_reg_info = (void**)(loaded_sig_regs->at(i_reg)->data);
		char* cur_loaded_reg_geno_sig = (char*)(cur_loaded_reg_info[0]);

		// Write the signal values.
		for (int i_s = 0; i_s < geno_sample_ids->size(); i_s++)
		{
			if (cur_reg_geno_sig[i_s] != cur_loaded_reg_geno_sig[i_s])
			{
				fprintf(stderr, "Non-match: %d. sample genotype: %d, %d\n",
						cur_reg_geno_sig[i_s], cur_loaded_reg_geno_sig[i_s]);
				exit(0);
			}
		} // i_s loop.
	} // i_reg loop.

	fprintf(stderr, "Check success!\n");
}

vector<t_annot_region*>* load_binarized_variant_genotype_signal_regions(const char* bin_geno_sig_bed_fp, vector<char*>* geno_sample_ids)
{
	FILE* f_bin_geno_sig_regs = open_f(bin_geno_sig_bed_fp, "rb");

	// Load the chromosomes.
	int n_chrs = 0;
	fread(&n_chrs, sizeof(int), 1, f_bin_geno_sig_regs);
	vector<char*>* chr_ids = new vector<char*>();
	for (int i_chr = 0; i_chr < n_chrs; i_chr++)
	{
		char cur_chr[1000];
		fread(cur_chr, sizeof(char), 1000, f_bin_geno_sig_regs);
		chr_ids->push_back(t_string::copy_me_str(cur_chr));
	} // i_chr loop.

	fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());
	int sample_size = 0;
	fread(&sample_size, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading sample size of %d.\n", sample_size);

	if (sample_size != geno_sample_ids->size())
	{
		fprintf(stderr, "Sanity check failed: Sample sizes do not match: %d, %d\n", sample_size, geno_sample_ids->size());
		exit(0);
	}

	int n_regs = 0;
	fread(&n_regs, sizeof(int), 1, f_bin_geno_sig_regs);
	fprintf(stderr, "Reading %d regions.\n", n_regs);

	vector<t_annot_region*>* geno_sig_regs = new vector<t_annot_region*>();
	for (int i_reg = 0; i_reg < n_regs; i_reg++)
	{
		//int i_chr = t_string::get_i_str(chr_ids, genotype_signal_regions->at(i_reg)->chrom);
		//int reg_BED_start = translate_coord(genotype_signal_regions->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base);
		//int reg_BED_end = translate_coord(genotype_signal_regions->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base);

		int i_chr = 0;
		int reg_BED_start = 0;
		int reg_BED_end = 0;
		fread(&i_chr, sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_start), sizeof(int), 1, f_bin_geno_sig_regs);
		fread(&(reg_BED_end), sizeof(int), 1, f_bin_geno_sig_regs);

		// Read the region's name.
		int l_reg_name_str = 0;
		fread(&l_reg_name_str, sizeof(int), 1, f_bin_geno_sig_regs);
		char* cur_reg_name = new char[l_reg_name_str + 2];
		memset(cur_reg_name, 0, sizeof(char) * (l_reg_name_str + 2));
		fread(cur_reg_name, sizeof(char), l_reg_name_str, f_bin_geno_sig_regs);

		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chr_ids->at(i_chr));
		reg->start = translate_coord(reg_BED_start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		reg->end = translate_coord(reg_BED_end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		reg->strand = '+';
		reg->name = cur_reg_name;

		void** cur_reg_info = new void*[2];
		char* cur_reg_geno_sig = new char[sample_size + 2];
		fread(cur_reg_geno_sig, sizeof(char), sample_size, f_bin_geno_sig_regs);

		cur_reg_info[0] = cur_reg_geno_sig;
		cur_reg_info[1] = NULL;

		reg->data = cur_reg_info;

		geno_sig_regs->push_back(reg);
	} // i_reg loop.

	// Close the file.
	close_f(f_bin_geno_sig_regs, bin_geno_sig_bed_fp);

	return(geno_sig_regs);
}

vector<t_annot_region*>* load_variant_genotype_signal_regions(char* link_variant_genotype_signal_fp, vector<char*>* geno_sample_ids)
{
	vector<t_annot_region*>* link_var_regs = new vector<t_annot_region*>();
	char tok_buff[1000];
	FILE* f_link_variant_genotype_signal = open_f(link_variant_genotype_signal_fp, "r");

	if (f_link_variant_genotype_signal == NULL)
	{
		fprintf(stderr, "Could not open %s\n", link_variant_genotype_signal_fp);
		exit(0);
	}

	int i_reg = 0;
	while (1)
	{
		char* cur_reg_line = getline(f_link_variant_genotype_signal);
		if (cur_reg_line == NULL)
		{
			break;
		}
		else
		{
			i_reg++;
		}

		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Parsing genotype signal for %d. region.                  \r", i_reg);
		}

		//char* cur_reg_line = (char*)(link_var_regs->at(i_reg)->data);
		//t_string_tokens* toks = t_string::tokenize_by_chars(cur_reg_line, "\t");
		char* cur_var_geno_signal = new char[geno_sample_ids->size()];

		t_annot_region* cur_reg = get_empty_region();

		// Copy the name.
		int i_cur_char = 0;
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->chrom = t_string::copy_me_str(tok_buff);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->start = translate_coord(atoi(tok_buff), BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->end = translate_coord(atoi(tok_buff), BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_reg->strand = '+';

		t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char);
		cur_reg->name = t_string::copy_me_str(tok_buff);

		//fprintf(stderr, "%s\t%d\t%d: %s            \r", link_var_regs->at(i_reg)->chrom, link_var_regs->at(i_reg)->start, link_var_regs->at(i_reg)->end, link_var_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < geno_sample_ids->size(); i_s++)
		{
			if (t_string::get_next_token(cur_reg_line, tok_buff, 1000, "\t", i_cur_char) == false)
			{
				fprintf(stderr, "Could not parse the %d. sample's genotype signal entry.\n", i_s);
				exit(0);
			}

			cur_var_geno_signal[i_s] = (char)(atoi(tok_buff));
		} // i_tok loop.

		void** link_var_reg_info = new void*[3];
		link_var_reg_info[0] = cur_var_geno_signal;
		link_var_reg_info[1] = NULL;

		cur_reg->data = link_var_reg_info;

		link_var_regs->push_back(cur_reg);
	} // i_reg loop.

	close_f(f_link_variant_genotype_signal, link_variant_genotype_signal_fp);

	return(link_var_regs);
}

void get_SNV_splice_effect(t_annot_region* var_reg,
							t_annot_region* intronic_reg,
							int& splice_effect_type)
{
	// Init the effect type.
	splice_effect_type = -1;

	t_variant_info* cur_var_info = (t_variant_info*)(var_reg->annotation_info);
	t_intron_annotation_info* cur_intron_info = (t_intron_annotation_info*)(intronic_reg->annotation_info);
	char* intron_seq = cur_intron_info->sequence;
	int relative_var_posn_in_intron = var_reg->start - intronic_reg->start;

	int l_intron = intronic_reg->end - intronic_reg->start + 1;

	char splice_donor_motif[] = "GT";
	char splice_acceptor_motif[] = "AG";
	char* donor_acceptor_motifs[2];
	donor_acceptor_motifs[0] = splice_donor_motif;
	donor_acceptor_motifs[1] = splice_acceptor_motif;

	if (__DUMP_VARIATION_TOOLS_MSGS__)
	{
		fprintf(stderr, "Var: %s:%d-%d (%c, %c); Intron: %d-%d; rel. pos.: %d\n", var_reg->chrom, var_reg->start, var_reg->end, cur_var_info->ref_allele[0], cur_var_info->alt_allele[0], intronic_reg->start, intronic_reg->end, relative_var_posn_in_intron);
	}

	// Depending on the variant type, evaluate effect on splicing.
	if(cur_var_info->var_type == VAR_TYPE_SNV)
	{
		char var_alternate_allele = cur_var_info->alt_allele[0];

		// Left side check.
		if (relative_var_posn_in_intron <= 1)
		{
			if(intronic_reg->strand == '+')
			{
				// 5' end: This is a splice donor site.
				char cur_motif_str[10];
				strcpy(cur_motif_str, splice_donor_motif);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "@ 5' end (+): %c%c (%s); %c <-> %c\n", intron_seq[0], intron_seq[1], cur_motif_str, intron_seq[relative_var_posn_in_intron], cur_var_info->alt_allele[0]);
				}

				if (intron_seq[relative_var_posn_in_intron] != cur_var_info->alt_allele[0] &&
					nuc_2_num(intron_seq[0]) == nuc_2_num(cur_motif_str[0]) &&
					nuc_2_num(intron_seq[1]) == nuc_2_num(cur_motif_str[1]))
				{
					splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
				}
			}
			else
			{
				// 3' end: This is a splice acceptor site.
				char cur_motif_str[10];
				strcpy(cur_motif_str, splice_acceptor_motif);
				reverse_complement_seq(cur_motif_str);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "@ 5' end (-): %c%c (%s); %c <-> %c\n", intron_seq[0], intron_seq[1], cur_motif_str, intron_seq[relative_var_posn_in_intron], cur_var_info->alt_allele[0]);
				}

				if (intron_seq[relative_var_posn_in_intron] != cur_var_info->alt_allele[0] &&
					nuc_2_num(intron_seq[0]) == nuc_2_num(cur_motif_str[0]) &&
					nuc_2_num(intron_seq[1]) == nuc_2_num(cur_motif_str[1]))
				{
					splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
				}
			}
		} // Left side of the intron check.
		else if (relative_var_posn_in_intron >= l_intron - 2)
		{
			if (intronic_reg->strand == '+')
			{
				// 5' end: This is a splice acceptor site.
				char cur_motif_str[10];
				strcpy(cur_motif_str, splice_acceptor_motif);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "@ 3' end (+): %c%c (%s); %c <-> %c\n", intron_seq[l_intron - 2], intron_seq[l_intron - 1], cur_motif_str, intron_seq[relative_var_posn_in_intron], cur_var_info->alt_allele[0]);
				}

				if (intron_seq[relative_var_posn_in_intron] != cur_var_info->alt_allele[0] &&
					nuc_2_num(intron_seq[l_intron - 2]) == nuc_2_num(cur_motif_str[0]) &&
					nuc_2_num(intron_seq[l_intron - 1]) == nuc_2_num(cur_motif_str[1]))
				{
					splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
				}
			}
			else
			{
				// 3' end: This is a splice acceptor site.
				char cur_motif_str[10];
				strcpy(cur_motif_str, splice_donor_motif);
				reverse_complement_seq(cur_motif_str);

				if (__DUMP_VARIATION_TOOLS_MSGS__)
				{
					fprintf(stderr, "@ 3' end (-): %c%c (%s); %c <-> %c\n", intron_seq[l_intron - 2], intron_seq[l_intron - 1], cur_motif_str, intron_seq[relative_var_posn_in_intron], cur_var_info->alt_allele[0]);
				}

				if (intron_seq[relative_var_posn_in_intron] != cur_var_info->alt_allele[0] &&
					nuc_2_num(intron_seq[l_intron - 2]) == nuc_2_num(cur_motif_str[0]) &&
					nuc_2_num(intron_seq[l_intron - 1]) == nuc_2_num(cur_motif_str[1]))
				{
					splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
				}
			}
		} // Right side of the intron check.

		// Check the gains.
		if (1 == 0)
		{
			char cur_motif_str[10];
			strcpy(cur_motif_str, splice_donor_motif);
			reverse_complement_seq(cur_motif_str);

			// Check gains.
			if (relative_var_posn_in_intron + 1 < l_intron &&
				toupper(intron_seq[relative_var_posn_in_intron]) != cur_motif_str[0] &&
				toupper(intron_seq[relative_var_posn_in_intron + 1]) == cur_motif_str[1] &&
				toupper(var_alternate_allele) == cur_motif_str[0])
			{
				splice_effect_type = EFFECT_INTRON_SPLICE_GAIN;
			}

			if (relative_var_posn_in_intron > 1 &&
				toupper(intron_seq[relative_var_posn_in_intron - 1]) == cur_motif_str[0] &&
				toupper(intron_seq[relative_var_posn_in_intron]) != cur_motif_str[1] &&
				toupper(var_alternate_allele) == cur_motif_str[1])
			{
				splice_effect_type = EFFECT_INTRON_SPLICE_GAIN;
			}
		}
	} // snv effect check.
	else if(cur_var_info->var_type == VAR_TYPE_INSERTION)
	{
		// Check splice motif loss: The insertion's start position can be -1, which points to the fact that the insertion is between exon and intron boundary.
		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "@ Insertion %s:%d-%d (%s); Intron: %d-%d\n", var_reg->chrom, var_reg->start, var_reg->end, 
					cur_var_info->alt_allele, 
					intronic_reg->start, intronic_reg->end);
		}

		if(relative_var_posn_in_intron <= 0)
		{
			splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
		}

		if(relative_var_posn_in_intron >= l_intron-2)
		{
			splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
		}

		// TODO::Check splicing motif gain, too.
	} // insertion effect check.
	else if (cur_var_info->var_type == VAR_TYPE_DELETION)
	{
		int l_deletion = t_string::string_length(cur_var_info->ref_allele);

		if (__DUMP_VARIATION_TOOLS_MSGS__)
		{
			fprintf(stderr, "@ Deletion %s:%d-%d (%s); Intron: %d-%d\n", var_reg->chrom, var_reg->start, var_reg->end,
				cur_var_info->ref_allele,
				intronic_reg->start, intronic_reg->end);
		}

		// Check splice motif loss: 
		if ((relative_var_posn_in_intron + l_deletion) <= 1)
		{
			splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
		}

		if ((relative_var_posn_in_intron + l_deletion) >= l_intron - 2)
		{
			splice_effect_type = EFFECT_INTRON_SPLICE_LOSS;
		}

		// TODO::Check splicing motif gain, too.
	}
}

char get_nuc_per_transcript_i(t_annot_region* transcript, int nuc_trans_i)
{
	vector<t_annot_region*>* cds_regs = ((t_transcript_annotation_info*)(transcript->annotation_info))->CDSs;
	sort(cds_regs->begin(), cds_regs->end(), sort_regions);

	for(int i_cds = 0; i_cds < cds_regs->size(); i_cds++)
	{
		t_CDS_annotation_info* cur_cds_info = (t_CDS_annotation_info*)(cds_regs->at(i_cds)->annotation_info);

		int l_cds = cds_regs->at(i_cds)->end - cds_regs->at(i_cds)->start + 1;

		if(cur_cds_info->l_transcript_so_far <= nuc_trans_i &&
			cur_cds_info->l_transcript_so_far + l_cds > nuc_trans_i)
		{
			int nuc_i_per_cur_cds = nuc_trans_i - cur_cds_info->l_transcript_so_far;
			return(cur_cds_info->sequence[nuc_i_per_cur_cds]);
		}
	} // i_cds loop.

	fprintf(stderr, "Could not find nucleotide:\n");
	exit(0);
	return(-1);
}

char** get_variant_effect_name_strings_list()
{
	char** effect_name_strings = new char*[N_VAR_EFFECT_TYPES + 2];
	effect_name_strings[EFFECT_INTERGENIC] = t_string::copy_me_str("EFFECT_INTERGENIC");
	effect_name_strings[EFFECT_CDS_SYNONYMOUS] = t_string::copy_me_str("EFFECT_CDS_SYNONYMOUS");
	effect_name_strings[EFFECT_CDS_NONSYNONYMOUS] = t_string::copy_me_str("EFFECT_CDS_NONSYNONYMOUS");
	effect_name_strings[EFFECT_CDS_EARLY_STOP] = t_string::copy_me_str("EFFECT_CDS_EARLY_STOP");
	effect_name_strings[EFFECT_CDS_STOP_DISRUPTION] = t_string::copy_me_str("EFFECT_CDS_STOP_DISRUPTION");
	effect_name_strings[EFFECT_CDS_FRAMESHIFT_INDEL] = t_string::copy_me_str("EFFECT_CDS_FRAMESHIFT_INDEL");
	effect_name_strings[EFFECT_CDS_INFRAME_INDEL] = t_string::copy_me_str("EFFECT_CDS_INFRAME_INDEL");
	effect_name_strings[EFFECT_CDS_POST_INFRAME_INSERTION] = t_string::copy_me_str("EFFECT_CDS_POST_INFRAME_INSERTION");
	effect_name_strings[EFFECT_CDS_CNV] = t_string::copy_me_str("EFFECT_CDS_CNV");
	effect_name_strings[EFFECT_TF_MOTIF_BREAK] = t_string::copy_me_str("EFFECT_TF_MOTIF_BREAK");
	effect_name_strings[EFFECT_TF_MOTIF_GAIN] = t_string::copy_me_str("EFFECT_TF_MOTIF_GAIN");
	effect_name_strings[EFFECT_5P_UTR_OVERLAP] = t_string::copy_me_str("EFFECT_5P_UTR_OVERLAP");
	effect_name_strings[EFFECT_3P_UTR_OVERLAP] = t_string::copy_me_str("EFFECT_3P_UTR_OVERLAP");
	effect_name_strings[EFFECT_PROMOTOR_OVERLAP] = t_string::copy_me_str("EFFECT_PROMOTOR_OVERLAP");
	effect_name_strings[EFFECT_INTRON_SPLICE_GAIN] = t_string::copy_me_str("EFFECT_INTRON_SPLICE_GAIN");
	effect_name_strings[EFFECT_INTRON_SPLICE_LOSS] = t_string::copy_me_str("EFFECT_INTRON_SPLICE_LOSS");
	effect_name_strings[EFFECT_INTRON_OVERLAP] = t_string::copy_me_str("EFFECT_INTRON_OVERLAP");

	return(effect_name_strings);
}

/*
TODO::Add intronic regions and splice donor creating variant?

Use the elements loaded with load_elements_set_annotation_information_per_GENCODE_GFF(...).
*/
void add_variant_effect_information(vector<t_annot_region*>* candidate_var_regs, 
									vector<t_annot_region*>* gencode_transcript_regs,
									vector<t_annot_region*>* gencode_cds_regs,
									vector<t_annot_region*>* gencode_5p_UTR_regs,
									vector<t_annot_region*>* gencode_3p_UTR_regs,
									vector<t_annot_region*>* gencode_promoter_regs,
									vector<t_annot_region*>* gencode_intronic_regs,
									vector<t_annot_region*>* encode_motif_regs)
{
	fprintf(stderr, "Adding variant effect information, using elements loaded with load_elements_set_annotation_information_per_GENCODE_GFF(...)\n");

	// Laod and initialize the gencode transcript regions with sequence information.
	fprintf(stderr, "Loaded %d transcripts and %d CDSs, assigning variant effects.\n", gencode_transcript_regs->size(), gencode_cds_regs->size());

	// Set all the variant effects to NULL.
	for(int i_var = 0; i_var < candidate_var_regs->size(); i_var++)
	{
		t_variant_info* cur_var_info = (t_variant_info*)(candidate_var_regs->at(i_var)->annotation_info);
		if(cur_var_info->variant_effects == NULL)
		{
			cur_var_info->variant_effects = new vector<t_variant_effect_info*>();
		}
	} // i_var loop.

	/******************************************************************************************
	TODO::Include promoters, and motifs in the following.
	*******************************************************************************************/
	// Intersect variants with CDSs and assign effect info.
	vector<t_annot_region*>* cds_intersects = intersect_annot_regions(gencode_cds_regs, candidate_var_regs, true);
	for(int i_int = 0; i_int < cds_intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(cds_intersects->at(i_int)->data);
		t_annot_region* cds_reg = cur_int_info->src_reg;
		t_annot_region* var_reg = cur_int_info->dest_reg;

		t_CDS_annotation_info* cur_cds_info = (t_CDS_annotation_info*)(cds_reg->annotation_info);
		t_variant_info* cur_var_info = (t_variant_info*)(var_reg->annotation_info);

		// If there is no variant information, then skip.
		if(cur_var_info == NULL)
		{
			fprintf(stderr, "Variant annotation requires annotation_info field to contain a t_variant_info entry: %d. entry has null\n", i_int);
			exit(0);
		}

		// First, make sure the transcript of this CDS has valid 3mer frame.
		t_transcript_annotation_info* cur_cds_transcript_info = ((t_transcript_annotation_info*)(cur_cds_info->transcript->annotation_info));
		if(cur_cds_transcript_info->valid_3mer_frame)
		{
			if(cur_var_info->var_type == VAR_TYPE_SNV)
			{
				// Get the CDS regions and sort them.
				vector<t_annot_region*>* cur_trans_cds_regs = ((t_transcript_annotation_info*)(cur_cds_info->transcript->annotation_info))->CDSs;
				sort(cur_trans_cds_regs->begin(), cur_trans_cds_regs->end(), sort_regions);

				// Check if the mutation overlaps with the first or the last codon.
				char reference_codon[10];
				char mutated_codon[10];
				int affected_nuc_i_in_transcript_coords = (var_reg->start - cds_reg->start) + cur_cds_info->l_transcript_so_far;
				int affected_codon_start_in_transcript_coords = affected_nuc_i_in_transcript_coords - (affected_nuc_i_in_transcript_coords % 3);

if(__DUMP_VARIATION_TOOLS_MSGS__)
{
				fprintf(stderr, "%s:%d-%d::%c<->%c @ %s::%s:%d-%d (%d, %d)\n", 
						var_reg->chrom, var_reg->start, var_reg->end, 
						cur_var_info->ref_allele[0], cur_var_info->alt_allele[0], 
						cur_cds_info->transcript->name, cds_reg->chrom, cds_reg->start, cds_reg->end, 
						((t_transcript_annotation_info*)(cur_cds_info->transcript->annotation_info))->l_transcript,
						((t_transcript_annotation_info*)(cur_cds_info->transcript->annotation_info))->valid_3mer_frame);

				fprintf(stderr, "affected_nuc_i_in_transcript_coords: %d, affected_codon_start_in_transcript_coords: %d\n", 
						affected_nuc_i_in_transcript_coords, affected_codon_start_in_transcript_coords);
}

				for(int nuc_i = affected_codon_start_in_transcript_coords; nuc_i < affected_codon_start_in_transcript_coords+3; nuc_i++)
				{
if(__DUMP_VARIATION_TOOLS_MSGS__)
{
					fprintf(stderr, "Getting nuc %d: ", nuc_i);
}

					char cur_trans_nuc = get_nuc_per_transcript_i(cur_cds_info->transcript, nuc_i);

					reference_codon[nuc_i - affected_codon_start_in_transcript_coords] = cur_trans_nuc;
					mutated_codon[nuc_i - affected_codon_start_in_transcript_coords] = cur_trans_nuc;

					if(nuc_i == affected_nuc_i_in_transcript_coords)
					{
						mutated_codon[nuc_i - affected_codon_start_in_transcript_coords] = cur_var_info->alt_allele[0];
					}

if(__DUMP_VARIATION_TOOLS_MSGS__)
{
					fprintf(stderr, "%c\n", cur_trans_nuc);
}
				} // nuc_i loop.

				mutated_codon[3] = 0;
				reference_codon[3] = 0;

				// If the CDS is on the negative strand, reverse complement it.
				if(cds_reg->strand == '-')
				{
					reverse_complement_seq(mutated_codon);
					reverse_complement_seq(reference_codon);
				}

				char* mutated_AA = get_AA_code_per_codon(mutated_codon);
				char* reference_AA = get_AA_code_per_codon(reference_codon);

				t_variant_effect_info* var_effect = new t_variant_effect_info();
				var_effect->effected_element_region = cds_reg;
				var_effect->effected_element_type = EFFECTED_ELEMENT_CDS;
				var_effect->effect_type = EFFECT_INTERGENIC;

				t_CDS_effect_info* CDS_var_effect = new t_CDS_effect_info();
				CDS_var_effect->alt_AA = t_string::copy_me_str(mutated_AA);
				CDS_var_effect->ref_AA = t_string::copy_me_str(reference_AA);
				var_effect->element_specific_effect_info = CDS_var_effect;

				// Add the current vairant effect to the list of variant effects for this variant.
				cur_var_info->variant_effects->push_back(var_effect);

				if(t_string::compare_strings(mutated_AA, "STOP"))
				{
					var_effect->effect_type = EFFECT_CDS_EARLY_STOP;
				}
				else if(t_string::compare_strings(reference_AA, "STOP"))
				{
					var_effect->effect_type = EFFECT_CDS_STOP_DISRUPTION;
				}
				else if(!t_string::compare_strings(reference_AA, mutated_AA))
				{
					var_effect->effect_type = EFFECT_CDS_NONSYNONYMOUS;
				}
				else
				{
					var_effect->effect_type = EFFECT_CDS_SYNONYMOUS;
				}
			} // snv effect type setting.
			else if (cur_var_info->var_type == VAR_TYPE_INSERTION)
			{
				// Insertion coordinates are described by the left and right coordinates of the breakpoint insertion site.
				// If the start is exactly at the end of a frame and end is the start of the next frame, this could be an inframe insertion if the length is multiple of 3.
				// Below aims to align the end of the variant to the start of the CDS's start position.
				int affected_nuc_i_in_transcript_coords = (var_reg->end - cds_reg->start) + cur_cds_info->l_transcript_so_far;
				int affected_codon_start_in_transcript_coords = affected_nuc_i_in_transcript_coords - (affected_nuc_i_in_transcript_coords % 3);

				int l_insert = strlen(cur_var_info->alt_allele);
				t_variant_effect_info* var_effect = new t_variant_effect_info();
				var_effect->effected_element_region = cds_reg;
				var_effect->effected_element_type = EFFECTED_ELEMENT_CDS;
				var_effect->effect_type = EFFECT_INTERGENIC; // Initialize the effect.

				t_CDS_effect_info* CDS_var_effect = new t_CDS_effect_info();
				CDS_var_effect->alt_AA = t_string::copy_me_str("INS_AA");
				CDS_var_effect->ref_AA = t_string::copy_me_str(".");
				var_effect->element_specific_effect_info = CDS_var_effect;

				// Add the current vairant effect to the list of variant effects for this variant.
				cur_var_info->variant_effects->push_back(var_effect);

				// The indel type is either inframe of frameshift indel.
				if (l_insert % 3 == 0 &&
					affected_nuc_i_in_transcript_coords % 3 == 0)
				{
					// This is an in-frame deletion.
					var_effect->effect_type = EFFECT_CDS_INFRAME_INDEL;
				}
				else if (affected_codon_start_in_transcript_coords + l_insert % 3 == 0)
				{
					var_effect->effect_type = EFFECT_CDS_POST_INFRAME_INSERTION;
				}
				else
				{
					var_effect->effect_type = EFFECT_CDS_FRAMESHIFT_INDEL;
				}
			} // insertion effect type setting.
			else if (cur_var_info->var_type == VAR_TYPE_DELETION)
			{
				// Deletion coordinates are described by the left and right coordinates of the deleted nucleotides in genomic coordinates.
				// Below aims to align the end of the variant to the start of the CDS's start position.
				int affected_nuc_i_in_transcript_coords = (var_reg->start - cds_reg->start) + cur_cds_info->l_transcript_so_far;
				int affected_codon_start_in_transcript_coords = affected_nuc_i_in_transcript_coords - (affected_nuc_i_in_transcript_coords % 3);

				int l_deletion = strlen(cur_var_info->ref_allele);
				t_variant_effect_info* var_effect = new t_variant_effect_info();
				var_effect->effected_element_region = cds_reg;
				var_effect->effected_element_type = EFFECTED_ELEMENT_CDS;
				var_effect->effect_type = EFFECT_INTERGENIC; // Initialize the effect.

				t_CDS_effect_info* CDS_var_effect = new t_CDS_effect_info();
				CDS_var_effect->alt_AA = t_string::copy_me_str(".");
				CDS_var_effect->ref_AA = t_string::copy_me_str("DEL_AA");
				var_effect->element_specific_effect_info = CDS_var_effect;

				// Add the current vairant effect to the list of variant effects for this variant.
				cur_var_info->variant_effects->push_back(var_effect);

				// The indel type is either inframe of frameshift indel.
				if (l_deletion % 3 == 0 &&
					affected_nuc_i_in_transcript_coords % 3 == 0)
				{
					// This is an in-frame deletion.
					var_effect->effect_type = EFFECT_CDS_INFRAME_INDEL;
				}
				else
				{
					var_effect->effect_type = EFFECT_CDS_FRAMESHIFT_INDEL;
				}
			} // deletion effect type setting.
			else if(cur_var_info->var_type == VAR_TYPE_CNV)
			{
				// Set the CNV effect type.
				int affected_nuc_i_in_transcript_coords = (var_reg->start - cds_reg->start) + cur_cds_info->l_transcript_so_far;
				int affected_codon_start_in_transcript_coords = affected_nuc_i_in_transcript_coords - (affected_nuc_i_in_transcript_coords % 3);

				int l_indel = var_reg->end - var_reg->start + 1;
				t_variant_effect_info* var_effect = new t_variant_effect_info();
				var_effect->effected_element_region = cds_reg;
				var_effect->effected_element_type = EFFECTED_ELEMENT_CDS;

				t_CDS_effect_info* CDS_var_effect = new t_CDS_effect_info();
				CDS_var_effect->alt_AA = t_string::copy_me_str("CNV_SEQ");
				CDS_var_effect->ref_AA = t_string::copy_me_str("CNV_SEQ");
				var_effect->element_specific_effect_info = CDS_var_effect;

				var_effect->effect_type = EFFECT_CDS_CNV;

				// Add the current vairant effect to the list of variant effects for this variant.
				cur_var_info->variant_effects->push_back(var_effect);
			}
		} // valid transcript frame check.

		delete cur_int_info;
	} // i_int loop.

	// Intersect the SNVs with the intronic regions; find the splice donor/acceptor gains.
	vector<t_annot_region*>* intronic_intersects = intersect_annot_regions(gencode_intronic_regs, candidate_var_regs, true);
	for(int i_int = 0; i_int < intronic_intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intronic_intersects->at(i_int)->data);
		t_annot_region* intronic_reg = cur_int_info->src_reg;
		t_annot_region* var_reg = cur_int_info->dest_reg;

		// Check if there is a splice acceptor/donor site addition or deletion.
		int splice_effect_type = -1;
		get_SNV_splice_effect(var_reg,
								intronic_reg,
								splice_effect_type);

		// If a splicing effect is assigned, add to the list.
		if(splice_effect_type != -1)
		{
			// If this is an SNV, check for acceptor site addition.
			t_variant_effect_info* var_effect = new t_variant_effect_info();

			var_effect->effected_element_region = intronic_reg;
			var_effect->effected_element_type = EFFECTED_ELEMENT_INTRON;

			var_effect->effect_type = splice_effect_type;

			t_intron_effect_info* intron_effect_info = new t_intron_effect_info();
			var_effect->element_specific_effect_info = intron_effect_info;

			// Add this variant effect.
			((t_variant_info*)(var_reg->annotation_info))->variant_effects->push_back(var_effect);
		}
		else
		{
			// If this is an SNV, check for acceptor site addition.
			t_variant_effect_info* var_effect = new t_variant_effect_info();

			var_effect->effected_element_region = intronic_reg;
			var_effect->effected_element_type = EFFECTED_ELEMENT_INTRON;

			var_effect->effect_type = EFFECT_INTRON_OVERLAP;

			t_intron_effect_info* intron_effect_info = new t_intron_effect_info();
			var_effect->element_specific_effect_info = intron_effect_info;

			// Add this variant effect.
			((t_variant_info*)(var_reg->annotation_info))->variant_effects->push_back(var_effect);
		}
	} // i_int loop.
	
	// 5' UTRs.
	vector<t_annot_region*>* UTR_5p_intersects = intersect_annot_regions(gencode_5p_UTR_regs, candidate_var_regs, true);
	fprintf(stderr, "%d 5' UTR intersects.\n", UTR_5p_intersects->size());
	for (int i_int = 0; i_int < UTR_5p_intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(UTR_5p_intersects->at(i_int)->data);
		t_annot_region* utr_reg = cur_int_info->src_reg;
		t_annot_region* var_reg = cur_int_info->dest_reg;

		// If this is an SNV, check for acceptor site addition.
		t_variant_effect_info* var_effect = new t_variant_effect_info();

		var_effect->effected_element_region = utr_reg;
		var_effect->effected_element_type = EFFECTED_ELEMENT_UTR;

		var_effect->effect_type = EFFECT_5P_UTR_OVERLAP;

		t_intron_effect_info* intron_effect_info = new t_intron_effect_info();
		var_effect->element_specific_effect_info = intron_effect_info;

		// Add this variant effect.
		((t_variant_info*)(var_reg->annotation_info))->variant_effects->push_back(var_effect);
	} // i_int loop.

	// 3' UTRs.
	vector<t_annot_region*>* UTR_3p_intersects = intersect_annot_regions(gencode_3p_UTR_regs, candidate_var_regs, true);
	for (int i_int = 0; i_int < UTR_3p_intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(UTR_3p_intersects->at(i_int)->data);
		t_annot_region* utr_reg = cur_int_info->src_reg;
		t_annot_region* var_reg = cur_int_info->dest_reg;

		// If this is an SNV, check for acceptor site addition.
		t_variant_effect_info* var_effect = new t_variant_effect_info();

		var_effect->effected_element_region = utr_reg;
		var_effect->effected_element_type = EFFECTED_ELEMENT_UTR;

		var_effect->effect_type = EFFECT_3P_UTR_OVERLAP;

		t_intron_effect_info* intron_effect_info = new t_intron_effect_info();
		var_effect->element_specific_effect_info = intron_effect_info;

		// Add this variant effect.
		((t_variant_info*)(var_reg->annotation_info))->variant_effects->push_back(var_effect);
	} // i_int loop.


	vector<t_annot_region*>* promotor_intersects = intersect_annot_regions(gencode_promoter_regs, candidate_var_regs, true);
	for (int i_int = 0; i_int < promotor_intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(promotor_intersects->at(i_int)->data);
		t_annot_region* promotor_reg = cur_int_info->src_reg;
		t_annot_region* var_reg = cur_int_info->dest_reg;

		// If this is an SNV, check for acceptor site addition.
		t_variant_effect_info* var_effect = new t_variant_effect_info();

		var_effect->effected_element_region = promotor_reg;
		var_effect->effected_element_type = EFFECTED_ELEMENT_PROMOTOR;

		var_effect->effect_type = EFFECT_PROMOTOR_OVERLAP;

		t_intron_effect_info* intron_effect_info = new t_intron_effect_info();
		var_effect->element_specific_effect_info = intron_effect_info;

		// Add this variant effect.
		((t_variant_info*)(var_reg->annotation_info))->variant_effects->push_back(var_effect);
	} // i_int loop.

}

vector<t_annot_region*>* load_VCF_regions(char* vcf_fp, bool parse_genotypes)
{
	vector<t_annot_region*>* vcf_regions = new vector<t_annot_region*>();

	FILE* f_vcf = open_f(vcf_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_vcf);
		if(cur_line == NULL)
		{
			break;
		}

		if(cur_line[0] == '#')
		{
			delete [] cur_line;
			continue;
		}

		// Parse the line and allocate the region.
		//t_string* cur_line_str = new t_string(cur_line);
		//t_string_tokens* cur_line_tokens = cur_line_str->tokenize_by_chars("\t");
		t_string_tokens* cur_line_tokens = t_string::tokenize_by_chars(cur_line, "\t");
	
		if (cur_line_tokens->size() < 8)
		{
			fprintf(stderr, "Could not find 8 mandatory fields @ %s(%d): %s\n", __FILE__, __LINE__, cur_line);
			delete[] cur_line;
			t_string::clean_tokens(cur_line_tokens);
			continue;
		}

		// 1       2337697 rs55744642;-/CGACA,venter,watson        GG      GCGACAG,GGGACAG 100     PASS    CA=0;DB;NF=16;NR=14;HP=3;DP=179 GT:GQ   2|1:100
		t_annot_region* new_region = get_empty_region();

		// 
		new_region->chrom = t_string::copy_me_str(cur_line_tokens->at(0)->str());
		new_region->start = translate_coord(atoi(cur_line_tokens->at(1)->str()), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
		new_region->end = new_region->start;

		// Copy name.
		new_region->name = t_string::copy_me_str(cur_line_tokens->at(2)->str());

		// The strand is positive.
		new_region->strand = '+';
	
		// Following are the fixed fields.
		t_vcf_info* new_vcf_info = new t_vcf_info();
		new_vcf_info->ref_pos = atoi(cur_line_tokens->at(1)->str());
		new_vcf_info->ref_allele_str = t_string::copy_me_str(cur_line_tokens->at(3)->str());
		new_vcf_info->alt_allele_str = t_string::copy_me_str(cur_line_tokens->at(4)->str());
		new_vcf_info->info_str = t_string::copy_me_str(cur_line_tokens->at(7)->str());

		// Set the vcf info to data.
		new_region->data = new_vcf_info;

		// Parse the genotypes?
		if(parse_genotypes)
		{
			fprintf(stderr, "Cannot parse the genotypes, yet.\n");
			exit(0);
		}

		vcf_regions->push_back(new_region);

		// Clean the line data.
		//delete(cur_line_str);
		t_string::clean_tokens(cur_line_tokens);
		delete [] cur_line;
	} // file reading loop.

	close_f(f_vcf, vcf_fp);

	return(vcf_regions);
}

t_VCF::t_VCF(char* vcf_fp)
{
	FILE* f_vcf = open_f(vcf_fp, "r");

	this->entries = new vector<t_vcf_entry*>();
	this->entries_per_chromosome = new vector<vector<t_vcf_entry*>*>();	

	// Allocate the chromsome id's.
	this->chrom_ids = new vector<char*>();

	//char* cur_line = new char[1000000];
	while(1)
	{
		//fprintf(stderr, "Reading a new line.\n");
		char* cur_line = getline(f_vcf);
		//fgets(cur_line, 1000000, f_vcf);
		
		//fprintf(stderr, "Read a new line: %d\n", strlen(cur_line));

		if(cur_line == NULL)
		{
			break;
		}

		// Skip comments.
		if(cur_line[0] == '#')
		{
			delete [] cur_line;
		}
		else
		{
			//fprintf(stderr, "Tokenizing a new line.\n");
			// #CHROM  POS     ID	   REF     ALT	   QUAL    FILTER  INFO							FORMAT  NA12878
			// 1       233713  .       CAAGTG  CG      100     PASS    DP=179;HP=2;CA=M;NR=6;NF=5	GT:GQ	0|1:58
			// 1       2337697 rs55744642;-/CGACA,venter,watson        GG      GCGACAG,GGGACAG 100     PASS    CA=0;DB;NF=16;NR=14;HP=3;DP=179 GT:GQ   2|1:100
			//t_string* cur_line_str = new t_string(cur_line);
			//t_string_tokens* line_str_tokens = cur_line_str->tokenize_by_chars(" \t");

			//fprintf(stderr, "%d: %d tokens.\n", this->entries->size(), line_str_tokens->size());

			//fprintf(stderr, "Copying things.\n");
			//int chr_id = atoi(line_str_tokens->at(0)->str());
			//int pos_ref = atoi(line_str_tokens->at(1)->str());
			//char* ref_allele = new char[line_str_tokens->at(3)->length() + 2];
			//char* alt_allele = new char[line_str_tokens->at(4)->length() + 2];
			//strcpy(ref_allele, line_str_tokens->at(3)->str());
			//strcpy(alt_allele, line_str_tokens->at(4)->str());

			t_string* cur_line_str = new t_string(cur_line);
			t_string_tokens* cur_line_tokens = cur_line_str->tokenize_by_chars("\n\t");

			//char* _cur_ref_allele_str = new char[strlen(cur_line) + 2];
			//char* _cur_alt_allele_str = new char[strlen(cur_line) + 2];
			//char* _cur_info = new char[strlen(cur_line) + 2];
			//char* _genotype_format = new char[strlen(cur_line) + 2];

			//char cur_chrom[100];
			//char cur_entry_id[1000];
			//int cur_ref_pos;
			//char cur_qual[100];
			//char cur_filter[100];
			//
			//if(sscanf(cur_line, "%s %d %s %s %s %s %s %s", 
			//	cur_chrom, 
			//	&cur_ref_pos, 
			//	cur_entry_id,
			//	_cur_ref_allele_str, 
			//	_cur_alt_allele_str, 
			//	cur_qual,
			//	cur_filter,
			//	_cur_info) != 8)
			//{
			//	printf("Could not parse VCF line:\n%s\n", cur_line);
			//	exit(0);
			//}

			if(this->entries->size() % 1000 == 0)
			{
				fprintf(stderr, "%d            \r", this->entries->size());
			}

			t_vcf_entry* new_entry = new t_vcf_entry();

			// Do chromosome id check.
			char chrom[100];
			if(!cur_line_tokens->at(0)->starts_with("chr"))
			{
				sprintf(chrom, "chr%s", cur_line_tokens->at(0)->str());				
			}
			else
			{
				strcpy(chrom, cur_line_tokens->at(0)->str());
			}

			new_entry->chrom = new char[strlen(chrom) + 2];
			strcpy(new_entry->chrom, chrom);

			new_entry->ref_pos = atoi(cur_line_tokens->at(1)->str());

			new_entry->id = new char[strlen(cur_line_tokens->at(2)->str()) + 2];
			strcpy(new_entry->id, cur_line_tokens->at(2)->str());

			//char* _cur_ref_allele_str = new char[strlen(cur_line_tokens->at(3)->str()) + 2];
			//strcpy(_cur_ref_allele_str, cur_line_tokens->at(3)->str());

			//char* _cur_alt_allele_str = new char[strlen(cur_line_tokens->at(4)->str()) + 2];
			//strcpy(_cur_alt_allele_str, cur_line_tokens->at(4)->str());

			//char* qual_str = new char[strlen(cur_line_tokens->at(5)->str()) + 2];
			//strcpy(qual_str, cur_line_tokens->at(5)->str());

			// Skip filter string.
			new_entry->info_str = new char[strlen(cur_line_tokens->at(7)->str()) + 2];
			strcpy(new_entry->info_str, cur_line_tokens->at(7)->str());

			// Parse the reference alleles and alternate alleles.
			new_entry->ref_allele = new t_string(cur_line_tokens->at(3));
			new_entry->alt_alleles = cur_line_tokens->at(4)->tokenize_by_chars(",");

/*
			// Dump the parsed information.
			fprintf(stderr, "%s: %s(%d)\n\
Ref: %s\n\
Alt: ", new_entry->id, new_entry->chrom, new_entry->ref_pos, new_entry->ref_allele->str());

			for(int i_alt = 0; i_alt < new_entry->alt_alleles->size(); i_alt++)
			{
				fprintf(stderr, "%s, ", new_entry->alt_alleles->at(i_alt)->str());
			} // i_alt loop.
			fprintf(stderr, "\nInfo: %s\n", new_entry->info_str);

			// If there is genotype data, store it.
			if(cur_line_tokens->size() > 8)
			{
				// Allocate the genotypes.
				new_entry->genotypes = new vector<t_genotype*>();

				// Store format.
				new_entry->genotype_format = new char[strlen(cur_line_tokens->at(8)->str()) + 2];
				strcpy(new_entry->genotype_format, cur_line_tokens->at(8)->str());

				fprintf(stderr, "Genotype format: %s\n", new_entry->genotype_format);

				// Assuming that the genotype is at the 0th entry.
				int i_gt_entry = 0;

				// Parse the genotype information if it is there.
				for(int i_ent = 9; i_ent < cur_line_tokens->size(); i_ent++)
				{
					t_string* gt_str = cur_line_tokens->at(i_ent);
					t_string_tokens* gt_str_tokens = gt_str->tokenize_by_chars(":");
					char* cur_genotype = gt_str_tokens->at(i_gt_entry)->str();

					// Process 3 characters in the genotype string. Note that the genotype can be called with for haplotype organism.
					if(strlen(cur_genotype) != 3)
					{
						fprintf(stderr, "The genotype string is not recognizable.\n", cur_genotype);
						exit(0);
					}

					t_genotype* new_gt = new t_genotype();
					new_gt->pat_allele = cur_genotype[0];
					new_gt->mat_allele = cur_genotype[2];
					
					// Check phasing information.
					if(cur_genotype[1] == '|')
					{
						new_gt->phased = true;
						//fprintf(stderr, "Phased(%c, %c)\n", new_gt->pat_allele, new_gt->mat_allele);
					}
					else if(cur_genotype[1] == '/')
					{
						new_gt->phased = false;
						//fprintf(stderr, "Unphased(%c, %c)\n", new_gt->pat_allele, new_gt->mat_allele);
					}

					// Add the new genotype to the list of genotypes.
					new_entry->genotypes->push_back(new_gt);

					// Clean tokens.
					t_string::clean_tokens(gt_str_tokens);
				} // i_ent loop.
			} // Check if the genotype data is includeded.
			else
			{
				new_entry->genotypes = NULL;
			}
*/

			t_string::clean_tokens(cur_line_tokens);
			delete(cur_line_str);

			// Add this entry.
			this->entries->push_back(new_entry);

			int i_cur_chr = t_string::get_i_str(this->chrom_ids, chrom);

			if(i_cur_chr == this->chrom_ids->size())
			{
				// Add the new chromosome and the new list of vcf entries.
				this->chrom_ids->push_back(t_string::copy_me_str(chrom));
				this->entries_per_chromosome->push_back(new vector<t_vcf_entry*>());

				i_cur_chr = t_string::get_i_str(this->chrom_ids, chrom);

				if(i_cur_chr == this->chrom_ids->size())
				{
					fprintf(stderr, "Chromosome id addition failed for %s\n", chrom);
					exit(0);
				}
			}

			// At this point, i_cur_chr points to a valid list of vcf entries.
			this->entries_per_chromosome->at(i_cur_chr)->push_back(new_entry);

			// Free memory.
			delete [] cur_line;
		} // comment check.
	} // file reading loop.

	printf("Loaded %d VCF entries.\n", this->entries->size());

	fclose(f_vcf);
}

void t_VCF::dump_VCF_BED(char* bed_fp)
{
	FILE* f_bed = open_f(bed_fp, "w");
	for(int i_ent = 0; i_ent < this->entries->size(); i_ent++)
	{
		int start = this->entries->at(i_ent)->ref_pos;
		int var_length = strlen(this->entries->at(i_ent)->ref_allele->str());
		
		//if(strlen(this->entries->at(i_ent)->ref_allele) > strlen(this->entries->at(i_ent)->alt_allele))
		//{
		//	var_length = strlen(this->entries->at(i_ent)->ref_allele) - strlen(this->entries->at(i_ent)->alt_allele) + 1;
		//}
		//else
		//{
		//	var_length = strlen(this->entries->at(i_ent)->alt_allele) - strlen(this->entries->at(i_ent)->ref_allele) + 1;
		//}

		int end = start + var_length;

		//fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\t%s\t%s\n", this->entries->at(i_ent)->i_chrom,
		//	start, end, this->entries->at(i_ent)->ref_allele, this->entries->at(i_ent)->alt_allele);

		fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\n", this->entries->at(i_ent)->chrom,
			start, end);
	} // i_ent loop.
	fclose(f_bed);
}

void t_VCF::dump_het_variants_in_ref_coords(char* bed_fp)
{
	FILE* f_bed = open_f(bed_fp, "w");
	for(int i_ent = 0; i_ent < this->entries->size(); i_ent++)
	{
		// Check if the variant is alternate in one of the parents and it is het, i.e., different in the parents.
		if(this->entries->at(i_ent)->genotypes != NULL &&
			(this->entries->at(i_ent)->genotypes->at(0)->pat_allele != '0' ||
			this->entries->at(i_ent)->genotypes->at(0)->mat_allele != '0') &&
			this->entries->at(i_ent)->genotypes->at(0)->mat_allele != this->entries->at(i_ent)->genotypes->at(0)->pat_allele)
		{
			int start = this->entries->at(i_ent)->ref_pos;
			int var_length = strlen(this->entries->at(i_ent)->ref_allele->str());
		
			//if(strlen(this->entries->at(i_ent)->ref_allele) > strlen(this->entries->at(i_ent)->alt_allele))
			//{
			//	var_length = strlen(this->entries->at(i_ent)->ref_allele) - strlen(this->entries->at(i_ent)->alt_allele) + 1;
			//}
			//else
			//{
			//	var_length = strlen(this->entries->at(i_ent)->alt_allele) - strlen(this->entries->at(i_ent)->ref_allele) + 1;
			//}

			int end = start + var_length;

			//fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\t%s\t%s\n", this->entries->at(i_ent)->i_chrom,
			//	start, end, this->entries->at(i_ent)->ref_allele, this->entries->at(i_ent)->alt_allele);

			fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\n", this->entries->at(i_ent)->chrom,
				start, end);
		}
	} // i_ent loop.
	fclose(f_bed);
}

void t_VCF::dump_paternal_alternate_sites_in_ref_coords(char* bed_fp)
{
	FILE* f_bed = open_f(bed_fp, "w");
	for(int i_ent = 0; i_ent < this->entries->size(); i_ent++)
	{
		// Check if the variant is in paternal haplotype.
		if(this->entries->at(i_ent)->genotypes != NULL &&
			this->entries->at(i_ent)->genotypes->at(0)->pat_allele != '0')
		{
			int start = this->entries->at(i_ent)->ref_pos;
			int var_length = strlen(this->entries->at(i_ent)->ref_allele->str());
		
			//if(strlen(this->entries->at(i_ent)->ref_allele) > strlen(this->entries->at(i_ent)->alt_allele))
			//{
			//	var_length = strlen(this->entries->at(i_ent)->ref_allele) - strlen(this->entries->at(i_ent)->alt_allele) + 1;
			//}
			//else
			//{
			//	var_length = strlen(this->entries->at(i_ent)->alt_allele) - strlen(this->entries->at(i_ent)->ref_allele) + 1;
			//}

			int end = start + var_length;

			//fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\t%s\t%s\n", this->entries->at(i_ent)->i_chrom,
			//	start, end, this->entries->at(i_ent)->ref_allele, this->entries->at(i_ent)->alt_allele);

			fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\n", this->entries->at(i_ent)->chrom,
				start, end);
		}
	} // i_ent loop.
	fclose(f_bed);
}

void t_VCF::dump_maternal_alternate_sites_in_ref_coords(char* bed_fp)
{
	FILE* f_bed = open_f(bed_fp, "w");
	for(int i_ent = 0; i_ent < this->entries->size(); i_ent++)
	{
		// Check if the variant is in paternal haplotype.
		if(this->entries->at(i_ent)->genotypes != NULL &&
			this->entries->at(i_ent)->genotypes->at(0)->mat_allele != '0')
		{
			int start = this->entries->at(i_ent)->ref_pos;
			int var_length = strlen(this->entries->at(i_ent)->ref_allele->str());
		
			//if(strlen(this->entries->at(i_ent)->ref_allele) > strlen(this->entries->at(i_ent)->alt_allele))
			//{
			//	var_length = strlen(this->entries->at(i_ent)->ref_allele) - strlen(this->entries->at(i_ent)->alt_allele) + 1;
			//}
			//else
			//{
			//	var_length = strlen(this->entries->at(i_ent)->alt_allele) - strlen(this->entries->at(i_ent)->ref_allele) + 1;
			//}

			int end = start + var_length;

			//fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\t%s\t%s\n", this->entries->at(i_ent)->i_chrom,
			//	start, end, this->entries->at(i_ent)->ref_allele, this->entries->at(i_ent)->alt_allele);

			fprintf(f_bed, "chr%s\t%d\t%d\t.\t1000\t+\n", this->entries->at(i_ent)->chrom,
				start, end);
		}
	} // i_ent loop.
	fclose(f_bed);
}


char* copy_nucs(char* seq)
{
	char* copy = new char[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(char));

	for (int i = 0; i < strlen(seq); i++)
	{
		copy[i] = seq[i];
	} // i loop.

	return(copy);
}

int* copy_bps(char* seq, int* bps)
{
	int* copy = new int[strlen(seq) + 2];
	memset(copy, 0, strlen(seq) + 2 * sizeof(int));

	for (int i = 0; i < strlen(seq); i++)
	{
		copy[i] = bps[i];
	} // i loop.

	return(copy);
}

void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	if (snp_vcf_entry->alt_alleles->size() > 1)
	{
		total_allele_count = 0;
		return;
	}

	t_string* info_string = new t_string(snp_vcf_entry->info_str);
	t_string_tokens* info_str_tokens = info_string->tokenize_by_chars("=;");

	int i_tok = 0;
	while (i_tok < info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = (info_str_tokens->at(i_tok + 1)->str())[0];
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	delete info_string;
	t_string::clean_tokens(info_str_tokens);
}

void parse_variants_per_info_string(t_vcf_entry* snp_vcf_entry,
	char* ancestral_allele,
	char* derived_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	if (snp_vcf_entry->alt_alleles->size() > 1)
	{
		total_allele_count = 0;
		return;
	}

	t_string* info_string = new t_string(snp_vcf_entry->info_str);
	t_string_tokens* info_str_tokens = info_string->tokenize_by_chars("=;");

	int i_tok = 0;
	while (i_tok < info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			//ancestral_allele = (info_str_tokens->at(i_tok+1)->str())[0];
			strcpy(ancestral_allele, info_str_tokens->at(i_tok + 1)->str());
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	if (t_string::compare_strings_ci(ancestral_allele, snp_vcf_entry->ref_allele->str()))
	{
		// Reference allele is the ancestral allele.
		strcpy(derived_allele, snp_vcf_entry->alt_alleles->at(0)->str());
	}
	else
	{
		// Reference allele is the derived allele.
		strcpy(derived_allele, snp_vcf_entry->ref_allele->str());
	}
}

bool verify_snp(t_vcf_entry* snp, char* ncrna_nucs, char ncrna_strand, int snp_i_seq)
{
	// If there are multiple alternate alleles, do not process this snp.
	if (snp->alt_alleles->size() > 1)
	{
		fprintf(stderr, "Invalid SNP @ %s:%d: Multiple alternate alleles.\n", snp->chrom, snp->ref_pos);
		return(false);
	}

	// Check if the reference allele matches the ncrna's nucleotide at the position.
	char ncrna_dna_nuc = get_transcribing_dna_nuc_per_rna_nuc(ncrna_nucs[snp_i_seq]);
	char pos_strand_nuc = toupper(ncrna_dna_nuc);
	if (ncrna_strand == '-')
	{
		pos_strand_nuc = get_complementary_dna_nuc_per_dna_nuc(pos_strand_nuc);
	}

	// Check whether the reference allele matches the sequence data. Note that the 
	char ref_allele_dna_nuc = (snp->ref_allele->str())[0];
	char alt_allele_dna_nuc = (snp->alt_alleles->at(0)->str())[0];

	if (!is_valid_nuc_char(ref_allele_dna_nuc) ||
		!is_valid_nuc_char(alt_allele_dna_nuc))
	{
		return(false);
	}

	if (ref_allele_dna_nuc == pos_strand_nuc)
	{
		// Check whether the ancestral allele matches one of variants or the reference allele.
		char reference_allele = (snp->ref_allele->str())[0];
		char alternate_allele = (snp->alt_alleles->at(0)->str())[0];

		char ancestral_allele = 0;
		int alt_all_cnt = 0;
		int total_all_cnt = 0;
		parse_variants_per_info_string(snp, ancestral_allele, alt_all_cnt, total_all_cnt);

		if (toupper(ancestral_allele) != toupper(reference_allele) &&
			toupper(ancestral_allele) != toupper(alternate_allele))
		{
			fprintf(stderr, "Invalid SNP @ %s:%d(%d): Ancestral allele does not match neither ref nor alternate allele. %c: %c, %c (%s).\n", snp->chrom, snp->ref_pos, snp_i_seq,
				ancestral_allele, reference_allele, alternate_allele, snp->info_str);
			return(false);
		}

		// The snp seems to be ok to process.
		return(true);
	}
	else
	{
		// The SNP's reference allele does not match the reference nucleotide in the ncrna. Something is messed up in snp entry.
		fprintf(stderr, "Invalid SNP @ %s:%d: SNP's reference allele does not match (%c, %c) the hg19 reference sequence:\n%s\n%d\n", snp->chrom, snp->ref_pos, ref_allele_dna_nuc, pos_strand_nuc, ncrna_nucs, snp_i_seq);
		return(false);
	}
}

void parse_variants_per_info_string(char* info_str,
	char& ancestral_allele,
	int& alternate_allele_count,
	int& total_allele_count)
{
	t_string_tokens* info_str_tokens = t_string::tokenize_by_chars(info_str, "=;");

	int i_tok = 0;
	while (i_tok < info_str_tokens->size())
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = (info_str_tokens->at(i_tok + 1)->str())[0];
		} // read the ancestral allele.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		} // ancestral count check.
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_allele_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}

		i_tok++;
	} // info string tokens loop.

	t_string::clean_tokens(info_str_tokens);
}


void get_mutation_matrices_per_DAF(char* snp_vcf_fp, double min_daf, double max_daf)
{
	int** mutation_cnts = new int*[4];
	for (int i = 0; i < 4; i++)
	{
		mutation_cnts[i] = new int[4];
		for (int j = 0; j < 4; j++)
		{
			mutation_cnts[i][j] = 0;
		} // j loop.
	} // i loop.

	  //t_VCF* vcf = new t_VCF(snp_vcf_fp);

	FILE* f_vcf = open_f(snp_vcf_fp, "r");

	//for(int i_snp = 0; i_snp < vcf->entries->size(); i_snp++)
	int i_line = 0;
	while (1)
	{
		if (i_line % 100000 == 0)
		{
			fprintf(stderr, "Processing %d. variant.       \r", i_line);
		}

		char* cur_line = getline(f_vcf);

		if (cur_line == NULL)
		{
			break;
		}

		//t_vcf_entry* cur_snp = vcf->entries->at(i_snp);
		// Get the reference and alternate alleles.
		// 10      68533   .       T       G       .       PASS    AA=.;AC=19;AN=120;DP=80 GT:DP:CB
		char chrom[100];
		char posn[100];
		char dot[100];
		char ref_allele[100];
		char alt_allele[100];
		char dot2[100];
		char pass_str[100];
		char info_str[100];
		if (sscanf(cur_line, "%s %s %s %s %s %s %s %s", chrom, posn, dot, ref_allele, alt_allele, dot2, pass_str, info_str) != 8)
		{
			fprintf(stderr, "Could not parse %s\n", cur_line);
			exit(0);
		}

		char anc_all = 0;
		char der_all = 0;
		int alt_all_cnt = 0;
		int anc_all_cnt = 0;
		int total_all_cnt = 0;

		parse_variants_per_info_string(info_str,
			anc_all,
			alt_all_cnt,
			total_all_cnt);

		// Set the derived allele.
		if (anc_all == ref_allele[0])
		{
			der_all = alt_allele[0];
		}
		else
		{
			der_all = ref_allele[0];
		}

		// Valid nucleotide check on the ancestral allele.
		if (is_valid_nuc_char(anc_all))
		{
			// Make sure that the ancestral allele is either the 
			if (toupper(anc_all) == toupper(ref_allele[0]) ||
				toupper(anc_all) == toupper(alt_allele[0]))
			{
				if (toupper(anc_all) == toupper(alt_allele[0]))
				{
					der_all = toupper(ref_allele[0]);

					// Ancestral allele is the alternate allele. The derived allele is the reference allele.
					anc_all_cnt = alt_all_cnt;
				}
				else
				{
					der_all = toupper(alt_allele[0]);

					anc_all_cnt = total_all_cnt - alt_all_cnt;
				}

				if (is_valid_nuc_char(der_all))
				{
					double der_all_freq = (double)(total_all_cnt - anc_all_cnt) / total_all_cnt;

					if (max_daf > der_all_freq &&
						min_daf < der_all_freq)
					{
						//fprintf(stderr, "%s:%d, %s: %c -> %c (%lf)\n", cur_snp->chrom, cur_snp->ref_pos, cur_snp->info_str, anc_all, der_all, der_all_freq);
						mutation_cnts[nuc_2_num(anc_all)][nuc_2_num(der_all)]++;
					}
				}
			}
		} // validity check on the ancestral allele.

		delete[] cur_line;

		i_line++;
	} // i_snp loop.

	fclose(f_vcf);

	int n_total_trans = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			n_total_trans += mutation_cnts[i][j];
		} // j loop.
	} // i loop.	

	fprintf(stderr, "\n%d nucleotide transitions:\n", n_total_trans);

	fprintf(stderr, "  ");
	for (int j = 0; j < 4; j++)
	{
		fprintf(stderr, "%4c ", num_2_nuc(j));
	}
	fprintf(stderr, "\n");
	for (int i = 0; i < 4; i++)
	{
		fprintf(stderr, "%c ", num_2_nuc(i));
		for (int j = 0; j < 4; j++)
		{
			//fprintf(stderr, "%4d", mutation_cnts[i][j]);
			fprintf(stderr, "%.2lf ", (double)mutation_cnts[i][j] / n_total_trans);
		} // j loop.

		fprintf(stderr, "\n");
	} // i loop.	
}


char mutate_nuc_per_transition_transversion(char nuc)
{
	if (toupper(nuc) == 'A')
	{
		return('G');
	}
	else if (toupper(nuc) == 'C')
	{
		return('T');
	}
	else if (toupper(nuc) == 'G')
	{
		return('A');
	}
	else if (toupper(nuc) == 'T')
	{
		return('C');
	}
	else
	{
		fprintf(stderr, "WTF??\n");
		exit(0);
	}
}

char mutate_nuc(char nuc, t_rng* rng)
{
	char cur_rand_nuc = toupper(nuc);
	while (cur_rand_nuc == toupper(nuc))
	{
		int rand_i = (int)(floor(4 * rng->random_double_ran3()));
		char rand_nucs[] = "ACGT";
		cur_rand_nuc = rand_nucs[rand_i];
	}

	return(cur_rand_nuc);
}

bool mutate_GC_2_AU(char* seq)
{
	// Count GC content.
	int n_gc = 0;
	for (int i = 0; i < strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			n_gc++;
		}
	} // i loop.

	if (n_gc == 0)
	{
		return(false);
	}

	int i_rand = rand() % n_gc;

	int i_gc = 0;
	for (int i = 0; i < strlen(seq); i++)
	{
		if (toupper(seq[i]) == 'G' ||
			toupper(seq[i]) == 'C')
		{
			if (i_gc == i_rand)
			{
				// Choose a nuc.
				int n_rand = rand() % 2;
				char rand_nuc = (n_rand == 0) ? ('A') : ('U');
				fprintf(stderr, "%c->%c @ %d\n", seq[i], rand_nuc, i);
				seq[i] = rand_nuc;
				break;
			}

			i_gc++;
		}
	} // i loop.

	return(true);
}

char mutate_nuc_2_GC(char nuc, t_rng* rng)
{
	if (toupper(nuc) != 'G' &&
		toupper(nuc) != 'C')
	{
		int val = (int)(floor(2 * rng->random_double_ran3()));
		if (val == 0)
		{
			return('G');
		}
		else
		{
			return('C');
		}
	}
	else
	{
		return(toupper(nuc));
	}
}

char* get_random_seq(int l)
{
	char* seq = new char[l + 2];
	memset(seq, 0, (l + 1) * sizeof(char));
	for (int i = 0; i < l; i++)
	{
		int random = rand() % 4;
		seq[i] = num_2_nuc(random);
	} // i loop

	return(seq);
}

///*
//Do a mutation of {G,C}->{A,U} in one sequence.
//*/
//void length_based_Duplex_deltaG_simulation()
//{
//	int l_match = 7;
//
//	// Do 1000 simulations.
//	FILE* f_diff_dup_energy_per_l_dup = open_f("diff_dup_energy_per_l_dup.txt", "w");
//	fclose(f_diff_dup_energy_per_l_dup);
//
//	f_diff_dup_energy_per_l_dup = open_f("diff_dup_energy_per_l_dup.txt", "a");
//	for(int l_duplex = 20; l_duplex < 70; l_duplex+=5)
//	{
//		vector<double>* current_energy_diffs = new vector<double>();
//		for(int i = 0; i < 100; i++)
//		{
//			fprintf(stderr, "%d. seq. pair (l_dup=%d).\n", i, l_duplex);
//
//			// Generate the random sequence.
//			char* cur_random_seq1 = get_random_seq(l_duplex); 
//
//			// Generate the second random sequence.
//			char* matching_seq = new char[l_match + 2];
//
//			for(int i = 0; i < l_match; i++)
//			{
//				matching_seq[i] = get_complementary_rna_nuc_per_rna_nuc(cur_random_seq1[i]);
//			} // i loop.
//
//			t_string* matching_seq_str = new t_string(matching_seq);
//			matching_seq_str->revert();
//
//			char* prefix = get_random_seq(l_duplex - l_match);
//
//			t_string* cur_random_seq2_str = new t_string(prefix);
//			cur_random_seq2_str->concat_string(matching_seq_str);
//
//			char* cur_random_seq2 = cur_random_seq2_str->str();
//
//			fprintf(stderr, "DuplexFolding:\n%s\n%s\n", cur_random_seq1, cur_random_seq2);
//			double cur_orig_dup_energy;
//			get_DuplexFold_energy(cur_random_seq1, cur_random_seq2, cur_orig_dup_energy);
//
//			// Mutate a G or C into A or U.
//			char* variant_seq = new char[strlen(cur_random_seq2) + 1];
//			strcpy(variant_seq, cur_random_seq2);
//
//			// Make sure that sequence is mutatable.
//			if(mutate_GC_2_AU(variant_seq))
//			{
//				fprintf(stderr, "DuplexFolding:\n%s\n%s\n", cur_random_seq1, variant_seq);
//				double cur_var_dup_energy;
//				get_DuplexFold_energy(cur_random_seq1, variant_seq, cur_var_dup_energy);
//
//				getc(stdin);
//
//				current_energy_diffs->push_back(cur_var_dup_energy - cur_orig_dup_energy);
//			}
//			// Free memory.
//			delete matching_seq_str;
//			delete [] prefix;
//			delete cur_random_seq2_str;
//			delete [] cur_random_seq1;
//			delete [] variant_seq;
//		} // i loop.
//
//		double total_diff_en = 0.0;
//		for(int i = 0; i < current_energy_diffs->size(); i++)
//		{
//			total_diff_en += current_energy_diffs->at(i);
//		} // i loop.
//
//		// Compute the average of differences.
//		fprintf(f_diff_dup_energy_per_l_dup, "%d\t%lf\t%d\n", l_duplex, total_diff_en, current_energy_diffs->size());
//
//		delete current_energy_diffs;
//	} // l_duplex loop.
//
//	fclose(f_diff_dup_energy_per_l_dup);
//}

bool get_DAF_per_snp_region(t_annot_region* snp_region, double& daf)
{
	// 1       2337697 rs55744642;-/CGACA,venter,watson        GG      GCGACAG,GGGACAG 100     PASS    CA=0;DB;NF=16;NR=14;HP=3;DP=179 GT:GQ   2|1:100
	t_vcf_info* snp_vcf_info = (t_vcf_info*)(snp_region->data);
	char* snp_info_str = snp_vcf_info->info_str;

	char ref_allele_snp = snp_vcf_info->ref_allele_str[0];
	t_string_tokens* info_str_tokens = t_string::tokenize_by_chars(snp_info_str, ";=");

	// Get the AA (ancestral allele) and AC (alternate count) entries.
	char ancestral_allele = 0;
	int alternate_count = 0;
	int total_count = 0;
	for (int i_tok = 0; i_tok < info_str_tokens->size(); i_tok++)
	{
		if (strcmp(info_str_tokens->at(i_tok)->str(), "AA") == 0)
		{
			ancestral_allele = ((char*)(info_str_tokens->at(i_tok + 1)->str()))[0];
		}
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AC") == 0)
		{
			alternate_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}
		else if (strcmp(info_str_tokens->at(i_tok)->str(), "AN") == 0)
		{
			total_count = atoi(info_str_tokens->at(i_tok + 1)->str());
		}
	} // i_tok loop.

	t_string::clean_tokens(info_str_tokens);

	// Is the ancestral allele resolved?
	if (ancestral_allele == '.' || ancestral_allele == 0)
	{
		return(false);
	}

	// Resolve the ancestral and derived counts.
	int reference_count = total_count - alternate_count;
	if (ref_allele_snp == ancestral_allele)
	{
		daf = (double)reference_count / total_count;
	}
	else
	{
		daf = (double)(total_count - reference_count) / total_count;
	}

	return(true);
}


