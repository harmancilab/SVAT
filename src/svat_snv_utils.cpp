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

bool __DUMP_CRYPTANNOT_SNV_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

// Write the ternary signal on the EOI=gene_exons.
void signalize_VCF_SNV_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
	char* per_chrom_VCF_dir,
	char* op_dir)
{
	// Load the EOIs.
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	double total_EOI_covg = coverage(EOI_regs);
	fprintf(stderr, "Loaded %d EOI regions covering %d nucleotides.\n", (int)(EOI_regs->size()), (int)total_EOI_covg);
	vector<char*>* all_EOI_ids = new vector<char*>();
	for (int i_e = 0; i_e < EOI_regs->size(); i_e++)
	{
		all_EOI_ids->push_back(t_string::copy_me_str(EOI_regs->at(i_e)->name));
	} // i_e loop.

	  // These are the gene id's that will be used in impact assignments.
	vector<char*>* factorized_EOI_ids = t_string::get_unique_entries(all_EOI_ids);

	fprintf(stderr, "%d factorized gene ids.\n", factorized_EOI_ids->size());
	int n_bits_per_element_id = (int)(ceil(log(factorized_EOI_ids->size()) / log(2)));
	fprintf(stderr, "%d bits per element.\n", n_bits_per_element_id);
	if (n_bits_per_element_id > 16)
	{
		fprintf(stderr, "ERROR: Cannot use %d bits per element; too many elements?\n", n_bits_per_element_id);
		exit(0);
	}

	// The EOI's determine the coordinate system that we will use. The impact value should include this.
	t_restr_annot_region_list* restr_EOI_regs = restructure_annot_regions(EOI_regs);

	// Set the chromosome id's.
	vector<char*>* chr_ids = restr_EOI_regs->chr_ids;

	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_EOI_regs->regions_per_chrom[i_chr];
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		fprintf(stderr, "Allocating per allele impact signal.\n");
		unsigned int** per_allele_VCF_signal = new unsigned int*[5];
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			per_allele_VCF_signal[i_allele] = new unsigned int[cur_chr_EOI_covg + 2];
			memset(per_allele_VCF_signal[i_allele], 0, sizeof(unsigned int) * cur_chr_EOI_covg);
		} // i_chr loop.

		  // Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
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
			fprintf(stderr, "Reading SNVs in %s\n", cur_chrom_vcf_fp);
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

			t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");

			// Parse the location.
			char* var_chr = t_string::copy_me_str(cur_line_toks->at(CHROM_col_i)->str());
			int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str());

			// Following takes only bi-allelic SNVs.
			int eoi_reg_chr_i = t_string::get_i_str(restr_EOI_regs->chr_ids, var_chr);
			if (eoi_reg_chr_i == restr_EOI_regs->chr_ids->size() ||
				t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str()) != 1 ||
				t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str()) != 1)

			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// We take the first allele in every SNP.
			char ref_allele = cur_line_toks->at(REF_ALLELE_col_i)->str()[0];
			char alt_allele = cur_line_toks->at(ALT_ALLELE_col_i)->str()[0];

			char ref_alt_allele[10];
			ref_alt_allele[0] = ref_allele;
			ref_alt_allele[1] = alt_allele;

			// Generate the per allele counts from genotype value.
			char* geno_str = cur_line_toks->at(GENOTYPE_col_i)->str();
			int per_allele_cnt[5];
			memset(per_allele_cnt, 0, sizeof(int) * 5);
			per_allele_cnt[nuc_2_num(ref_alt_allele[geno_str[0] - '0'])] = 1;
			per_allele_cnt[nuc_2_num(ref_alt_allele[geno_str[2] - '0'])] = 1;

			if (geno_str[0] - '0' > 1 ||
				geno_str[2] - '0' > 1)
			{
				fprintf(stderr, "Invalid allele code: %s: %c, %c; %d, %d",
						cur_line,
						geno_str[0],
						geno_str[2],
						geno_str[0] - '0',
						geno_str[2] - '0');
				exit(0);
			}

			// Search the current variant.
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			// For every overlap, update the allelic signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				bool overlap = (cur_chr_EOI_regs->at(cur_EOI_reg_i)->start <= var_posn &&
					cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn);
				if (overlap)
				{
					// For each overlap, update ref and alternate alleles.
					int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

					for (int allele_i = 0; allele_i < 4; allele_i++)
					{
						per_allele_VCF_signal[nuc_2_num(alt_allele)][track_posn] = 1;
					} // allele_i loop.

					found_overlap = true;

					if (__DUMP_CRYPTANNOT_SNV_MSGS__)
					{
						fprintf(stderr, "Overlap %s: %s:%d-%d:%s\n", cur_line,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);
					}
				} // overlap check.

				cur_EOI_reg_i++;
			} // region loop.

			if (!found_overlap && __DUMP_CRYPTANNOT_SNV_MSGS__)
			{
				fprintf(stderr, "Could not find an overlap for %s\n", cur_line);
			}

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			delete[] cur_line;
		} // file reading loop.
		close_f(f_vcf, cur_chrom_vcf_fp);

		fprintf(stderr, "\nFinished reading VCF file.\n");

		// Save the files.
		char cur_chr_VCF_signal_op_fp[1000];
		sprintf(cur_chr_VCF_signal_op_fp, "%s/variant_signal_%s.bin.gz", op_dir, chr_ids->at(i_chr));
		FILE* f_cur_chr_VCF_signal = open_f(cur_chr_VCF_signal_op_fp, "wb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fwrite(per_allele_VCF_signal[i_allele], sizeof(unsigned int), cur_chr_EOI_covg, f_cur_chr_VCF_signal);
			delete[] per_allele_VCF_signal[i_allele];
		} // i_chr loop.
		delete[] per_allele_VCF_signal;
		close_f(f_cur_chr_VCF_signal, cur_chr_VCF_signal_op_fp);

		// Write the BED file with sorted regions that match the signal coordinates in signal file.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, cur_chr_EOI_regs);
	} // i_chr loop.

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", op_dir);
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.
	fclose(f_chrs);
} // signalize_VCF_SNV_genotypes_per_EOI_regs


  // https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
  // EOI setting from the VEP output.
  // We have the EOIs; each EOI has an id.
  // EOI id's must match with VEP element id's.
void signalize_VEP_annotated_SNVs_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
	char* genome_dir,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir)
{
	// Load the EOIs.
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	double total_EOI_covg = coverage(EOI_regs);

	// Load the impact values from VEP's output.
	//vector<char*>* sorted_impact_vals = buffer_file(sorted_impact_value_strings_list_fp);
	vector<char*>* sorted_impact_vals = new vector<char*>();
	vector<t_VEP_term_ctx*>* vep_term_ctx = load_VEP_annotation_term_context(sorted_impact_value_strings_list_fp, sorted_impact_vals);
	fprintf(stderr, "Loaded %d impact strings.\n", sorted_impact_vals->size());

	// The EOI's determine the coordinate system that we will use. The impact value should include this.
	t_restr_annot_region_list* restr_EOI_regs = restructure_annot_regions(EOI_regs);

	// Set the chromosome id's.
	vector<char*>* chr_ids = restr_EOI_regs->chr_ids;

	int n_bits_per_nuc = 3;
	int n_neigh_nucs = 6;

	fprintf(stderr, "Allocating impact signals per chromosome.\n");
	char unmatched_annotations_fp[] = "unmatched_vep_impacts.op.gz";
	FILE* f_unmatched = open_f(unmatched_annotations_fp, "w");
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", chr_ids->at(i_chr));

		char chr_chrom_seq_fp[1000];
		sprintf(chr_chrom_seq_fp, "%s/%s.bin.gz", genome_dir, chr_ids->at(i_chr));
		fprintf(stderr, "Loading the chromosome sequence from %s.\n", chr_chrom_seq_fp);
		int l_seq = 0;
		char* cur_chrom_seq = load_binary_sequence_file(chr_chrom_seq_fp, l_seq);
		fprintf(stderr, "Loaded %d long chromosome sequence for %s\n", l_seq, chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_EOI_regs->regions_per_chrom[i_chr];
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		fprintf(stderr, "Allocating per allele impact signal.\n");
		unsigned long long** per_allele_VEP_signal = new unsigned long long*[5];
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			per_allele_VEP_signal[i_allele] = new unsigned long long[cur_chr_EOI_covg + 2];
			memset(per_allele_VEP_signal[i_allele], 0, sizeof(unsigned long long) * cur_chr_EOI_covg);
		} // i_chr loop.

		  // Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		  // Sanity check: Last should be last.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// Initialize the SNV annotation vector with the neighboring nucleotide signal.
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			for (int var_posn = cur_chr_EOI_regs->at(i_reg)->start;
				var_posn <= cur_chr_EOI_regs->at(i_reg)->end;
				var_posn++)
			{
				int track_posn = (cur_chr_EOI_regs->at(i_reg)->score + var_posn - cur_chr_EOI_regs->at(i_reg)->start);

				// Copy the surrounding tri-nucleotides: We store the 6 nucleotides in the neighborhood, 32 bit value is more than enough.
				unsigned int neigh_seq_signal = 0; // 32-bit neighbor value.
				int cur_base = 0;
				for (int i = var_posn - 3; i < var_posn; i++)
				{
					unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
					neigh_seq_signal = (neigh_seq_signal | (cur_nuc_val << cur_base));
					cur_base += n_bits_per_nuc;
				} // i_l loop.

				for (int i = var_posn + 1; i <= var_posn + 3; i++)
				{
					unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
					neigh_seq_signal = (neigh_seq_signal | (cur_nuc_val << cur_base));
					cur_base += n_bits_per_nuc;
				} // i_l loop.

				// Assign score to the EOI impact score track; set the coding frame and the impact values to be none as they do not mean anything for this position.
				// This entry contains only the sequence information.
				unsigned long long cur_impact_val = 0;
				int coding_frame = 0;
				unsigned long long cur_signal = pack_set_impact_signal_values(coding_frame,
					cur_impact_val,
					neigh_seq_signal, n_bits_per_nuc, n_neigh_nucs);

				per_allele_VEP_signal[0][track_posn] = cur_signal;
				per_allele_VEP_signal[1][track_posn] = cur_signal;
				per_allele_VEP_signal[2][track_posn] = cur_signal;
				per_allele_VEP_signal[3][track_posn] = cur_signal;
			} // var_posn loop.
		} // i_reg loop.

		char vep_op_fp[1000];
		sprintf(vep_op_fp, "%s/%s.vep.gz", vep_op_dir, chr_ids->at(i_chr));
		FILE* f_vep_op = open_f(vep_op_fp, "r");
		if (!vep_op_fp)
		{
			fprintf(stderr, "Could not find vep output file @ %s\n", vep_op_fp);
			exit(0);
		}

		unsigned long long LLONE = 1;

		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		int EOI_element_VEP_col_i = -1;
		if (EOI_element_type == VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION)
		{
			EOI_element_VEP_col_i = TRANSCRIPT_ID_col_i;
		}
		else if (EOI_element_type == VEP_SIGNALIZE_GENE_SUMMARIZATION)
		{
			EOI_element_VEP_col_i = GENE_ID_col_i;
		}

		fprintf(stderr, "Reading VEP output file from %s.\n", vep_op_fp);
		int n_lines = 0;
		while (1)
		{
			char* cur_line = getline(f_vep_op);
			if (cur_line == NULL)
			{
				break;
			}

			n_lines++;
			if (n_lines % 10000 == 0)
			{
				fprintf(stderr, "Processed %d lines.            \r", n_lines);
			}

			// Skip comment lines.
			if (cur_line[0] == '#')
			{
				delete[] cur_line;
				continue;
			}

			t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");

			char* vep_element_id = cur_line_toks->at(EOI_element_VEP_col_i)->str();

			// Parse the location.
			t_string_tokens* location_str_toks = cur_line_toks->at(LOCATION_col_i)->tokenize_by_chars(":-");
			if (location_str_toks->size() != 2)
			{
				fprintf(stderr, "Sanity check failed: Expecting 2 strings in location column: %s\n", cur_line);
				exit(0);
			}

			char* var_name = cur_line_toks->at(VAR_ID_col_i)->str();
			char* var_chr = t_string::copy_me_str(location_str_toks->at(0)->str());
			int var_posn = atoi(location_str_toks->at(1)->str());
			t_string::clean_tokens(location_str_toks);

			// Set the coding frame.
			int coding_frame = 0;
			char* cds_posn_str = cur_line_toks->at(CDS_POSN_col_i)->str();
			if (t_string::is_number(cds_posn_str) &&
				!t_string::compare_strings(cds_posn_str, "-"))
			{
				int cds_posn = atoi(cds_posn_str);
				coding_frame = ((cds_posn - 1) % 3);
			}

			int eoi_reg_chr_i = t_string::get_i_str(restr_EOI_regs->chr_ids, var_chr);
			if (eoi_reg_chr_i == restr_EOI_regs->chr_ids->size())
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			if (!t_string::compare_strings(var_chr, chr_ids->at(i_chr)))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			vector<t_annot_region*>* cur_chr_EOI_regs = restr_EOI_regs->regions_per_chrom[eoi_reg_chr_i];

			int allele_i = nuc_2_num(cur_line_toks->at(ALLELE_col_i)->str()[0]);

			if (__DUMP_CRYPTANNOT_SNV_MSGS__)
			{
				fprintf(stderr, "allele_i: %d for %s\n", allele_i, cur_line);
			}

			// Go over the consequence tokens and choose the highest impact.
			t_string_tokens* conseq_toks = cur_line_toks->at(CONSEQUENCE_col_i)->tokenize_by_chars(", ");

			unsigned long long cur_impact_val = 0;
			for (int i_tok = 0; i_tok < conseq_toks->size(); i_tok++)
			{
				int cur_tok_impact_val = t_string::get_i_str(sorted_impact_vals, conseq_toks->at(i_tok)->str());

				if (cur_tok_impact_val == sorted_impact_vals->size())
				{
					fprintf(stderr, "Could not parse impact string: %s\n", cur_line_toks->at(CONSEQUENCE_col_i)->str());
					exit(0);
				}

				cur_impact_val = cur_impact_val | (LLONE << cur_tok_impact_val);
			} // i_tok loop.

			// Clean tokens.
			t_string::clean_tokens(conseq_toks);

			// Search the current variant.
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			// Copy the surrounding tri-nucleotides: We store the 6 nucleotides in the neighborhood, 32 bit value is more than enough.
			unsigned int neigh_seq_signal = 0; // 32-bit neighbor value.
			int cur_base = 0;
			for (int i = var_posn - 3; i < var_posn; i++)
			{
				unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
				neigh_seq_signal = (neigh_seq_signal | (cur_nuc_val << cur_base));
				cur_base += n_bits_per_nuc;
			} // i_l loop.

			for (int i = var_posn + 1; i <= var_posn + 3; i++)
			{
				unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
				neigh_seq_signal = (neigh_seq_signal | (cur_nuc_val << cur_base));
				cur_base += n_bits_per_nuc;
			} // i_l loop.

			// For every overlap, update the allelic signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				bool overlap = (cur_chr_EOI_regs->at(cur_EOI_reg_i)->start <= var_posn &&
					cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn);

				if (overlap)
				{
					int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);
					if (track_posn > cur_chr_EOI_covg)
					{
						fprintf(stderr, "The position is further from the end of track: %s:%d, Current element start: %d, track length: %d\n",
							chr_ids->at(i_chr), var_posn,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							cur_chr_EOI_covg);

						exit(0);
					}

					// Make sure the impact assigned element is the same element that we are setting the signal at.
					if (t_string::compare_strings(vep_element_id, cur_chr_EOI_regs->at(cur_EOI_reg_i)->name))
					{
						//unsigned int prev_impact_val = (unsigned int)((deletion_VEP_signal[track_posn] >> 18));
						char prev_neigh_seq_buff[100];
						unsigned int prev_coding_frame = 0;
						unsigned long long prev_impact_val = 0;
						unpack_impact_signal_values(per_allele_VEP_signal[allele_i][track_posn], prev_coding_frame, prev_impact_val, prev_neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						// If the new impact is higher, assign it.
						cur_impact_val = cur_impact_val | prev_impact_val;

						if (__DUMP_CRYPTANNOT_SNV_MSGS__)
						{
							fprintf(stderr, "%s/%d::%s:%d-%d:: Coding Frame @ SNV: %d (track posn: %d;; allele impact: %llX)\n",
								var_name, allele_i,
								var_chr,
								var_posn,
								var_posn,
								coding_frame,
								track_posn,
								cur_impact_val);
						}

						// Assign score to the EOI impact score track.	
						unsigned long long cur_signal = pack_set_impact_signal_values(coding_frame,
							cur_impact_val,
							neigh_seq_signal, n_bits_per_nuc, n_neigh_nucs);

						per_allele_VEP_signal[allele_i][track_posn] = cur_signal;

						// Sanity check to make sure packing worked.
						char neigh_seq_buff[100];
						unsigned int s_coding_frame = 0;
						unsigned long long s_impact_val = 0;
						unpack_impact_signal_values(cur_signal, s_coding_frame, s_impact_val, neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						if (s_impact_val != cur_impact_val ||
							s_coding_frame != coding_frame)
						{
							fprintf(stderr, "Sanity check failed: Coding frame or the impact value is not coded correctly:\n\
%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s (%u)\n",
cur_line,
s_impact_val,
s_coding_frame,
neigh_seq_buff, neigh_seq_signal);
							exit(0);
						}

						if (__DUMP_CRYPTANNOT_SNV_MSGS__)
						{
							fprintf(stderr, "%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s (%u)\n",
								cur_line,
								s_impact_val,
								s_coding_frame,
								neigh_seq_buff, neigh_seq_signal);
						}

						found_overlap = true;

						if (__DUMP_CRYPTANNOT_SNV_MSGS__)
						{
							fprintf(stderr, "Overlap %s: %s:%d-%d:%s\n", cur_line,
								cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
								cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
								cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
								cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);
						}
					} // gene id comparison.
				} // overlap check.

				cur_EOI_reg_i++;
			} // region loop.

			// If overlap not found, write to the unmatched impact entries.
			if (!found_overlap)
			{
				fprintf(stderr, "Could not find an overlap for %s; %s:%d\n", cur_line, vep_element_id, var_posn);
				fprintf(f_unmatched, "%s\n", cur_line);
			}

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			delete[] cur_line;
		} // file reading loop.
		close_f(f_vep_op, vep_op_fp);

		fprintf(stderr, "Finished reading, saving the VEP impact signals on %s.\n", chr_ids->at(i_chr));

		// Save the files.
		char cur_chr_impacts_op_fp[1000];
		sprintf(cur_chr_impacts_op_fp, "%s/impact_signal_%s.bin.gz", op_dir, chr_ids->at(i_chr));
		FILE* f_vep_impacts = open_f(cur_chr_impacts_op_fp, "wb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Writing allele %d\n", i_allele);
			fwrite(per_allele_VEP_signal[i_allele], sizeof(unsigned long long), cur_chr_EOI_covg, f_vep_impacts);
			delete[] per_allele_VEP_signal[i_allele];
		} // i_allele loop.
		close_f(f_vep_impacts, cur_chr_impacts_op_fp);
		delete[] per_allele_VEP_signal;

		// Write the BED file with sorted regions that match the signal coordinates in signal file.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, cur_chr_EOI_regs);

		delete[] cur_chrom_seq;
	} // i_chr loop.
	close_f(f_unmatched, unmatched_annotations_fp);

	// Write the impact values.
	char op_sorted_impact_value_strings_list_fp[1000];
	sprintf(op_sorted_impact_value_strings_list_fp, "%s/impact_value_strings.list", op_dir);
	copy_file(sorted_impact_value_strings_list_fp, op_sorted_impact_value_strings_list_fp);

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", op_dir);
	FILE* f_chrs = open_f(chr_ids_fp, "w");
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		fprintf(f_chrs, "%s\n", chr_ids->at(i_chr));
	} // i_chr loop.
	fclose(f_chrs);
}

void extract_annotated_SNVs_from_signal(char* annotated_variant_signal_dir, char* op_fp)
{
	fprintf(stderr, "Extracting the annotated SNV signals from %s to %s.\n", annotated_variant_signal_dir, op_fp);

	// Load the annotation signal.
	char anno_var_sig_chr_ids_fp[1000];
	sprintf(anno_var_sig_chr_ids_fp, "%s/chr_ids.txt", annotated_variant_signal_dir);
	vector<char*>* anno_var_chr_ids = buffer_file(anno_var_sig_chr_ids_fp);
	if (anno_var_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotated variants chromosome id's from %s\n", anno_var_sig_chr_ids_fp);
		exit(0);
	}

	FILE* f_op = open_f(op_fp, "w");
	for (int i_chr = 0; i_chr < anno_var_chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", anno_var_chr_ids->at(i_chr));

		// First load the EOI regions.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));

		if (!check_file(cur_sorted_BED_fp))
		{
			fprintf(stderr, "Skipping %s, no impact signals for this chromosome.\n", anno_var_chr_ids->at(i_chr));
		}

		vector<t_annot_region*>* sorted_cur_chr_EOI_regs = load_BED(cur_sorted_BED_fp);
		fprintf(stderr, "Loaded %d EOI regions.\n", sorted_cur_chr_EOI_regs->size());

		// Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < sorted_cur_chr_EOI_regs->size(); i_reg++)
		{
			sorted_cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (sorted_cur_chr_EOI_regs->at(i_reg)->end - sorted_cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.
		fprintf(stderr, "Loaded %d EOI regions on %s\n", sorted_cur_chr_EOI_regs->size(), anno_var_chr_ids->at(i_chr));

		// Load the impact signals.
		fprintf(stderr, "Loading per allele impact signals.\n");
		int** per_allele_annotated_variant_signal = new int*[4];

		char cur_chr_impacts_op_fp[1000];
		sprintf(cur_chr_impacts_op_fp, "%s/annotated_variant_signal_%s.bin.gz", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));
		FILE* f_vep_impacts = open_f(cur_chr_impacts_op_fp, "rb");
		int cur_chr_EOI_covg = coverage(sorted_cur_chr_EOI_regs);
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Reading impacts per allele %d\n", i_allele);
			per_allele_annotated_variant_signal[i_allele] = new int[cur_chr_EOI_covg + 2];
			if(fread(per_allele_annotated_variant_signal[i_allele], sizeof(int), cur_chr_EOI_covg, f_vep_impacts) != cur_chr_EOI_covg)
			{
				fprintf(stderr, "Could not read allele_i=%d from %s\n", i_allele, cur_chr_impacts_op_fp);
				exit(0);
			}
		} // i_allele loop.

		  // Load the impact values from VEP's output.
		char sorted_impact_value_strings_list_fp[1000];
		sprintf(sorted_impact_value_strings_list_fp, "%s/impact_value_strings.list", annotated_variant_signal_dir);
		vector<char*>* sorted_impact_vals = buffer_file(sorted_impact_value_strings_list_fp);
		fprintf(stderr, "Loaded %d impact strings.\n", sorted_impact_vals->size());

		// Go over all the regions.
		for (int cur_sig_i = 1; cur_sig_i < cur_chr_EOI_covg; cur_sig_i++)
		{
			// Analyze all the alleles for this position.
			for (int i_allele = 0; i_allele < 4; i_allele++)
			{
				if (per_allele_annotated_variant_signal[i_allele][cur_sig_i] != 0)
				{
					// Find the region that corresponds to this position.
					for (int i_reg = 0; i_reg < sorted_cur_chr_EOI_regs->size(); i_reg++)
					{
						t_annot_region* cur_EOI_region = sorted_cur_chr_EOI_regs->at(i_reg);

						// Check if the signal index is within this EOI region's mapped coordinates.
						if (cur_sig_i >= cur_EOI_region->score &&
							cur_sig_i <= (cur_EOI_region->score + cur_EOI_region->end - cur_EOI_region->start))
						{
							// This region is what we overlap with.
							int genome_posn = sorted_cur_chr_EOI_regs->at(i_reg)->start + cur_sig_i - sorted_cur_chr_EOI_regs->at(i_reg)->score;

							fprintf(f_op, "%s\t%d\t%c\t%s\t%s\n",
								sorted_cur_chr_EOI_regs->at(i_reg)->chrom, genome_posn,
								num_2_nuc(i_allele),
								sorted_cur_chr_EOI_regs->at(i_reg)->name,
								sorted_impact_vals->at(per_allele_annotated_variant_signal[i_allele][cur_sig_i]));

							break;
						}
					} // i_reg loop.
				} // annotated variant signal check.
			} // i_allele loop.
		} // i loop.
	} // i_chr loop.

	  // Close output file.
	close_f(f_op, op_fp);
}

void multiply_SNV_variant_and_annotation_signals(char* annotation_signal_dir, char* variant_signal_dir, char* op_dir)
{
	// Load the annotation signal.
	char anno_sig_chr_ids_fp[1000];
	sprintf(anno_sig_chr_ids_fp, "%s/chr_ids.txt", annotation_signal_dir);
	vector<char*>* anno_chr_ids = buffer_file(anno_sig_chr_ids_fp);
	if (anno_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotation chromosome id's from %s\n", anno_sig_chr_ids_fp);
		exit(0);
	}

	// Load the variant signal.
	char var_sig_chr_ids_fp[1000];
	sprintf(var_sig_chr_ids_fp, "%s/chr_ids.txt", variant_signal_dir);
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
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", variant_signal_dir, var_chr_ids->at(i_chr));
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
		sprintf(cur_chr_annotation_signal_fp, "%s/impact_signal_%s.bin.gz", annotation_signal_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_annotation_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_annotation_signal_fp, var_chr_ids->at(i_chr));
			exit(0);
		}

		char cur_chr_variant_signal_fp[1000];
		sprintf(cur_chr_variant_signal_fp, "%s/variant_signal_%s.bin.gz", variant_signal_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_variant_signal_fp, var_chr_ids->at(i_chr));
			exit(0);
		}

		// Load the annotation signals.
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

		// Load the variant signals.
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
	sprintf(impact_string_fp, "%s/impact_value_strings.list", annotation_signal_dir);
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

// 
void translate_annotated_SNVs_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp)
{
	fprintf(stderr, "Interpreting the annotated deletions signals from %s to %s.\n", annotated_variant_signal_dir, op_fp);

	// Load the annotation signal.
	char anno_var_sig_chr_ids_fp[1000];
	sprintf(anno_var_sig_chr_ids_fp, "%s/chr_ids.txt", annotated_variant_signal_dir);
	vector<char*>* anno_var_chr_ids = buffer_file(anno_var_sig_chr_ids_fp);
	if (anno_var_chr_ids == NULL)
	{
		fprintf(stderr, "Could not load annotated variants chromosome id's from %s\n", anno_var_sig_chr_ids_fp);
		exit(0);
	}

	FILE* f_op = open_f(op_fp, "w");
	fprintf(f_op, "#CHROM\tDEL_START\tDEL_END\tDEL_NAME\tELEMENT_NAME\tELEMENT_STRAND\tREF_AA\tDEL_AA\tDEL_AA_alt\tIMPACT\n");

	for (int i_chr = 0; i_chr < anno_var_chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", anno_var_chr_ids->at(i_chr));

		// First load the EOI regions.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));

		if (!check_file(cur_sorted_BED_fp))
		{
			fprintf(stderr, "Skipping %s, no impact signals for this chromosome.\n", anno_var_chr_ids->at(i_chr));
		}

		vector<t_annot_region*>* sorted_cur_chr_EOI_regs = load_BED(cur_sorted_BED_fp);
		sort_set_sorting_info(sorted_cur_chr_EOI_regs, sort_regions_per_start_end_name);

		int cur_chr_EOI_covg = coverage(sorted_cur_chr_EOI_regs);
		fprintf(stderr, "Loaded %d EOI regions with %d nucleotide coverage on %s.\n",
			sorted_cur_chr_EOI_regs->size(),
			cur_chr_EOI_covg,
			anno_var_chr_ids->at(i_chr));

		// Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < sorted_cur_chr_EOI_regs->size(); i_reg++)
		{
			sorted_cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (sorted_cur_chr_EOI_regs->at(i_reg)->end - sorted_cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.
		fprintf(stderr, "Loaded %d EOI regions on %s\n", sorted_cur_chr_EOI_regs->size(), anno_var_chr_ids->at(i_chr));

		// Load the impact signals.
		fprintf(stderr, "Loading per allele impact signals.\n");
		unsigned long long** per_allele_SNV_impact_signal = new unsigned long long*[4];

		char cur_chr_anno_variant_signal_fp[1000];
		sprintf(cur_chr_anno_variant_signal_fp, "%s/annotated_variant_signal_%s.bin.gz", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_anno_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_anno_variant_signal_fp, anno_var_chr_ids->at(i_chr));
			exit(0);
		}

		FILE* f_annotation_signal = open_f(cur_chr_anno_variant_signal_fp, "rb");
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			fprintf(stderr, "Reading impacts per allele %d\n", i_allele);
			per_allele_SNV_impact_signal[i_allele] = new unsigned long long[cur_chr_EOI_covg + 2];
			if (fread(per_allele_SNV_impact_signal[i_allele], sizeof(unsigned long long), cur_chr_EOI_covg, f_annotation_signal) != cur_chr_EOI_covg)
			{
				fprintf(stderr, "Could not read allele_i=%d from %s\n", i_allele, cur_chr_anno_variant_signal_fp);
				exit(0);
			}
		} // i_allele loop.

		// Load the impact values from VEP's output.
		char sorted_impact_value_strings_list_fp[1000];
		sprintf(sorted_impact_value_strings_list_fp, "%s/impact_value_strings.list", annotated_variant_signal_dir);
		//vector<char*>* sorted_impact_vals = buffer_file(sorted_impact_value_strings_list_fp);
		vector<char*>* sorted_impact_vals = new vector<char*>();
		vector<t_VEP_term_ctx*>* vep_term_ctx = load_VEP_annotation_term_context(sorted_impact_value_strings_list_fp, sorted_impact_vals);
		fprintf(stderr, "Loaded %d impact strings.\n", sorted_impact_vals->size());

		// 1       9374375 rs926250        G       A       100     PASS    AC=3207;AF=0.640375;AN=5008;NS=2504;DP=13999;EAS_AF=0.2698;AMR_AF=0.6023;AFR_AF=0.8759;EUR_AF=0.6889;SAS_AF=0.681;AA=A|||;VT=SNP        GT      1|0
		char cur_chrom_vcf_fp[1000];
		sprintf(cur_chrom_vcf_fp, "%s/%s.vcf.gz", per_chrom_VCF_dir, anno_var_chr_ids->at(i_chr));
		if (check_file(cur_chrom_vcf_fp))
		{
			fprintf(stderr, "Reading deletes in %s\n", cur_chrom_vcf_fp);
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

		int n_bits_per_nuc = 3;
		int n_neigh_nucs = 6;

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

			t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");
			if (cur_line_toks->size() < (GENOTYPE_col_i + 1))
			{
				fprintf(stderr, "Expecting at least %d columns in VCF file.\n", GENOTYPE_col_i + 1);
				exit(0);
			}

			// Parse the location: Note that, for deletions, the position is set to be 1 plus the recorded position in the VCF file as that is the first deleted nucleotide.
			char* var_chr = t_string::copy_me_str(cur_line_toks->at(CHROM_col_i)->str());
			int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str());
			char* var_name = cur_line_toks->at(VAR_ID_col_i)->str();

			// Exclude SNVs.
			int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());
			char* ref_allele_seq = cur_line_toks->at(REF_ALLELE_col_i)->str();
			char* alt_allele_seq = cur_line_toks->at(ALT_ALLELE_col_i)->str();

			int eoi_reg_chr_i = t_string::get_i_str(anno_var_chr_ids, var_chr);
			if (eoi_reg_chr_i == anno_var_chr_ids->size() ||
				(l_alt_allele != 1 || l_ref_allele != 1))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// Generate the per allele counts from genotype value.
			char* geno_str = cur_line_toks->at(GENOTYPE_col_i)->str();
			//int per_allele_cnt[5];
			//memset(per_allele_cnt, 0, sizeof(int) * 5);
			//per_allele_cnt[nuc_2_num(ref_alt_allele[geno_str[0] - '0'])]++;
			//per_allele_cnt[nuc_2_num(ref_alt_allele[geno_str[2] - '0'])]++;

			// Analyze the genotype: Based on the specified genotype, assign the 
			if (geno_str[0] - '0' > 1 ||
				geno_str[2] - '0' > 1)
			{
				fprintf(stderr, "Invalid allele code: %s: %c, %c; %d, %d",
					cur_line,
					geno_str[0],
					geno_str[2],
					geno_str[0] - '0',
					geno_str[2] - '0');
				exit(0);
			}

			// Search the current variant.
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, sorted_cur_chr_EOI_regs, 0, sorted_cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_EOI_reg_i < sorted_cur_chr_EOI_regs->size() &&
				sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			bool found_overlap = false;
			while (cur_EOI_reg_i < sorted_cur_chr_EOI_regs->size() &&
				sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				int overlap_start = MAX(sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start, var_posn);
				int overlap_end = MIN(sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn);

				// Make sure the deletion is engulfed.
				if (overlap_start <= overlap_end &&
					overlap_start <= var_posn &&
					overlap_end >= var_posn)
				{
					fprintf(stderr, "----------------- Tracking SNV: %s:%d-%d (%s) -----------------\n",
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);


					if (__DUMP_CRYPTANNOT_SNV_MSGS__)
					{
						fprintf(stderr, "Adding SNV Overlap %s: %s:%d-%d:%s ; [%d-%d]\n", cur_line,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							overlap_start, overlap_end);
					}

					// For every overlap, update the allelic signal.
					vector<unsigned int>* cur_var_sorted_impacts = new vector<unsigned int>();

					// Copy the impact signal values, set the first nucleotide's coding frame, and the junction sequence.
					char* junction_cdna_seq = new char[n_neigh_nucs + 3];
					char* original_cdna_seq = new char[n_neigh_nucs + 3];
					memset(junction_cdna_seq, 0, (n_neigh_nucs + 3));
					memset(original_cdna_seq, 0, (n_neigh_nucs + 3));
					vector<int>* all_impacts_values = new vector<int>();
					int first_nuc_frame = 0;

					// Copy the annotation signal value.
					int track_posn = sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

					if (__DUMP_CRYPTANNOT_SNV_MSGS__)
					{
						fprintf(stderr, "Setting SNV Overlap @ %d (%d)\n",
							var_posn, track_posn);
						getc(stdin);
					}

					int alt_ref_nuc = nuc_2_num(alt_allele_seq[0]);
					unsigned long long cur_snv_annotation_signal_val = per_allele_SNV_impact_signal[alt_ref_nuc][track_posn];

					// Unpack and translate the indel.
					char s_neigh_seq_buff[100];
					unsigned int s_coding_frame = 0;
					unsigned long long s_impact_val = 0;
					unpack_impact_signal_values(cur_snv_annotation_signal_val,
												s_coding_frame, s_impact_val, s_neigh_seq_buff,
												n_bits_per_nuc, n_neigh_nucs);

					first_nuc_frame = s_coding_frame;
					unsigned long long cur_snv_allele_impact = s_impact_val;

					// Copy the original sequence.
					for (int nuc_i = 0; nuc_i < n_neigh_nucs / 2; nuc_i++)
					{
						junction_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
					} // nuc_i loop
					junction_cdna_seq[n_neigh_nucs / 2] = alt_allele_seq[0];
					for (int nuc_i = n_neigh_nucs / 2 + 1; nuc_i <= n_neigh_nucs; nuc_i++)
					{
						junction_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i - 1];
					} // nuc_i loop

					// Copy the original sequence.
					for (int nuc_i = 0; nuc_i < n_neigh_nucs / 2; nuc_i++)
					{
						original_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
					} // nuc_i loop
					original_cdna_seq[n_neigh_nucs / 2] = ref_allele_seq[0];
					for (int nuc_i = n_neigh_nucs / 2 + 1; nuc_i <= n_neigh_nucs; nuc_i++)
					{
						original_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i - 1];
					} // nuc_i loop

					// Concatenate the first and last sequences.
					// If the annotation is on the negative strand, reverse-complement the junction sequence.
					if (sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand == '-')
					{
						reverse_complement_seq(junction_cdna_seq);
						reverse_complement_seq(original_cdna_seq);
					}

					//////////////////////////////////////////////////////////////////////////////
					// We have setup the sequences and frames, interpret the impact.
					//////////////////////////////////////////////////////////////////////////////
					char* new_codon = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);
					char* prev_codon = (char*)(&original_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);
					char* new_AA = get_AA_code_per_codon(new_codon);
					char* prev_AA = get_AA_code_per_codon(prev_codon);

					fprintf(stderr, "%s::%s:%d-%d: %s->%s :: %s->%s @ Frame: %d (track posn: %d;; allele impact: %llX)\n",
						var_name,
						anno_var_chr_ids->at(i_chr),
						var_posn,
						var_posn,
						original_cdna_seq, junction_cdna_seq,
						prev_AA, new_AA,
						first_nuc_frame, 
						track_posn,
						cur_snv_allele_impact);

					// Depending on where the nucleotide is inserted, report the impact.
					t_string* vep_impact_str = new t_string();

					// Get the list of VEP strings and copy them.
					for (int term_i = 0; term_i < vep_term_ctx->size(); term_i++)
					{
						unsigned long long LLONE = 1;

						if (((LLONE << term_i) & cur_snv_allele_impact) > 0)
						{
							vep_impact_str->concat_string(vep_term_ctx->at(term_i)->id);
							vep_impact_str->concat_string(";");
						}
					} // term_i loop.s

					if (vep_impact_str->length() == 0)
					{
						vep_impact_str->concat_string("NOMINAL_IMPACT;");
					}

					// Write the translated annotation.
					fprintf(f_op, "%s\t%d\t%d\t%s\t%s\t%c\t%s\t%s\t%s\t%s\n",
						anno_var_chr_ids->at(i_chr),
						var_posn,
						var_posn,
						var_name,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand,
						prev_AA, new_AA, new_AA,
						vep_impact_str->str());

					if (__DUMP_CRYPTANNOT_SNV_MSGS__)
					{
						fprintf(stderr, "Sorted Impacts @ %s:%d-%d on %s (%c):\n",
								anno_var_chr_ids->at(i_chr),
								var_posn,
								var_posn,
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand);

						//for (int i_impact = 0; i_impact < cur_var_sorted_impacts->size(); i_impact++)
						{
							// Unpack and translate the indel.
							char s_neigh_seq_buff[100];
							unsigned int s_coding_frame = 0;
							unsigned long long s_impact_val = 0;
							unpack_impact_signal_values(cur_snv_annotation_signal_val,
														s_coding_frame, s_impact_val, s_neigh_seq_buff,
														n_bits_per_nuc, n_neigh_nucs);

							fprintf(stderr, "%llX -- %s @ frame %d\n",
									s_impact_val,
									s_neigh_seq_buff,
									s_coding_frame);
						} // i_impact loop.

						fprintf(stderr, "Junction cDNA Seq: %s\n", junction_cdna_seq);

						getc(stdin);
					} // messaging check.

					found_overlap = true;

					fprintf(stderr, "----------------- End New SNV Tracking -----------------\n");
				} // overlap check.

				cur_EOI_reg_i++;
			} // Forward searching loop.

			// FRee vcf memory.
			t_string::clean_tokens(cur_line_toks);
			delete[] cur_line;
		} // VCF file path.

		// Free memory.
		for (int i_allele = 0; i_allele < 4; i_allele++)
		{
			delete[] per_allele_SNV_impact_signal[i_allele];
		} // i_allele loop.
	} // i_chr loop.

	fclose(f_op);
}

void generate_EOI_randomized_SNV_VCF(char* EOI_bed_fp, double per_posn_var_prob, char* bin_seq_dir, char* op_vcf_fp)
{
	fprintf(stderr, "Simulating SNVs with per position probability of %.4f using sequence from %s.\n", per_posn_var_prob, bin_seq_dir);

	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	vector<t_annot_region*>* eoi_regs = load_BED(EOI_bed_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", eoi_regs->size());
	t_restr_annot_region_list* restr_eoi_regs = restructure_annot_regions(eoi_regs);

	FILE* f_op_vcf = open_f(op_vcf_fp, "w");
	for (int i_chr = 0; i_chr < restr_eoi_regs->chr_ids->size(); i_chr++)
	{
		// Load the sequence of the chromosome.
		fprintf(stderr, "Loading sequence of %s\n", restr_eoi_regs->chr_ids->at(i_chr));
		char cur_chr_seq_fp[1000];
		sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, restr_eoi_regs->chr_ids->at(i_chr));
		if (!check_file(cur_chr_seq_fp))
		{
			sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, restr_eoi_regs->chr_ids->at(i_chr));

			if (!check_file(cur_chr_seq_fp))
			{
				fprintf(stderr, "Could not find %s.\n", cur_chr_seq_fp);
				continue;
			}
		}

		int l_chrom_seq = 0;
		char* cur_chr_seq = load_binary_sequence_file(cur_chr_seq_fp, l_chrom_seq);
		fprintf(stderr, "Loaded %d nucleotides.\n", l_chrom_seq);

		// Start writing the SNV VCF.
		char nucs[] = "ACGT";

		char* reference_allele = new char[5];
		char* alternate_allele = new char[5];

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_eoi_regs->regions_per_chrom[i_chr];
		fprintf(stderr, "Writing the SNVs for %d regions.\n", (int)(cur_chr_EOI_regs->size()));
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			if (i_reg % 500 == 0)
			{
				fprintf(stderr, "Processing %d. region.            \r", i_reg);
			}

			for (int i = cur_chr_EOI_regs->at(i_reg)->start; i <= cur_chr_EOI_regs->at(i_reg)->end; i++)
			{
				if (rng->random_double_ran3() > per_posn_var_prob)
				{
					continue;
				}

				char alt_allele = cur_chr_seq[i];
				while (nuc_2_num(alt_allele) == nuc_2_num(cur_chr_seq[i]))
				{
					alt_allele = num_2_nuc(rng->random_double_ran3() * 4);
				}

				alternate_allele[0] = alt_allele;
				alternate_allele[1] = 0;

				reference_allele[0] = cur_chr_seq[i];
				reference_allele[1] = 0;

				// #CHROM POS     ID        REF ALT    QUAL FILTER
				fprintf(f_op_vcf, "%s\t%d\t%s_%d_SNV_%s\t%s\t%s\t60\tPASS\t.\tGT\t0|1\n",
					restr_eoi_regs->chr_ids->at(i_chr),
					translate_coord(i, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
					restr_eoi_regs->chr_ids->at(i_chr), i, cur_chr_EOI_regs->at(i_reg)->name,
					reference_allele,
					alternate_allele);
			} // i loop.
		} // i_reg loop.

		delete[] cur_chr_seq;
	} // i_chr loop.
	close_f(f_op_vcf, op_vcf_fp);
}


// This generates the indicateor signal.
void generate_EOI_enumerating_SNV_VCF(char* EOI_bed_fp, char* bin_seq_dir, char* op_vcf_fp)
{
	vector<t_annot_region*>* eoi_regs = load_BED(EOI_bed_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", eoi_regs->size());
	t_restr_annot_region_list* restr_eoi_regs = restructure_annot_regions(eoi_regs);

	FILE* f_op_vcf = open_f(op_vcf_fp, "w");
	for (int i_chr = 0; i_chr < restr_eoi_regs->chr_ids->size(); i_chr++)
	{
		// Load the sequence of the chromosome.
		fprintf(stderr, "Loading sequence of %s\n", restr_eoi_regs->chr_ids->at(i_chr));
		char cur_chr_seq_fp[1000];
		sprintf(cur_chr_seq_fp, "%s/%s.bin", bin_seq_dir, restr_eoi_regs->chr_ids->at(i_chr));
		if (!check_file(cur_chr_seq_fp))
		{
			sprintf(cur_chr_seq_fp, "%s/%s.bin.gz", bin_seq_dir, restr_eoi_regs->chr_ids->at(i_chr));

			if (!check_file(cur_chr_seq_fp))
			{
				fprintf(stderr, "Could not find %s.\n", cur_chr_seq_fp);
				continue;
			}
		}

		int l_chrom_seq = 0;
		char* cur_chr_seq = load_binary_sequence_file(cur_chr_seq_fp, l_chrom_seq);
		fprintf(stderr, "Loaded %d nucleotides.\n", l_chrom_seq);

		// Start writing the SNV VCF.
		char nucs[] = "ACGT";

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_eoi_regs->regions_per_chrom[i_chr];
		fprintf(stderr, "Writing the SNVs for %d regions.\n", (int)(cur_chr_EOI_regs->size()));
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			for (int i = cur_chr_EOI_regs->at(i_reg)->start; i <= cur_chr_EOI_regs->at(i_reg)->end; i++)
			{
				for (int nuc_i = 0; nuc_i < 4; nuc_i++)
				{
					if (toupper(nucs[nuc_i]) != toupper(cur_chr_seq[i]))
					{
						// #CHROM POS     ID        REF ALT    QUAL FILTER 
						fprintf(f_op_vcf, "%s\t%d\t%s_%d_%c_%s\t%c\t%c\t60\tPASS\n",
							restr_eoi_regs->chr_ids->at(i_chr),
							translate_coord(i, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
							restr_eoi_regs->chr_ids->at(i_chr), i, nucs[nuc_i], cur_chr_EOI_regs->at(i_reg)->name,
							cur_chr_seq[i], nucs[nuc_i]);
					}
				} // nuc_i loop.
			} // i loop.
		} // i_reg loop.

		delete[] cur_chr_seq;
	} // i_chr loop.
	close_f(f_op_vcf, op_vcf_fp);
}

