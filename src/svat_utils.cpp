#include <stdio.h>
#include <stdlib.h>
#include "file_utils.h"
#include "svat_ansi_string.h"
#include "svat_annot_region_tools.h"
#include "svat_genomics_coords.h"
#include "svat_nomenclature.h"
#include "svat_genome_sequence_tools.h"
#include "svat_nucleotide.h"
#include "svat_variation_tools.h"
#include "svat_utils.h"
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <algorithm>
#include <vector>
using namespace std;

bool __DUMP_CRYPTANNOT_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))


bool check_element_context_per_impact_bitmap(vector<t_VEP_term_ctx*>* vep_term_ctx, int context_id, unsigned long long impact_entry_bitmap)
{
	unsigned long long val = 0;
	unsigned long long LLONE = 1;
	for (int vep_term_i = 0; vep_term_i < vep_term_ctx->size(); vep_term_i++)
	{
		// Check if this term is included in the current 
		if ((impact_entry_bitmap & (LLONE << vep_term_i)) > 0)
		{
			if (vep_term_ctx->at(vep_term_i)->context_array[context_id])
			{
				return(true);
			}
		}
	} // vep_term_i loop.

	return(false);
}

void unpack_neigh_nucs(unsigned int neigh_seq_signal, char* buff, int n_bits_per_nuc, int n_nucs)
{
	unsigned long long LLONE = 1;
	unsigned int mask = (1 << n_bits_per_nuc) - 1;
	for (int nuc_i = n_nucs - 1; nuc_i >= 0; nuc_i--)
	{
		int nuc_value = (neigh_seq_signal >> (nuc_i * n_bits_per_nuc)) & mask;
		buff[nuc_i] = num_2_nuc(nuc_value);
	} // nuc_i loop.

	buff[n_nucs] = 0;
}

unsigned long long pack_set_impact_signal_values_w_CDS_neighborhood(unsigned int coding_frame,
	unsigned long long impact_value,
	unsigned int geno_neigh_seq_signal,
	unsigned int cds_neigh_seq_signal,
	int n_bits_per_nuc, int n_neigh_nucs)
{
	unsigned long long LLONE = 1;
	if (coding_frame >= (LLONE << 2) ||
		geno_neigh_seq_signal >= (LLONE << (n_bits_per_nuc*n_neigh_nucs)) ||
		cds_neigh_seq_signal >= (LLONE << (n_bits_per_nuc*n_neigh_nucs)))
	{
		fprintf(stderr, "The signal values are not within range: %u, %u, %u\n", coding_frame, impact_value, geno_neigh_seq_signal);
		exit(0);
	}

	unsigned long long packed_impact_signal_value = 0;
	packed_impact_signal_value = (impact_value);
	packed_impact_signal_value = (packed_impact_signal_value << (n_bits_per_nuc*n_neigh_nucs)) | cds_neigh_seq_signal;
	packed_impact_signal_value = (packed_impact_signal_value << (n_bits_per_nuc*n_neigh_nucs)) | geno_neigh_seq_signal;
	packed_impact_signal_value = (packed_impact_signal_value << 2) | coding_frame;

	return(packed_impact_signal_value);
}

unsigned long long pack_set_impact_signal_values(unsigned int coding_frame,
	unsigned long long impact_value,
	unsigned int neigh_seq_signal, 
	int n_bits_per_nuc, int n_neigh_nucs)
{
	unsigned long long LLONE = 1;
	if (coding_frame >= (LLONE << 2) ||
		neigh_seq_signal >= (LLONE << (n_bits_per_nuc*n_neigh_nucs)))
	{
		fprintf(stderr, "The signal values are not within range: %u, %u, %u\n", coding_frame, impact_value, neigh_seq_signal);
		exit(0);
	}

	unsigned long long packed_impact_signal_value = 0;
	packed_impact_signal_value = (impact_value);
	packed_impact_signal_value = (packed_impact_signal_value << (n_bits_per_nuc*n_neigh_nucs)) | neigh_seq_signal;
	packed_impact_signal_value = (packed_impact_signal_value << 2) | coding_frame;

	return(packed_impact_signal_value);
}

void unpack_impact_signal_values_w_CDS_neighborhood(unsigned long long packed_impact_signal_value,
	unsigned int& coding_frame,
	unsigned long long& impact_value,
	char* geno_neigh_seq_nucs, 
	char* CDS_neigh_seq_nucs, int n_bits_per_nuc, int n_neigh_nucs)
{
	// Read the coding frame.
	coding_frame = packed_impact_signal_value & 3;

	packed_impact_signal_value = packed_impact_signal_value >> 2;

	unsigned int neigh_seq_mask = ((1 << (n_bits_per_nuc * n_neigh_nucs)) - 1);

	// Read the genomic-neighborhood.
	unsigned int geno_neigh_seq_signal = (packed_impact_signal_value & neigh_seq_mask);
	unpack_neigh_nucs(geno_neigh_seq_signal, geno_neigh_seq_nucs, n_bits_per_nuc, n_neigh_nucs);

	packed_impact_signal_value = packed_impact_signal_value >> (n_bits_per_nuc*n_neigh_nucs);

	// Read the CDS-neighborhood.
	unsigned int CDS_neigh_seq_signal = (packed_impact_signal_value & neigh_seq_mask);
	unpack_neigh_nucs(geno_neigh_seq_signal, geno_neigh_seq_nucs, n_bits_per_nuc, n_neigh_nucs);

	packed_impact_signal_value = packed_impact_signal_value >> (n_bits_per_nuc*n_neigh_nucs);
	impact_value = packed_impact_signal_value & 0xffffffffffffffff;
}

void unpack_impact_signal_values(unsigned long long packed_impact_signal_value,
	unsigned int& coding_frame,
	unsigned long long& impact_value,
	char* neigh_seq_nucs, int n_bits_per_nuc, int n_neigh_nucs)
{
	coding_frame = packed_impact_signal_value & 3;
	packed_impact_signal_value = packed_impact_signal_value >> 2;

	unsigned int neigh_seq_mask = ((1 << (n_bits_per_nuc * n_neigh_nucs)) - 1);
	unsigned int neigh_seq_signal = (packed_impact_signal_value & neigh_seq_mask);
	unpack_neigh_nucs(neigh_seq_signal, neigh_seq_nucs, n_bits_per_nuc, n_neigh_nucs);

	packed_impact_signal_value = packed_impact_signal_value >> (n_bits_per_nuc*n_neigh_nucs);
	impact_value = packed_impact_signal_value & 0xffffffffffffffff;
}

vector<t_VEP_term_ctx*>* load_VEP_annotation_term_context(char* vep_annotation_term_ctx_fp, vector<char*>* sorted_vep_terms)
{
	if (!check_file(vep_annotation_term_ctx_fp))
	{
		fprintf(stderr, "Could not find VEP annotation term file @ %s\n", vep_annotation_term_ctx_fp);
		exit(0);
	}
	
	vector<t_VEP_term_ctx*>* vep_annotation_term_ctx = new vector<t_VEP_term_ctx*>();

	FILE* f_vep_annotation_terms_ctx = open_f(vep_annotation_term_ctx_fp, "r");

	// Load the header.
	char* header_line = getline(f_vep_annotation_terms_ctx);
	vector<char*>* tokens = t_string::copy_tokens_2_strs(t_string::tokenize_by_chars(header_line, "\t"));
	int cds_col_i = t_string::get_i_str_per_substr_search(tokens, "CDS");
	int splice_region_col_i = t_string::get_i_str_per_substr_search(tokens, "Splice_Region");
	int splice_acceptor_col_i = t_string::get_i_str_per_substr_search(tokens, "Splice_Acceptor");
	int splice_donor_col_i = t_string::get_i_str_per_substr_search(tokens, "Splice_Donor");
	int intron_col_i = t_string::get_i_str_per_substr_search(tokens, "Intron");
	int fp_UTR_col_i = t_string::get_i_str_per_substr_search(tokens, "5p_UTR");
	int tp_UTR_col_i = t_string::get_i_str_per_substr_search(tokens, "3p_UTR");
	int start_col_i = t_string::get_i_str_per_substr_search(tokens, "Start_Codon");
	int stop_col_i = t_string::get_i_str_per_substr_search(tokens, "Stop_Codon");

	if (cds_col_i == tokens->size() || 
		splice_region_col_i == tokens->size() ||
		splice_donor_col_i == tokens->size() ||
		splice_acceptor_col_i == tokens->size() ||
		intron_col_i == tokens->size() ||
		fp_UTR_col_i == tokens->size() ||
		tp_UTR_col_i == tokens->size() ||
		start_col_i == tokens->size() ||
		stop_col_i == tokens->size())
	{
		fprintf(stderr, "Could not find CDS/splicing/start_stop context columns: %s\n", header_line);
		exit(0);
	}

	while (1)
	{
		char* cur_line = getline(f_vep_annotation_terms_ctx);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* toks = t_string::tokenize_by_chars(cur_line, "\t");
		if (toks->size() < cds_col_i ||
			toks->size() < start_col_i ||
			toks->size() < stop_col_i ||
			toks->size() < splice_region_col_i ||
			toks->size() < splice_donor_col_i ||
			toks->size() < splice_acceptor_col_i || 
			toks->size() < fp_UTR_col_i || 
			toks->size() < tp_UTR_col_i || 
			toks->size() < intron_col_i)
		{
			fprintf(stderr, "Could not find enough columns in %s while loading VEP annotation terms.\n", cur_line);
			exit(0);
		}

		t_VEP_term_ctx* cur_term_ctx = new t_VEP_term_ctx();
		cur_term_ctx->id = t_string::copy_me_str(toks->at(0)->str());
		cur_term_ctx->intron_ctx = toks->at(intron_col_i)->str()[0] == '1';
		cur_term_ctx->CDS_ctx = toks->at(cds_col_i)->str()[0] == '1';		
		cur_term_ctx->start_codon_ctx = toks->at(start_col_i)->str()[0] == '1';
		cur_term_ctx->stop_codon_ctx = toks->at(stop_col_i)->str()[0] == '1';
		cur_term_ctx->tp_UTR_ctx = toks->at(tp_UTR_col_i)->str()[0] == '1';
		cur_term_ctx->fp_UTR_ctx = toks->at(fp_UTR_col_i)->str()[0] == '1';
		cur_term_ctx->splice_donor_ctx = toks->at(splice_donor_col_i)->str()[0] == '1';
		cur_term_ctx->splice_acceptor_ctx = toks->at(splice_acceptor_col_i)->str()[0] == '1';
		cur_term_ctx->splice_region_ctx = toks->at(splice_region_col_i)->str()[0] == '1';

		cur_term_ctx->context_array = new bool[N_VEP_CTX];
		cur_term_ctx->context_array[intron_ctx] = toks->at(intron_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[CDS_ctx] = toks->at(cds_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[start_codon_ctx] = toks->at(start_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[stop_codon_ctx] = toks->at(stop_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[tp_UTR_ctx] = toks->at(tp_UTR_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[fp_UTR_ctx] = toks->at(fp_UTR_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[splice_donor_ctx] = toks->at(splice_donor_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[splice_acceptor_ctx] = toks->at(splice_acceptor_col_i)->str()[0] == '1';
		cur_term_ctx->context_array[splice_region_ctx] = toks->at(splice_region_col_i)->str()[0] == '1';

		sorted_vep_terms->push_back(t_string::copy_me_str(cur_term_ctx->id));

		vep_annotation_term_ctx->push_back(cur_term_ctx);
	} // file reading loop.
	close_f(f_vep_annotation_terms_ctx, vep_annotation_term_ctx_fp);

	return(vep_annotation_term_ctx);
}

vector<t_annot_region*>* load_VEP_as_regions(char* vep_op_fp)
{
	vector<t_annot_region*>* vep_regs = new vector<t_annot_region*>();

	FILE* f_vep_op = open_f(vep_op_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_vep_op);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			continue;
		}

		char cur_var_id[1000];
		sscanf(cur_line, "%s", cur_var_id);

		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str("XX");
		reg->start = 1;
		reg->end = 1;
		reg->strand = '+';
		reg->name = t_string::copy_me_str(cur_var_id);
		reg->data = cur_line;

		vep_regs->push_back(reg);
	}
	fclose(f_vep_op);

	return(vep_regs);
}

void compare_VEP_annotations_with_SVAT_annotations(char* vep_annotated_op_fp, 
													char* VEP_impact_string_ctx_fp,
													char* vep_impacts_2_focus_list_fp,
													char* high_impacts_list_fp,
													char* svat_annotated_vcf_fp)
{
	vector<char*>* sorted_impact_vals = new vector<char*>();
	vector<t_VEP_term_ctx*>* vep_term_ctx = load_VEP_annotation_term_context(VEP_impact_string_ctx_fp, sorted_impact_vals);
	fprintf(stderr, "Loaded %d VEP impact strings.\n", sorted_impact_vals->size());

	vector<t_annot_region*>* svat_regs = load_BED_with_line_information(svat_annotated_vcf_fp);
	fprintf(stderr, "Loaded %d SVAT annotated regions.\n", svat_regs->size());
	
	vector<t_annot_region*>* vep_regs = load_VEP_as_regions(vep_annotated_op_fp);
	fprintf(stderr, "Loaded %d VEP annotated regions.\n", vep_regs->size());

	vector<char*>* vep_impacts_2_focus = buffer_file(vep_impacts_2_focus_list_fp);
	fprintf(stderr, "Focusing on %d vep impacts.\n", vep_impacts_2_focus->size());

	vector<char*>* high_impacts = buffer_file(high_impacts_list_fp);
	fprintf(stderr, "Evaluating %d high impacts.\n", high_impacts->size());

	for (int i_reg = 0; i_reg < vep_regs->size(); i_reg++)
	{
		vep_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d VEP annotations, %d signal annotations.\n", vep_regs->size(), svat_regs->size());

	// Intersect the annotated regions.
	vector<t_annot_region*>* intersects = intersect_regions_per_names(svat_regs, vep_regs, true);
	FILE* f_vep_specific_annotations = open_f("vep_specific_annotations.txt", "w");
	FILE* f_svat_specific_annotations = open_f("svat_specific_annotations.txt", "w");

	FILE* f_vep_specific_high_impact_annotations = open_f("vep_specific_HI_annotations.txt", "w");
	FILE* f_svat_specific_high_impact_annotations = open_f("svat_specific_HI_annotations.txt", "w");

	FILE* f_matching_svat = open_f("matched_svat_annotations.txt", "w");
	FILE* f_matching_vep = open_f("matched_vep_annotations.txt", "w");

	FILE* f_non_existing_annotations = open_f("non_existing_annotations.txt", "w");

	FILE* f_any_HI_mismatches = open_f("any_HI_mismatches.txt", "w");

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "Processing %d/%d intersect.                \r", i_int, intersects->size());
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* svat_reg = int_info->src_reg;
		t_annot_region* vep_reg = int_info->dest_reg;		

		// Parse the VEP line.
		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		int SVAT_IMPACT_col_i = 9;
		int SVAT_element_col_i = 4;

		// Lock to transcript specific analysis.
		int EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		int EOI_element_VEP_col_i = -1;
		if (EOI_element_type == VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION)
		{
			EOI_element_VEP_col_i = TRANSCRIPT_ID_col_i;
		}
		else if (EOI_element_type == VEP_SIGNALIZE_GENE_SUMMARIZATION)
		{
			EOI_element_VEP_col_i = GENE_ID_col_i;
		}

		char* cur_vep_line = (char*)(vep_reg->data);
		if (cur_vep_line == NULL)
		{
			break;
		}

		t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_vep_line, "\t");

		char* vep_element_id = cur_line_toks->at(EOI_element_VEP_col_i)->str();

		// Go over the consequence tokens and choose the highest impact.
		t_string_tokens* VEP_conseq_toks = cur_line_toks->at(CONSEQUENCE_col_i)->tokenize_by_chars(", ");
		vector<char*>* VEP_conseq_strs = t_string::copy_tokens_2_strs(VEP_conseq_toks);
		t_string::clean_tokens(VEP_conseq_toks);

		// Extract the SVAT annotated impacts.
		char* cur_svat_line = (char*)(svat_reg->data);

		// 1       11479173        11479176        1_11479172_DEL_ENST00000294484  ENST00000423056 +       0       1       -       -       IMPACT_NOMINAL
		t_string_tokens* svat_line_toks = t_string::tokenize_by_chars(cur_svat_line, "\t");
		t_string_tokens* svat_impact_toks = svat_line_toks->at(SVAT_IMPACT_col_i)->tokenize_by_chars(";");
		vector<char*>* svat_impact_strs = t_string::copy_tokens_2_strs(svat_impact_toks);
		t_string::clean_tokens(svat_impact_toks);

		// The impacted elements must match.
		if (t_string::compare_strings(vep_element_id, svat_line_toks->at(SVAT_element_col_i)->str()))
		{
			vep_reg->score = 1;

			// Pairwise compare.
			bool found_vep_focus_term = false;
			bool cur_VEP_HI = false;
			for (int i_tok = 0; i_tok < VEP_conseq_strs->size(); i_tok++)
			{
				if (t_string::get_i_str(high_impacts, VEP_conseq_strs->at(i_tok)) < high_impacts->size())
				{
					cur_VEP_HI = true;
				}

				// Make sure that this is an impact we want to focus on.
				if (t_string::get_i_str(vep_impacts_2_focus, VEP_conseq_strs->at(i_tok)) < vep_impacts_2_focus->size())
				{
					found_vep_focus_term = true;

					bool found_cur_impact = 0;
					if (t_string::get_i_str(svat_impact_strs, VEP_conseq_strs->at(i_tok)) < svat_impact_strs->size())
					{
						found_cur_impact = true;
					}

					if (!found_cur_impact)
					{
						fprintf(f_vep_specific_annotations, "%s\t%s\t%s\t%s\t%s\n",
							VEP_conseq_strs->at(i_tok),
							cur_line_toks->at(CONSEQUENCE_col_i)->str(),
							svat_line_toks->at(SVAT_IMPACT_col_i)->str(),
							cur_vep_line,
							cur_svat_line);

						if (t_string::get_i_str(high_impacts, VEP_conseq_strs->at(i_tok)) < high_impacts->size())
						{
							fprintf(f_vep_specific_high_impact_annotations, "%s\t%s\t%s\t%s\t%s\n",
								VEP_conseq_strs->at(i_tok),
								cur_line_toks->at(CONSEQUENCE_col_i)->str(),
								svat_line_toks->at(SVAT_IMPACT_col_i)->str(),
								cur_vep_line,
								cur_svat_line);
						}
					}
				}
			} // i_tok loop.

			if (found_vep_focus_term)
			{
				fprintf(f_matching_vep, "%s\t%s\n", cur_vep_line, cur_svat_line);
			}

			// Pairwise compare.
			bool found_svat_focus_term = false;
			bool cur_svat_HI = false;
			for (int i_tok = 0; i_tok < svat_impact_strs->size(); i_tok++)
			{
				if (t_string::get_i_str(high_impacts, svat_impact_strs->at(i_tok)) < high_impacts->size())
				{
					cur_svat_HI = true;
				}

				// Make sure that this is an impact we want to focus on.
				if (t_string::get_i_str(vep_impacts_2_focus, svat_impact_strs->at(i_tok)) < vep_impacts_2_focus->size())
				{
					found_svat_focus_term = true;

					bool found_cur_impact = 0;
					if (t_string::get_i_str(VEP_conseq_strs, svat_impact_strs->at(i_tok)) < VEP_conseq_strs->size())
					{
						found_cur_impact = true;
					}

					if (!found_cur_impact)
					{
						fprintf(f_svat_specific_annotations, "%s\t%s\t%s\t%s\t%s\n",
								svat_impact_strs->at(i_tok), 
								svat_line_toks->at(SVAT_IMPACT_col_i)->str(),
								cur_line_toks->at(CONSEQUENCE_col_i)->str(),
								cur_vep_line, 
								cur_svat_line);

						if (t_string::get_i_str(high_impacts, svat_impact_strs->at(i_tok)) < high_impacts->size())
						{
							fprintf(f_svat_specific_high_impact_annotations, "%s\t%s\t%s\t%s\t%s\n",
								svat_impact_strs->at(i_tok),
								svat_line_toks->at(SVAT_IMPACT_col_i)->str(),
								cur_line_toks->at(CONSEQUENCE_col_i)->str(),
								cur_vep_line,
								cur_svat_line);
						}
					}
				}
			} // i_tok loop.

			if (cur_svat_HI != cur_VEP_HI)
			{
				fprintf(f_any_HI_mismatches, "%d\t%d\t%s\t%s\n", cur_svat_HI, cur_VEP_HI, cur_vep_line, cur_svat_line);
			}

			if (found_svat_focus_term)
			{
				fprintf(f_matching_svat, "%s\t%s\n", cur_vep_line, cur_svat_line);
			}
		} // element_id comparison.

		// Clean tokens.
		t_string::clean_tokens(cur_line_toks);
		t_string::clean_string_list(VEP_conseq_strs);
		t_string::clean_tokens(svat_line_toks);
		t_string::clean_string_list(svat_impact_strs);
	} // i_int loop.

	fclose(f_vep_specific_annotations);
	fclose(f_svat_specific_annotations);

	fclose(f_vep_specific_high_impact_annotations);
	fclose(f_svat_specific_high_impact_annotations);

	fclose(f_any_HI_mismatches);

	// Write the unmatched VEP annotations.
	for (int i_reg = 0; i_reg < vep_regs->size(); i_reg++)
	{
		if (vep_regs->at(i_reg)->score == 0)
		{
			char* cur_line = (char*)(vep_regs->at(i_reg)->data);
			fprintf(f_non_existing_annotations, "NOT_MATCHED::%s\n", cur_line);
		}
	} // i_reg loop.

	fclose(f_non_existing_annotations);

	fclose(f_matching_vep);
	fclose(f_matching_svat);
}

void separate_VCF_per_chromosome(char* vcf_fp, char* chr_ids_lengths_list_fp, char* per_chrom_op_dir)
{
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_list_fp, chr_ids, chr_lengths);
	fprintf(stderr, "Separating %s into %d chromosomes.\n", vcf_fp, chr_ids->size());

	vector<FILE*>* f_ptr_per_chrom = new vector<FILE*>();
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chrom_op_fp[1000];
		sprintf(cur_chrom_op_fp, "%s/%s.vcf.gz", per_chrom_op_dir, chr_ids->at(i_chr));
		FILE* f_chrom = open_f(cur_chrom_op_fp, "w");

		f_ptr_per_chrom->push_back(f_chrom);
	} // i_chr loop.

	int n_processed_lines = 0;
	FILE* f_vcf = open_f(vcf_fp, "r");
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

		n_processed_lines++;

		if (n_processed_lines % 10000 == 0)
		{
			fprintf(stderr, "@ %d. line           \r", n_processed_lines);
		}

		char cur_chrom_id[100];
		sscanf(cur_line, "%s", cur_chrom_id);
		normalize_chr_id(cur_chrom_id);
		int i_next_char = 0;

		int chr_i = t_string::get_i_str(chr_ids, cur_chrom_id);
		if (chr_i < chr_ids->size())
		{
			fprintf(f_ptr_per_chrom->at(chr_i), "%s\n", cur_line);
		}

		delete[] cur_line;
	} // vep output file reading loop.
	fclose(f_vcf);

	// Close the files.
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chrom_op_fp[1000];
		sprintf(cur_chrom_op_fp, "%s/%s.vcf.gz", per_chrom_op_dir, chr_ids->at(i_chr));
		close_f(f_ptr_per_chrom->at(i_chr), cur_chrom_op_fp);
	} // i_chr loop.
}

void separate_VEP_output_per_chromosome(char* vep_op_fp, char* chr_ids_lengths_list_fp, char* per_chrom_op_dir)
{
	vector<char*>* chr_ids = new vector<char*>();
	vector<int>* chr_lengths = new vector<int>();
	load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_list_fp, chr_ids, chr_lengths);
	fprintf(stderr, "Separating %s into %d chromosomes.\n", vep_op_fp, chr_ids->size());

	vector<FILE*>* f_ptr_per_chrom = new vector<FILE*>();
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chrom_op_fp[1000];
		sprintf(cur_chrom_op_fp, "%s/%s.vep.gz", per_chrom_op_dir, chr_ids->at(i_chr));
		FILE* f_chrom = open_f(cur_chrom_op_fp, "w");

		f_ptr_per_chrom->push_back(f_chrom);
	} // i_chr loop.

	int n_processed_lines = 0;
	FILE* f_vep = open_f(vep_op_fp, "r");
	while (1)
	{
		char* cur_line = getline(f_vep);
		if (cur_line == NULL)
		{
			break;
		}

		n_processed_lines++;

		if (n_processed_lines % 10000 == 0)
		{
			fprintf(stderr, "@ %d. line           \r", n_processed_lines);
		}

		char cur_location_str[100];
		sscanf(cur_line, "%*s %s", cur_location_str);
		int i_next_char = 0;
		t_string_tokens* toks = t_string::get_first_n_tokens(cur_location_str, 1, ":", i_next_char);
		char* chrom_id = toks->at(0)->str();

		int chr_i = t_string::get_i_str(chr_ids, chrom_id);
		if (chr_i < chr_ids->size())
		{
			fprintf(f_ptr_per_chrom->at(chr_i), "%s\n", cur_line);
		}

		// free memory.
		delete[] cur_line;
	} // vep output file reading loop.
	fclose(f_vep);

	// Close the files.
	for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
	{
		char cur_chrom_op_fp[1000];
		sprintf(cur_chrom_op_fp, "%s/%s.vep.gz", per_chrom_op_dir, chr_ids->at(i_chr));
		close_f(f_ptr_per_chrom->at(i_chr), cur_chrom_op_fp);
	} // i_chr loop.
}

bool sort_gene_transcript_ids_per_gene_id(char** ids1, char** ids2)
{
	return(t_string::sort_strings_per_prefix(ids1[0], ids2[0]));
}

void get_frame_gene_unique_EOIs_BED_per_CDS_regions(char* gene_ids_list_fp, char* CDS_exon_selector, char* GENCODE_gff_fp, int l_exon_ext, char* op_fp)
{
	if (!t_string::compare_strings(CDS_exon_selector, "CDS") &&
		!t_string::compare_strings(CDS_exon_selector, "exon"))
	{
		fprintf(stderr, "Need CDS/exon as the subelement selector.\n");
		exit(0);
	}

	// Read GFF file, load genes/transcripts and write exons for all the transcripts and genes.
	FILE* f_gff = open_f(GENCODE_gff_fp, "r");

	fprintf(stderr, "Generating exonic EOI regions using %s using \"%s\" as subelements\n", GENCODE_gff_fp, CDS_exon_selector);

	vector<char*>* gene_ids = buffer_file(gene_ids_list_fp);
	fprintf(stderr, "Loaded %d genes, setting up the exon region lists.\n", gene_ids->size());

	//vector<char**>* gene_transcript_ids = new vector<char**>();
	//for (int i_el = 0; i_el < gene_transcript_id_lines->size(); i_el++)
	//{
	//	char cur_gene_id[100];
	//	char cur_trans_id[100];
	//	if (sscanf(gene_transcript_id_lines->at(i_el), "%s %s", cur_gene_id, cur_trans_id) != 2)
	//	{
	//		fprintf(stderr, "Could not parse %s\n", gene_transcript_id_lines->at(i_el));
	//		exit(0);
	//	}

	//	char** cur_gene_transcript_ids = new char*[2];
	//	cur_gene_transcript_ids[0] = t_string::copy_me_str(cur_gene_id);
	//	cur_gene_transcript_ids[1] = t_string::copy_me_str(cur_trans_id);

	//	gene_transcript_ids->push_back(cur_gene_transcript_ids);
	//} // i_el loop.

	// Sort element id's.
	//sort(gene_transcript_ids->begin(), gene_transcript_ids->end(), sort_gene_transcript_ids_per_gene_id);
	sort(gene_ids->begin(), gene_ids->end(), t_string::sort_strings_per_prefix);

	//vector<char*>* gene_sorted_gene_ids = new vector<char*>();
	//vector<char*>* gene_sorted_trans_ids = new vector<char*>();

	vector<vector<t_annot_region*>*>* per_gene_exons = new vector<vector<t_annot_region*>*>();
	for (int i_el = 0; i_el < gene_ids->size(); i_el++)
	{
		per_gene_exons->push_back(new vector<t_annot_region*>());
	} // i_el loop.

	  // GFF file reading.
	int n_entries = 0;
	while (1)
	{
		char* cur_line = getline(f_gff);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if ((n_entries > 0) &&
			(n_entries % 10000) == 0)
		{
			fprintf(stderr, "Processing %d. entry: %d elements.            \r", n_entries, (int)(gene_ids->size()));
		}

		n_entries++;

		// If this is an exonic entry, extract its element id.
		int i_cur_char = 0;
		char strand_char = '+';
		char region_name[1000];
		strcpy(region_name, ".");
		int region_score = 0;

		char chrom[100];
		t_string::get_next_token(cur_line, chrom, 100, "\t", i_cur_char);
		char source[100];
		t_string::get_next_token(cur_line, source, 100, "\t", i_cur_char);
		char element_type[100];
		t_string::get_next_token(cur_line, element_type, 100, "\t", i_cur_char);

		if (t_string::compare_strings(element_type, CDS_exon_selector))
		{
			char start_str[100];
			char end_str[100];
			char score_str[100];
			char strand_str[100];
			char frame_str[100];
			char* attribute = new char[strlen(cur_line) + 2];
			t_string::get_next_token(cur_line, start_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, end_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, score_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, strand_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, frame_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, attribute, strlen(cur_line), "\t", i_cur_char);

			int cur_element_frame = atoi(frame_str);

			t_string_tokens* attribute_toks = t_string::tokenize_by_chars(attribute, "; =");
			char* gene_id = NULL;
			char* transcript_id = NULL;
			for (int i_t = 0; i_t + 1 < attribute_toks->size(); i_t++)
			{
				if (t_string::compare_strings(attribute_toks->at(i_t)->str(), "gene_id"))
				{
					gene_id = t_string::copy_me_str(attribute_toks->at(i_t + 1)->str());
				}

				if (t_string::compare_strings(attribute_toks->at(i_t)->str(), "transcript_id"))
				{
					transcript_id = t_string::copy_me_str(attribute_toks->at(i_t + 1)->str());
				}
			} // i_t loop.

			if (transcript_id == NULL || gene_id == NULL)
			{
				fprintf(stderr, "Could not parse gene/transcript id's: %s\n", cur_line);
				exit(0);
			}

			delete[] attribute;

			if (__DUMP_CRYPTANNOT_MSGS__)
			{
				fprintf(stderr, "%s, %s: %s\n", gene_id, transcript_id, cur_line);
			}

			// Free memory.
			t_string::clean_tokens(attribute_toks);

			// Allocate the exon region.
			t_annot_region* exon_reg = get_empty_region();
			exon_reg->chrom = t_string::copy_me_str(chrom);

			// Make sure to extend the exon coordinates.
			exon_reg->start = translate_coord(atoi(start_str) - l_exon_ext, GFF_COORDS::start_base, CODEBASE_COORDS::start_base);
			exon_reg->end = translate_coord(atoi(end_str) + l_exon_ext, GFF_COORDS::end_base, CODEBASE_COORDS::end_base);
			exon_reg->strand = strand_str[0];
			exon_reg->score = cur_element_frame;
			exon_reg->name = t_string::copy_me_str(transcript_id);

			// Assign the exon to the gene/transcript.
			char* element_id = gene_id;

			int i_element = t_string::fast_search_string_per_prefix(element_id, gene_ids, 0, gene_ids->size());
			while (i_element > 0 &&
				t_string::sort_strings_per_prefix(element_id, gene_ids->at(i_element)))
			{
				i_element--;
			} // i_element loop.

			bool found_element = false;
			while (i_element < gene_ids->size() &&
				(t_string::sort_strings_per_prefix(gene_ids->at(i_element), element_id) || t_string::compare_strings(gene_ids->at(i_element), element_id)))
			{
				if (t_string::compare_strings(gene_ids->at(i_element), element_id))
				{
					found_element = true;
					break;
				}

				i_element++;
			} // i_element loop.

			if (found_element)
			{
				per_gene_exons->at(i_element)->push_back(exon_reg);
			}
			else
			{
				fprintf(stderr, "Could not find %s        \r", element_id);
				//exit(0);
			}
		} // exon check.

		delete[] cur_line;
	} // file reading loop.

	close_f(f_gff, GENCODE_gff_fp);

	// Save the regions.
	fprintf(stderr, "\nSaving the exons for each element.\n");
	FILE* f_op = open_f(op_fp, "w");
	FILE* f_redundant_exons = open_f("redundant_exons.bed", "w");
	for (int i_el = 0; i_el < gene_ids->size(); i_el++)
	{
		vector<t_annot_region*>* cur_gene_exons = per_gene_exons->at(i_el);
		fprintf(stderr, "%s: %d exons          \r", gene_ids->at(i_el), cur_gene_exons->size());

		for (int i_ex = 0; i_ex < cur_gene_exons->size(); i_ex++)
		{
			bool cur_reg_new = true;
			for (int j_ex = 0; j_ex < i_ex; j_ex++)
			{
				if (cur_gene_exons->at(i_ex)->start == cur_gene_exons->at(j_ex)->start &&
					cur_gene_exons->at(i_ex)->end == cur_gene_exons->at(j_ex)->end &&
					cur_gene_exons->at(i_ex)->score == cur_gene_exons->at(j_ex)->score)
				{
					cur_reg_new = false;
					break;
				}
			} // j_ex loop.

			if (cur_reg_new)
			{
				fprintf(f_op, "%s\t%d\t%d\t%s\t.\t%c\n",
					cur_gene_exons->at(i_ex)->chrom,
					translate_coord(cur_gene_exons->at(i_ex)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(cur_gene_exons->at(i_ex)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					cur_gene_exons->at(i_ex)->name,
					cur_gene_exons->at(i_ex)->strand);
			}
			else
			{
				fprintf(f_redundant_exons, "%s\t%d\t%d\t%s\t.\t%c\n",
					cur_gene_exons->at(i_ex)->chrom,
					translate_coord(cur_gene_exons->at(i_ex)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(cur_gene_exons->at(i_ex)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					cur_gene_exons->at(i_ex)->name,
					cur_gene_exons->at(i_ex)->strand);
			}
		} // i_ex loop.
	} // i_el loop.
	fclose(f_redundant_exons);
	fclose(f_op);
} // get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF

// **** THE MOST IMPORTANT IDEA BEHIND EOIs IS THAT THEY CAN OVERLAP. **** // 
// TODO: SANITY CHECKS FOR THE EOI DATA STRUCTURES.
// TODO: SIGNALIZING INDELs FROM VCFs?
// TODO: STREAMING ENCODING/ENCRYPTION OF VCF FILES AT EOIs!
// TODO: UNENCRYPTED ANNOTATION USING MULTIPLICATION IDEA: MULTIPLY, RE-GENERATE THE ANNOTATIONS.
void get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF(char* element_ids_list_fp, char* GENCODE_gff_fp, char* CDS_exon_selector, int eoi_element_type, int l_exon_ext, char* op_fp)
{
	// Read GFF file, load genes/transcripts and write exons for all the transcripts and genes.
	FILE* f_gff = open_f(GENCODE_gff_fp, "r");

	if (!t_string::compare_strings(CDS_exon_selector, "CDS") &&
		!t_string::compare_strings(CDS_exon_selector, "exon"))
	{
		fprintf(stderr, "CDS/exon selector (%s) should be one of \"CDS\" or \"exon\".\n", CDS_exon_selector);
		exit(0);
	}

	fprintf(stderr, "Generating exonic EOI regions using %s\n", GENCODE_gff_fp);

	vector<char*>* element_id_list = buffer_file(element_ids_list_fp);
	fprintf(stderr, "Loaded %d elements, setting up the exon region lists.\n", element_id_list->size());

	// Sort element id's.
	sort(element_id_list->begin(), element_id_list->end(), t_string::sort_strings_per_prefix);

	vector<vector<t_annot_region*>*>* per_element_exons = new vector<vector<t_annot_region*>*>();
	for (int i_el = 0; i_el < element_id_list->size(); i_el++)
	{
		per_element_exons->push_back(new vector<t_annot_region*>());
	} // i_el loop.

	  // GFF file reading.
	int n_entries = 0;
	while (1)
	{
		char* cur_line = getline(f_gff);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if ((n_entries > 0) &&
			(n_entries % 10000) == 0)
		{
			fprintf(stderr, "Processing %d. entry: %d elements.            \r", n_entries, (int)(element_id_list->size()));
		}

		n_entries++;

		// If this is an exonic entry, extract its element id.
		int i_cur_char = 0;
		char strand_char = '+';
		char region_name[1000];
		strcpy(region_name, ".");
		int region_score = 0;

		char chrom[100];
		t_string::get_next_token(cur_line, chrom, 100, "\t", i_cur_char);
		char source[100];
		t_string::get_next_token(cur_line, source, 100, "\t", i_cur_char);
		char element_type[100];
		t_string::get_next_token(cur_line, element_type, 100, "\t", i_cur_char);

		if (t_string::compare_strings(element_type, CDS_exon_selector))
		{
			char start_str[100];
			char end_str[100];
			char score_str[100];
			char strand_str[100];
			char frame_str[100];
			char* attribute = new char[strlen(cur_line) + 2];
			t_string::get_next_token(cur_line, start_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, end_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, score_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, strand_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, frame_str, 100, "\t", i_cur_char);
			t_string::get_next_token(cur_line, attribute, strlen(cur_line), "\t", i_cur_char);

			t_string_tokens* attribute_toks = t_string::tokenize_by_chars(attribute, "; =");
			char* gene_id = NULL;
			char* transcript_id = NULL;
			for (int i_t = 0; i_t + 1 < attribute_toks->size(); i_t++)
			{
				if (t_string::compare_strings(attribute_toks->at(i_t)->str(), "gene_id"))
				{
					gene_id = t_string::copy_me_str(attribute_toks->at(i_t + 1)->str());
				}

				if (t_string::compare_strings(attribute_toks->at(i_t)->str(), "transcript_id"))
				{
					transcript_id = t_string::copy_me_str(attribute_toks->at(i_t + 1)->str());
				}
			} // i_t loop.

			if (transcript_id == NULL || gene_id == NULL)
			{
				fprintf(stderr, "Could not parse gene/transcript id's: %s\n", cur_line);
				exit(0);
			}

			delete[] attribute;

			if (__DUMP_CRYPTANNOT_MSGS__)
			{
				fprintf(stderr, "%s, %s: %s\n", gene_id, transcript_id, cur_line);
			}

			// Free memory.
			t_string::clean_tokens(attribute_toks);

			// Allocate the exon region.
			t_annot_region* exon_reg = get_empty_region();
			exon_reg->chrom = t_string::copy_me_str(chrom);

			// Make sure to extend the exon coordinates.
			exon_reg->start = translate_coord(atoi(start_str) - l_exon_ext, GFF_COORDS::start_base, CODEBASE_COORDS::start_base);
			exon_reg->end = translate_coord(atoi(end_str) + l_exon_ext, GFF_COORDS::end_base, CODEBASE_COORDS::end_base);
			exon_reg->strand = strand_str[0];
			exon_reg->name = NULL;

			// Assign the exon to the gene/transcript.
			char* element_id = NULL;
			if (eoi_element_type == EOI_GENE_ELEMENTS)
			{
				element_id = gene_id;
			}
			else if (eoi_element_type == EOI_TRANSCRIPT_ELEMENTS)
			{
				element_id = transcript_id;
			}

			int i_element = t_string::fast_search_string_per_prefix(element_id, element_id_list, 0, element_id_list->size());
			while (i_element > 0 &&
				t_string::sort_strings_per_prefix(element_id, element_id_list->at(i_element)))
			{
				i_element--;
			} // i_element loop.

			bool found_element = false;
			while (i_element < element_id_list->size() &&
				(t_string::sort_strings_per_prefix(element_id_list->at(i_element), element_id) || t_string::compare_strings(element_id_list->at(i_element), element_id)))
			{
				if (t_string::compare_strings(element_id_list->at(i_element), element_id))
				{
					found_element = true;
					break;
				}

				i_element++;
			} // i_element loop.

			if (found_element)
			{
				per_element_exons->at(i_element)->push_back(exon_reg);
			}
			else
			{
				fprintf(stderr, "Could not find %s        \r", element_id);
				//exit(0);
			}
		} // exon check.

		delete[] cur_line;
	} // file reading loop.

	close_f(f_gff, GENCODE_gff_fp);

	// Save the regions.
	fprintf(stderr, "\nSaving the exons for each element.\n");
	FILE* f_op = open_f(op_fp, "w");
	for (int i_el = 0; i_el < element_id_list->size(); i_el++)
	{
		vector<t_annot_region*>* merged_exons = merge_annot_regions(per_element_exons->at(i_el), 1);
		for (int i_reg = 0; i_reg < merged_exons->size(); i_reg++)
		{
			fprintf(f_op, "%s\t%d\t%d\t%s\t.\t%c\n",
				merged_exons->at(i_reg)->chrom,
				translate_coord(merged_exons->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(merged_exons->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				element_id_list->at(i_el),
				merged_exons->at(i_reg)->strand);
		} // i_reg loop.

		delete_annot_regions(merged_exons);
	} // i_el loop.
	fclose(f_op);
} // get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF
