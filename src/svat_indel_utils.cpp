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

char* track_original_del_neigh_sequence_per_overlap_posns(vector<t_annot_region*>* sorted_cur_chr_EOI_regs, int cur_EOI_reg_i, 
	int& first_nuc_frame,
	char* ref_allele_seq,
	unsigned long long* annotated_deletion_signal,
	int overlap_start, int overlap_end, 
	int n_bits_per_nuc, int n_neigh_nucs,
	vector<t_VEP_term_ctx*>* vep_term_ctx)
{
	char* deleted_CDS_nucs = new char[overlap_end - overlap_start + 8];
	memset(deleted_CDS_nucs, 0, sizeof(char) * (overlap_end - overlap_start + 8));
	int del_cds_nuc_i = 0;

	// This must be initialized to -1 otherwise is not set correctly for the positive strands.
	first_nuc_frame = -1;

	unsigned long long int last_cds_nuc_annotation_signal = -1;

	for (int del_nuc_posn = overlap_start;
		del_nuc_posn <= overlap_end;
		del_nuc_posn++)
	{
		// Copy the annotation signal value.
		int track_posn = sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (del_nuc_posn - sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

		if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
		{
			fprintf(stderr, "Setting Deletion Overlap @ %d (%d)\n",
				del_nuc_posn, track_posn);
			getc(stdin);
		}

		unsigned long long cur_del_annotation_signal_val = annotated_deletion_signal[track_posn];

		// Unpack and translate the indel.
		char s_neigh_seq_buff[100];
		unsigned int s_coding_frame = 0;
		unsigned long long int s_impact_val = 0;
		unpack_impact_signal_values(cur_del_annotation_signal_val,
									s_coding_frame, s_impact_val, s_neigh_seq_buff,
									n_bits_per_nuc, n_neigh_nucs);

		// Is this annotation signal pointing to a CDS context?
		if(check_element_context_per_impact_bitmap(vep_term_ctx, CDS_ctx, s_impact_val))
		{
			last_cds_nuc_annotation_signal = cur_del_annotation_signal_val;

			if (first_nuc_frame == -1 &&
				sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand == '+')
			{
				first_nuc_frame = s_coding_frame;
			}

			// If the target region is on the negative strand, set the nucleotide frame at every position where there is a coding nucleotide.
			if (sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand == '-')
			{
				first_nuc_frame = s_coding_frame;
			}

			// Append the buffering nucleotides at the beginning.
			if (del_cds_nuc_i == 0)
			{				
				for (int nuc_i = 0; 
					nuc_i < n_neigh_nucs / 2; 
					nuc_i++)
				{
					deleted_CDS_nucs[del_cds_nuc_i] = s_neigh_seq_buff[nuc_i];
					del_cds_nuc_i++;
				} // nuc_i loop
			}

			int del_nuc_i = del_nuc_posn - overlap_start + 1;

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "Adding actual del nuc: %c (@%d of %s)\n", ref_allele_seq[del_nuc_i], del_nuc_i, ref_allele_seq);
			}

			deleted_CDS_nucs[del_cds_nuc_i] = ref_allele_seq[del_nuc_i];
			del_cds_nuc_i++;
		} // cds context check.
	} // overlap position.

	// Add the buffer if there were coding nucleotides.
	if (del_cds_nuc_i > 0)
	{
		// Pad this sequence with the right side.
		char s_neigh_seq_buff[100];
		unsigned int s_coding_frame = 0;
		unsigned long long int s_impact_val = 0;
		unpack_impact_signal_values(last_cds_nuc_annotation_signal,
			s_coding_frame, s_impact_val, s_neigh_seq_buff,
			n_bits_per_nuc, n_neigh_nucs);

		fprintf(stderr, "Last CDS nuc's neighborhood: %s\n", s_neigh_seq_buff);

		for (int nuc_i = n_neigh_nucs / 2; nuc_i < n_neigh_nucs; nuc_i++)
		{
			deleted_CDS_nucs[del_cds_nuc_i] = s_neigh_seq_buff[nuc_i];
			del_cds_nuc_i++;
		} // nuc_i loop
	}

	return(deleted_CDS_nucs);
}

vector<t_annot_region*>* load_transcript_CDS_sequences(char* per_transcript_CDS_regs, char* genome_seq_dir)
{
	vector<t_annot_region*>* transcript_CDS_regs = load_BED(per_transcript_CDS_regs);

	t_restr_annot_region_list* restr_CDS_regs = restructure_annot_regions(transcript_CDS_regs);
	for (int i_chr = 0; i_chr < restr_CDS_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Loading CDS seqeunces for transcripts on %s\n", restr_CDS_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_cds_regs = restr_CDS_regs->regions_per_chrom[i_chr];
		sort(cur_chr_cds_regs->begin(), cur_chr_cds_regs->end(), sort_regions_per_name);

		char cur_chrom_seq_fp[1000];
		sprintf(cur_chrom_seq_fp, "%s/%s.bin.gz", genome_seq_dir, restr_CDS_regs->chr_ids->at(i_chr));
		int l_seq = 0;
		char* cur_chrom_seq = load_binary_sequence_file(cur_chrom_seq_fp, l_seq);
		fprintf(stderr, "Loaded %d nucleotides.\n", l_seq);

		// Go over all the regions.
		int i_reg = 0;
		while (i_reg < cur_chr_cds_regs->size())
		{
			int start_i = i_reg;
			int end_i = i_reg;
			vector<t_annot_region*>* cur_gene_CDS_regs = new vector<t_annot_region*>();
			int j_reg = i_reg;
			for (; j_reg < cur_chr_cds_regs->size(); j_reg++)
			{
				if (!t_string::compare_strings(cur_chr_cds_regs->at(j_reg)->name, cur_chr_cds_regs->at(i_reg)->name))
				{
					break;
				}

				cur_gene_CDS_regs->push_back(cur_chr_cds_regs->at(j_reg));
			} // j_reg loop.

			// Sort the CDS regions.
			sort(cur_gene_CDS_regs->begin(), cur_gene_CDS_regs->end(), sort_regions);

			fprintf(stderr, "%s: %d CDSs\n", cur_chr_cds_regs->at(i_reg)->name, cur_gene_CDS_regs->size());

			t_string* cur_gene_CDS = new t_string();
			for (int i = 0; i < cur_gene_CDS_regs->size();i++)
			{
				for (int j = cur_gene_CDS_regs->at(i)->start; j <= cur_gene_CDS_regs->at(i)->end; j++)
				{
					cur_gene_CDS->concat_char(cur_chrom_seq[j]);
				} // j loop. 
			} // i loop.

			t_annot_region* cur_transcript_reg = get_empty_region();
			cur_transcript_reg->chrom = t_string::copy_me_str(cur_gene_CDS_regs->at(0)->chrom);
			cur_transcript_reg->start = cur_gene_CDS_regs->at(0)->start;
			cur_transcript_reg->end = cur_gene_CDS_regs->back()->end;
			cur_transcript_reg->name = t_string::copy_me_str(cur_gene_CDS_regs->at(0)->name);
			char* cur_transcript_CDS_seq = t_string::copy_me_str(cur_gene_CDS->str());

			if (cur_transcript_reg->start == '-')
			{
				char* rev_cur_transcript_CDS_seq = get_reverse_complement(cur_transcript_CDS_seq, 0, strlen(cur_transcript_CDS_seq));
				delete[] cur_transcript_CDS_seq;
				cur_transcript_CDS_seq = rev_cur_transcript_CDS_seq;
			}
			cur_transcript_reg->data = cur_transcript_CDS_seq;

			transcript_CDS_regs->push_back(cur_transcript_reg);

			delete cur_gene_CDS_regs;

			i_reg = j_reg;
		} // i_reg loop.

		delete[] cur_chrom_seq;
	} // i_chr loop.

	return(transcript_CDS_regs);
}

// This must be updated to handle 
void count_deleted_elements_left2right_impact_states(vector<unsigned long long>* all_impact_index_bitmaps,
													vector<t_VEP_term_ctx*>* vep_term_ctx,
													int& n_deleted_coding_nucs, int& n_impacted_CDSs,
													int& n_del_splice_region_nucs, int& n_impacted_splice_regs,
													int& n_del_splice_acc_nucs, int& n_impacted_splice_acc_sites,
													int& n_del_splice_don_nucs, int& n_impacted_splice_don_sites,
													int& n_del_start_codon_nucs, int& n_impacted_start_codon_sites,
													int& n_del_stop_codon_nucs, int& n_impacted_stop_codon_sites,
													int& n_del_intron_nucs, int& n_impacted_introns,
													int& n_del_fp_utr_nucs, int& n_impacted_fp_utrs,
													int& n_del_tp_utr_nucs, int& n_impacted_tp_utrs)
{
	for (int del_nuc_i = 0; del_nuc_i < all_impact_index_bitmaps->size(); del_nuc_i++)
	{
		if (check_element_context_per_impact_bitmap(vep_term_ctx, CDS_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_deleted_coding_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, CDS_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_CDSs++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_acceptor_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_splice_acc_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, splice_acceptor_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_splice_acc_sites++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_donor_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_splice_don_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, splice_donor_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_splice_don_sites++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_region_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_splice_region_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, splice_region_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_splice_regs++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, intron_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_intron_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, intron_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_introns++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, fp_UTR_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_fp_utr_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, fp_UTR_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_fp_utrs++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, tp_UTR_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_tp_utr_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, tp_UTR_ctx, all_impact_index_bitmaps->at(del_nuc_i-1)))
			{
				n_impacted_tp_utrs++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, start_codon_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_start_codon_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, start_codon_ctx, all_impact_index_bitmaps->at(del_nuc_i - 1)))
			{
				n_impacted_start_codon_sites++;
			}
		}

		if (check_element_context_per_impact_bitmap(vep_term_ctx, stop_codon_ctx, all_impact_index_bitmaps->at(del_nuc_i)))
		{
			n_del_stop_codon_nucs++;

			if (del_nuc_i == 0 ||
				!check_element_context_per_impact_bitmap(vep_term_ctx, stop_codon_ctx, all_impact_index_bitmaps->at(del_nuc_i - 1)))
			{
				n_impacted_stop_codon_sites++;
			}
		}
	} // del_nuc_i loop.
}

void generate_EOI_randomized_Indel_VCF(char* EOI_bed_fp, int l_max_indel, double per_posn_var_prob, double del_prob, char* bin_seq_dir, char* op_vcf_fp)
{
	fprintf(stderr, "Simulating indels of maximum length %d.\n", l_max_indel);

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

		char* reference_allele = new char[l_max_indel + 2];
		char* alternate_allele = new char[l_max_indel + 2];

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_eoi_regs->regions_per_chrom[i_chr];
		fprintf(stderr, "Writing the Indels for %d regions.\n", (int)(cur_chr_EOI_regs->size()));
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

				int l_var = (int)(rng->random_double_ran3() * l_max_indel);
				if (l_var == 0)
				{
					l_var = 1;
				}

				// Select deletion/insertion.
				if (rng->random_double_ran3() > del_prob ||
					del_prob < 0.001)
				{
					// Insertions
					// Setup the reference allele.
					//memset(reference_allele, 0, sizeof(char) * (l_max_indel + 2));
					//memset(alternate_allele, 0, sizeof(char) * (l_max_indel + 2));

					// Set the inserted nucleotides to N's.
					for (int nuc_i = i - 1; nuc_i < i + l_var; nuc_i++)
					{
						char rand_nuc = num_2_nuc((int)(floor(rng->random_double_ran3() * 4)));
						alternate_allele[nuc_i - i + 1] = rand_nuc;
					} // nuc_i loop.
					alternate_allele[l_var + 1] = 0;
					alternate_allele[0] = cur_chr_seq[i - 1];

					reference_allele[0] = cur_chr_seq[i - 1];
					reference_allele[1] = 0;

					// #CHROM POS     ID        REF ALT    QUAL FILTER
					fprintf(f_op_vcf, "%s\t%d\t%s_%d_INS_%s\t%s\t%s\t60\tPASS\t.\tGT\t0|1\n",
						restr_eoi_regs->chr_ids->at(i_chr),
						translate_coord(i - 1, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
						restr_eoi_regs->chr_ids->at(i_chr), i - 1, cur_chr_EOI_regs->at(i_reg)->name,
						reference_allele,
						alternate_allele);
				}
				else
				{
					// Setup the reference allele.
					for (int nuc_i = i - 1; nuc_i < i + l_var; nuc_i++)
					{
						reference_allele[nuc_i - i + 1] = cur_chr_seq[nuc_i];
					} // nuc_i loop.
					reference_allele[l_var + 1] = 0;
					reference_allele[0] = cur_chr_seq[i - 1];

					alternate_allele[0] = cur_chr_seq[i - 1];
					alternate_allele[1] = 0;

					// #CHROM POS     ID        REF ALT    QUAL FILTER
					fprintf(f_op_vcf, "%s\t%d\t%s_%d_DEL_%s\t%s\t%s\t60\tPASS\t.\tGT\t0|1\n",
						restr_eoi_regs->chr_ids->at(i_chr),
						translate_coord(i - 1, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
						restr_eoi_regs->chr_ids->at(i_chr), i - 1, cur_chr_EOI_regs->at(i_reg)->name,
						reference_allele,
						alternate_allele);
				} // deletion check.
			} // i loop.
		} // i_reg loop.

		delete[] cur_chr_seq;
	} // i_chr loop.
	close_f(f_op_vcf, op_vcf_fp);
}

void generate_EOI_enumerating_Indel_VCF(char* EOI_bed_fp, int l_max_indel, char* bin_seq_dir, char* op_vcf_fp)
{
	fprintf(stderr, "Enumerating indels of maximum length %d.\n", l_max_indel);

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

		char* reference_allele = new char[l_max_indel + 2];
		char* alternate_allele = new char[l_max_indel + 2];

		vector<t_annot_region*>* cur_chr_EOI_regs = restr_eoi_regs->regions_per_chrom[i_chr];
		fprintf(stderr, "Writing the Indels for %d regions.\n", (int)(cur_chr_EOI_regs->size()));
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			if (i_reg % 500 == 0)
			{
				fprintf(stderr, "Processing %d. region.            \r", i_reg);
			}

			for (int i = cur_chr_EOI_regs->at(i_reg)->start; i <= cur_chr_EOI_regs->at(i_reg)->end; i++)
			{
				// Insertions
				for (int l_var = 1; l_var <= l_max_indel; l_var++)
				{
					// Setup the reference allele.
					//memset(reference_allele, 0, sizeof(char) * (l_max_indel + 2));
					//memset(alternate_allele, 0, sizeof(char) * (l_max_indel + 2));

					// Set the inserted nucleotides to N's.
					for (int nuc_i = i - 1; nuc_i < i + l_var; nuc_i++)
					{
						alternate_allele[nuc_i - i + 1] = 'N';
					} // nuc_i loop.
					alternate_allele[l_var + 1] = 0;
					alternate_allele[0] = cur_chr_seq[i - 1];

					reference_allele[0] = cur_chr_seq[i - 1];
					reference_allele[1] = 0;

					// #CHROM POS     ID        REF ALT    QUAL FILTER
					fprintf(f_op_vcf, "%s\t%d\t%s_%d_INS_%s\t%s\t%s\t60\tPASS\n",
						restr_eoi_regs->chr_ids->at(i_chr),
						translate_coord(i - 1, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
						restr_eoi_regs->chr_ids->at(i_chr), i - 1, cur_chr_EOI_regs->at(i_reg)->name,
						reference_allele,
						alternate_allele);
				} // nuc_i loop.

				  // Deletions:
				for (int l_var = 1; l_var <= l_max_indel; l_var++)
				{
					// Setup the reference allele.
					/*memset(reference_allele, 0, sizeof(char) * (l_max_indel + 2));
					memset(alternate_allele, 0, sizeof(char) * (l_max_indel + 2));*/
					for (int nuc_i = i - 1; nuc_i < i + l_var; nuc_i++)
					{
						reference_allele[nuc_i - i + 1] = cur_chr_seq[nuc_i];
					} // nuc_i loop.
					reference_allele[l_var + 1] = 0;
					reference_allele[0] = cur_chr_seq[i - 1];

					alternate_allele[0] = cur_chr_seq[i - 1];
					alternate_allele[1] = 0;

					// #CHROM POS     ID        REF ALT    QUAL FILTER
					fprintf(f_op_vcf, "%s\t%d\t%s_%d_DEL_%s\t%s\t%s\t60\tPASS\n",
						restr_eoi_regs->chr_ids->at(i_chr),
						translate_coord(i - 1, CODEBASE_COORDS::start_base, VCF_COORDS::start_base),
						restr_eoi_regs->chr_ids->at(i_chr), i - 1, cur_chr_EOI_regs->at(i_reg)->name,
						reference_allele,
						alternate_allele);
				} // nuc_i loop.
			} // i loop.
		} // i_reg loop.

		delete[] cur_chr_seq;
	} // i_chr loop.
	close_f(f_op_vcf, op_vcf_fp);
}

void signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization_with_CDS_vicinity(char* EOI_regs_BED_fp,
	char* genome_dir,
	char* CDS_regs_fp,
	int EOI_element_type,
	char* vep_op_dir,
	char* sorted_impact_value_strings_list_fp,
	char* op_dir)
{
	fprintf(stderr, "***THIS CODE IS NOT YET COMPLETE, RUN AT YOUR OWN RISK. MAKE SURE THE DATA THAT IS REPLACED IS CORRECTED.***");
	getc(stdin);

	// Load the EOIs.
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	for (int i_reg = 0; i_reg < EOI_regs->size(); i_reg++)
	{
		EOI_regs->at(i_reg)->data = NULL;
	} // i_Reg loop.

	double total_EOI_covg = coverage(EOI_regs);

	vector<t_annot_region*>* transcript_regs = load_transcript_CDS_sequences(CDS_regs_fp, genome_dir);

	// Assign the CDS sequences to the transcripts.
	fprintf(stderr, "Assigning the CDS sequences to the EOI regions.\n");
	vector<t_annot_region*>* intersects = intersect_regions_per_names(EOI_regs, transcript_regs, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* eoi_reg = int_info->src_reg;
		t_annot_region* trans_reg = int_info->dest_reg;

		// Copy the sequence.
		eoi_reg->data = trans_reg->data;
	} // i_int loop.

	// Load the impact values from VEP's output.
	vector<char*>* sorted_impact_vals = new vector<char*>();
	vector<t_VEP_term_ctx*>* vep_term_ctx = load_VEP_annotation_term_context(sorted_impact_value_strings_list_fp, sorted_impact_vals);
	fprintf(stderr, "Loaded %d impact strings.\n", sorted_impact_vals->size());

	// The EOI's determine the coordinate system that we will use. The impact value should include this.
	t_restr_annot_region_list* restr_EOI_regs = restructure_annot_regions(EOI_regs);

	// Set the chromosome id's.
	vector<char*>* chr_ids = restr_EOI_regs->chr_ids;

	fprintf(stderr, "Allocating impact signals per chromosome.\n");
	FILE* f_unmatched = open_f("unmatched_vep_impacts.op.gz", "w");
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

		fprintf(stderr, "Allocating the deletion impact signal.\n");
		unsigned int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		unsigned long long* deletion_VEP_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		memset(deletion_VEP_signal, 0, sizeof(unsigned long long) * cur_chr_EOI_covg);

		// Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "%s:%d-%d @ %d\n",
					cur_chr_EOI_regs->at(i_reg)->name,
					cur_chr_EOI_regs->at(i_reg)->start, cur_chr_EOI_regs->at(i_reg)->end,
					cumul_signal_index);
			}
		} // i_reg loop.

		  // Sanity check: Last should be last.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

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
			t_string_tokens* location_str_toks = cur_line_toks->at(LOCATION_col_i)->tokenize_by_chars(":");
			char* var_chr = t_string::copy_me_str(location_str_toks->at(0)->str());
			int var_posn = atoi(location_str_toks->at(1)->str());
			char* var_name = cur_line_toks->at(VAR_ID_col_i)->str();
			t_string::clean_tokens(location_str_toks);

			// Set the coding frame: This is the coding frame of the deleted nucleotide; we convert it into 0-based coordinates and move on.
			int coding_frame = 0;
			int cds_posn = 0;
			char* cds_posn_str = cur_line_toks->at(CDS_POSN_col_i)->str();
			if (t_string::is_number(cds_posn_str) &&
				!t_string::compare_strings(cds_posn_str, "-"))
			{
				cds_posn = atoi(cds_posn_str);
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

			// Go over the consequence tokens and choose the highest impact.
			t_string_tokens* conseq_toks = cur_line_toks->at(CONSEQUENCE_col_i)->tokenize_by_chars(", ");

			//int cur_impact_val = -1;
			unsigned long long cur_impact_val = 0;
			for (int i_tok = 0; i_tok < conseq_toks->size(); i_tok++)
			{
				int cur_tok_impact_val = t_string::get_i_str(sorted_impact_vals, conseq_toks->at(i_tok)->str());

				if (cur_tok_impact_val == sorted_impact_vals->size())
				{
					fprintf(stderr, "Could not parse impact string: %s\n", cur_line_toks->at(CONSEQUENCE_col_i)->str());
					exit(0);
				}

				// Add this impact value to the list of impact values.				
				cur_impact_val = cur_impact_val | (LLONE << cur_tok_impact_val);

				if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
				{
					fprintf(stderr, "Adding %d. impact: %s (%d) :: Aggregate Impact: %llX\n", i_tok, conseq_toks->at(i_tok)->str(), cur_tok_impact_val, cur_impact_val);
				}
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

			// Copy the genomic neighborhood: Surrounding tri-nucleotides: We store the 6 nucleotides in the neighborhood, 32 bit value is more than enough.
			unsigned int geno_neigh_seq_signal = 0; // 32-bit neighbor value.
			int n_bits_per_nuc = 2;
			int n_neigh_nucs = 6;
			int cur_base = 0;
			for (int i = var_posn - (n_neigh_nucs/2); i < var_posn; i++)
			{
				unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
				geno_neigh_seq_signal = (geno_neigh_seq_signal | (cur_nuc_val << cur_base));
				cur_base += n_bits_per_nuc;
			} // i_l loop.

			for (int i = var_posn + 1; i <= var_posn + (n_neigh_nucs/2); i++)
			{
				unsigned int cur_nuc_val = nuc_2_num(cur_chrom_seq[i]);
				geno_neigh_seq_signal = (geno_neigh_seq_signal | (cur_nuc_val << cur_base));
				cur_base += n_bits_per_nuc;
			} // i_l loop.

			  // For every overlap, update the annotation signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				bool overlap = (cur_chr_EOI_regs->at(cur_EOI_reg_i)->start <= var_posn &&
					cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn);

				if (overlap)
				{
					// Make sure the impact assigned element is the same element that we are setting the signal at.
					if (t_string::compare_strings(vep_element_id, cur_chr_EOI_regs->at(cur_EOI_reg_i)->name))
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

						// Copy the CDS neighborhood.
						char* cur_EOI_cds_seq = (char*)(cur_chr_EOI_regs->at(cur_EOI_reg_i)->data);
						unsigned int cds_neigh_seq_signal = 0;
						if (cur_EOI_cds_seq != NULL)
						{
							fprintf(stderr, "%s: Left CDS neighborhood: ", var_name);
							cur_base = 0;
							int cds_posn_0_based = cds_posn - 1;
							int i = cds_posn_0_based;
							for (int i = MAX(0, cds_posn_0_based - (n_neigh_nucs / 2)); i < cds_posn_0_based; i++)
							{
								unsigned int cur_nuc_val = nuc_2_num(cur_EOI_cds_seq[i]);
								cds_neigh_seq_signal = (cds_neigh_seq_signal | (cur_nuc_val << cur_base));
								cur_base += n_bits_per_nuc;

								fprintf(stderr, "%c", cur_EOI_cds_seq[i]);
							} // i loop.

							fprintf(stderr, "\n%s: Right CDS neighborhood: ", var_name);
							for (int i = cds_posn_0_based + 1; i <= cds_posn_0_based + (n_neigh_nucs/2); i++)
							{
								if (cur_EOI_cds_seq[i] == 0)
								{
									break;
								}
								else
								{
									unsigned int cur_nuc_val = nuc_2_num(cur_EOI_cds_seq[i]);
									cds_neigh_seq_signal = (cds_neigh_seq_signal | (cur_nuc_val << cur_base));
									cur_base += n_bits_per_nuc;
								}

								fprintf(stderr, "%c", cur_EOI_cds_seq[i]);
							} // i loop.

							fprintf(stderr, "\n");
						} // CDS neighborhood generation.

						// Extract the previously set impact value.
						char prev_geno_neigh_seq_buff[100];
						char prev_cds_neigh_seq_buff[100];
						unsigned int prev_coding_frame = 0;
						unsigned long long prev_impact_val = 0;
						unpack_impact_signal_values_w_CDS_neighborhood(deletion_VEP_signal[track_posn], 
							prev_coding_frame, prev_impact_val, 
							prev_geno_neigh_seq_buff, prev_cds_neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						// Add the previous impact value to this.
						cur_impact_val = cur_impact_val | prev_impact_val;

						// If the new impact is higher, assign it.
						// Assign score to the EOI impact score track.	
						unsigned long long cur_signal = pack_set_impact_signal_values_w_CDS_neighborhood(coding_frame,
							cur_impact_val,
							geno_neigh_seq_signal, cds_neigh_seq_signal, n_bits_per_nuc, n_neigh_nucs);

						deletion_VEP_signal[track_posn] = cur_signal;

						// Sanity check to make sure packing worked.
						char geno_neigh_seq_buff[100];
						char cds_neigh_seq_buff[100];
						unsigned int s_coding_frame = 0;
						unsigned long long s_impact_val = 0;
						unpack_impact_signal_values_w_CDS_neighborhood(cur_signal, s_coding_frame, s_impact_val, geno_neigh_seq_buff, cds_neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						if (s_impact_val != cur_impact_val ||
							s_coding_frame != coding_frame)
						{
							fprintf(stderr, "Sanity check failed: Coding frame or the impact value is not coded correctly:\n\
%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s // %s (%u, %u)\n",
cur_line,
s_impact_val,
s_coding_frame,
geno_neigh_seq_buff, cds_neigh_seq_buff, geno_neigh_seq_signal, cds_neigh_seq_signal);
							exit(0);
						}

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s (%u)\n",
								cur_line,
								s_impact_val,
								s_coding_frame,
								geno_neigh_seq_buff, cds_neigh_seq_signal);
						}

						found_overlap = true;

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
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
		fprintf(stderr, "Writing insertion impact signal.\n");
		fwrite(deletion_VEP_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_vep_impacts);
		delete[] deletion_VEP_signal;
		close_f(f_vep_impacts, cur_chr_impacts_op_fp);

		// Write the BED file with sorted regions that match the signal coordinates in signal file.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, cur_chr_EOI_regs);

		delete[] cur_chrom_seq;
	} // i_chr loop.
	fclose(f_unmatched);

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

void signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
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

	char unmatched_impacts_fp[] = "unmatched_vep_impacts.op.gz";
	FILE* f_unmatched = open_f(unmatched_impacts_fp, "w");
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

		fprintf(stderr, "Allocating the deletion impact signal.\n");
		unsigned int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		unsigned long long* deletion_VEP_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		memset(deletion_VEP_signal, 0, sizeof(unsigned long long) * cur_chr_EOI_covg);		

		// Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "%s:%d-%d @ %d\n", 
						cur_chr_EOI_regs->at(i_reg)->name, 
						cur_chr_EOI_regs->at(i_reg)->start, cur_chr_EOI_regs->at(i_reg)->end,
						cumul_signal_index);
			}
		} // i_reg loop.

		// Sanity check: Last should be last.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// Initialize the deletion signal with the neighboring signal.
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

				deletion_VEP_signal[track_posn] = cur_signal;
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
			t_string_tokens* location_str_toks = cur_line_toks->at(LOCATION_col_i)->tokenize_by_chars(":");
			char* var_chr = t_string::copy_me_str(location_str_toks->at(0)->str());
			int var_posn = atoi(location_str_toks->at(1)->str());			
			t_string::clean_tokens(location_str_toks);

			// Set the coding frame: This is the coding frame of the deleted nucleotide; we convert it into 0-based coordinates and move on.
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

			// Go over the consequence tokens and choose the highest impact.
			t_string_tokens* conseq_toks = cur_line_toks->at(CONSEQUENCE_col_i)->tokenize_by_chars(", ");

			//int cur_impact_val = -1;
			unsigned long long cur_impact_val = 0;
			for (int i_tok = 0; i_tok < conseq_toks->size(); i_tok++)
			{
				int cur_tok_impact_val = t_string::get_i_str(sorted_impact_vals, conseq_toks->at(i_tok)->str());

				if (cur_tok_impact_val == sorted_impact_vals->size())
				{
					fprintf(stderr, "Could not parse impact string: %s\n", cur_line_toks->at(CONSEQUENCE_col_i)->str());
					exit(0);
				}

				// Add this impact value to the list of impact values.				
				cur_impact_val = cur_impact_val | (LLONE << cur_tok_impact_val);

				if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
				{
					fprintf(stderr, "Adding %d. impact: %s (%d) :: Aggregate Impact: %llX\n", i_tok, conseq_toks->at(i_tok)->str(), cur_tok_impact_val, cur_impact_val);
				}
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

			// For every overlap, update the annotation signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				bool overlap = (cur_chr_EOI_regs->at(cur_EOI_reg_i)->start <= var_posn &&
					cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn);

				if (overlap)
				{
					// Deletion is engulfed, we should be able to set a vectorized "track" position.
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
						// Extract the previously set impact value.
						char prev_neigh_seq_buff[100];
						unsigned int prev_coding_frame = 0;
						unsigned long long prev_impact_val = 0;
						unpack_impact_signal_values(deletion_VEP_signal[track_posn], prev_coding_frame, prev_impact_val, prev_neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						// Add the previous impact value to this.
						cur_impact_val = cur_impact_val | prev_impact_val;

						// If the new impact is higher, assign it.
						// Assign score to the EOI impact score track.	
						unsigned long long cur_signal = pack_set_impact_signal_values(coding_frame,
																				cur_impact_val,
																				neigh_seq_signal, n_bits_per_nuc, n_neigh_nucs);

						deletion_VEP_signal[track_posn] = cur_signal;

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

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s (%u)\n", 
								cur_line, 
								s_impact_val, 
								s_coding_frame,
								neigh_seq_buff, neigh_seq_signal);
						}

						found_overlap = true;

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
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
		fprintf(stderr, "Writing deletion impact signal.\n");
		fwrite(deletion_VEP_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_vep_impacts);
		delete[] deletion_VEP_signal;
		close_f(f_vep_impacts, cur_chr_impacts_op_fp);

		// Write the BED file with sorted regions that match the signal coordinates in signal file.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, cur_chr_EOI_regs);

		delete[] cur_chrom_seq;
	} // i_chr loop.
	close_f(f_unmatched, unmatched_impacts_fp);

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

void signalize_VEP_annotated_Inserts_per_EOI_regs_element_summarization(char* EOI_regs_BED_fp,
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
	char unmatched_impacts_fp[] = "unmatched_vep_impacts.op.gz";
	FILE* f_unmatched = open_f(unmatched_impacts_fp, "w");
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

		fprintf(stderr, "Allocating the insertion impact signal.\n");
		unsigned int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		unsigned long long* insertion_VEP_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		memset(insertion_VEP_signal, 0, sizeof(unsigned long long) * cur_chr_EOI_covg);

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

		// Initialize the deletion signal with the neighboring signal.
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			for (int var_posn = cur_chr_EOI_regs->at(i_reg)->start;
				var_posn <= cur_chr_EOI_regs->at(i_reg)->end;
				var_posn++)
			{
				int track_posn = (cur_chr_EOI_regs->at(i_reg)->score + var_posn - cur_chr_EOI_regs->at(i_reg)->start);

				// Copy the surrounding tri-nucleotides: We store the 6 nucleotides in the neighborhood, 32 bit value is more than enough.
				// Note that the left sequence (left 3 nucleotides) contain the current nucleotide.
				// Note that this should be coding sequences.
				unsigned int neigh_seq_signal = 0; // 32-bit neighbor value.
				int cur_base = 0;
				for (int i = var_posn - 2; i <= var_posn; i++)
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

				insertion_VEP_signal[track_posn] = cur_signal;
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

		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		unsigned long long LLONE = 1;

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

			char* var_name = cur_line_toks->at(VAR_ID_col_i)->str();

			// Parse the location.
			t_string_tokens* location_str_toks = cur_line_toks->at(LOCATION_col_i)->tokenize_by_chars(":");
			char* var_chr = t_string::copy_me_str(location_str_toks->at(0)->str());
			int var_posn = atoi(location_str_toks->at(1)->str());
			t_string::clean_tokens(location_str_toks);

			// Set the coding frame.
			int coding_frame = 0;
			char* cds_posn_str = cur_line_toks->at(CDS_POSN_col_i)->str();			

			t_string_tokens* CDS_posn_tokens = t_string::tokenize_by_chars(cds_posn_str, "-");
			if(CDS_posn_tokens->size() == 2)
			{
				int cds_posn = atoi(CDS_posn_tokens->at(0)->str());
				coding_frame = ((cds_posn - 1) % 3);

				t_string::clean_tokens(CDS_posn_tokens);
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

				// Update the impact value, if necessary.
				cur_impact_val = cur_impact_val | (LLONE << cur_tok_impact_val);
				
				if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
				{
					fprintf(stderr, "Adding %d. impact: %s (%d) :: Aggregate Impact: %llX\n", i_tok, conseq_toks->at(i_tok)->str(), cur_tok_impact_val, cur_impact_val);
				}
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
			// Note that the left sequence (left 3 nucleotides) contain the current nucleotide.
			// Note that this should be coding sequences.
			unsigned int neigh_seq_signal = 0; // 32-bit neighbor value.
			int cur_base = 0;
			for (int i = var_posn - 2; i <= var_posn; i++)
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

			// For every overlap, update the annotation signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				int overlap_start = MAX(cur_chr_EOI_regs->at(cur_EOI_reg_i)->start, var_posn);
				int overlap_end = MIN(cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn + 1);

				if (overlap_start <= overlap_end &&
					overlap_start <= var_posn &&
					overlap_end >= var_posn + 1)
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
						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "%s:%d / %s @ %d : %llX\n", var_name, var_posn, cur_chr_EOI_regs->at(cur_EOI_reg_i)->name, track_posn, cur_impact_val);
						}

						char prev_neigh_seq_buff[100];
						unsigned int prev_coding_frame = 0;
						unsigned long long prev_impact_val = 0;
						unpack_impact_signal_values(insertion_VEP_signal[track_posn], prev_coding_frame, prev_impact_val, prev_neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						// Update the current impact value with the previous one.
						cur_impact_val  = cur_impact_val | prev_impact_val;

						// Assign score to the EOI impact score track.	
						unsigned long long cur_signal = pack_set_impact_signal_values(coding_frame,
																					cur_impact_val,
																					neigh_seq_signal, n_bits_per_nuc, n_neigh_nucs);

						insertion_VEP_signal[track_posn] = cur_signal;

						// Sanity check to make sure packing worked.
						char neigh_seq_buff[100];
						unsigned int s_coding_frame = 0;
						unsigned long long s_impact_val = 0;
						unpack_impact_signal_values(cur_signal, s_coding_frame, s_impact_val, neigh_seq_buff, n_bits_per_nuc, n_neigh_nucs);

						if (s_impact_val != cur_impact_val ||
							s_coding_frame != coding_frame)
						{
							fprintf(stderr, "Sanity check failed: Coding frame or the impact value is not coded correctly:\n\
%s\nimpact_val: %llX (%llX);; coding_frame: %u\nNeighbor Seq: %s (%u)\n",
								cur_line,
								s_impact_val, cur_impact_val,
								s_coding_frame,
								neigh_seq_buff, neigh_seq_signal);
							exit(0);
						}

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "%s\nimpact_val: %llX;; coding_frame: %u\nNeighbor Seq: %s (%u)\n",
								cur_line,
								s_impact_val,
								s_coding_frame,
								neigh_seq_buff, neigh_seq_signal);
						}

						found_overlap = true;

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
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
		fprintf(stderr, "Writing insertion impact signal.\n");
		fwrite(insertion_VEP_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_vep_impacts);
		delete[] insertion_VEP_signal;
		close_f(f_vep_impacts, cur_chr_impacts_op_fp);

		// Write the BED file with sorted regions that match the signal coordinates in signal file.
		char cur_sorted_BED_fp[1000];
		sprintf(cur_sorted_BED_fp, "%s/coord_matching_sorted_regions_%s.bed", op_dir, chr_ids->at(i_chr));
		dump_BED(cur_sorted_BED_fp, cur_chr_EOI_regs);

		delete[] cur_chrom_seq;
	} // i_chr loop.
	close_f(f_unmatched, unmatched_impacts_fp);

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

// Write the ternary signal on the EOI=gene_exons.
void signalize_VCF_Insertion_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
	char* per_chrom_VCF_dir,
	char* op_dir)
{
	fprintf(stderr, "Signalizing Inserts under %s over the target regions in %s and saving to %s\n", per_chrom_VCF_dir, EOI_regs_BED_fp, op_dir);

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
		
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));

		// Allocat the per indel per l_indel selector.
		fprintf(stderr, "Allocating and initializing the %d long insertion VCF signal.\n", cur_chr_EOI_covg);
		unsigned int* insertion_genotype_signal = new unsigned int[cur_chr_EOI_covg + 2];
		memset(insertion_genotype_signal, 0, sizeof(unsigned int) * cur_chr_EOI_covg);

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
			fprintf(stderr, "Reading inserts in %s\n", cur_chrom_vcf_fp);
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
			if (cur_line_toks->size() < (GENOTYPE_col_i + 1))
			{
				fprintf(stderr, "Expecting at least %d columns in VCF file.\n", GENOTYPE_col_i + 1);
				exit(0);
			}

			// Parse the location: Insert location is used as is, unlike deletions.
			char* var_chr = t_string::copy_me_str(cur_line_toks->at(CHROM_col_i)->str());
			int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str());

			// Exclude SNVs.
			int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());

			int eoi_reg_chr_i = t_string::get_i_str(restr_EOI_regs->chr_ids, var_chr);
			if (eoi_reg_chr_i == restr_EOI_regs->chr_ids->size() ||
				(l_alt_allele == 1 &&
				l_ref_allele == 1) ||
				(l_ref_allele >= l_alt_allele)) // Deletion check.

			{
				fprintf(stderr, "Skipping non-ins variant: %s\n", cur_line);
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// We take the first allele in every SNP.
			int l_insertion = -1;
			if (l_ref_allele < l_alt_allele)
			{
				l_insertion = l_alt_allele - l_ref_allele;
			}
			else
			{
				fprintf(stderr, "We were not supposed to be here.\n");
				exit(0);
			}

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "l_insertion: %d; %s\n", l_insertion, cur_line);
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
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
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
								cur_chr_EOI_regs->at(cur_EOI_reg_i)->end >= var_posn+1);

				if (overlap)
				{
					// For each overlap, update ref and alternate alleles.
					int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

					// Set all the positions to the right of this deletion's start to 1.
					insertion_genotype_signal[track_posn] = 1;

					found_overlap = true;

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding Insertion Overlap %s: %s:%d-%d:%s\n", cur_line,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);
					}
				}

				cur_EOI_reg_i++;
			} // region loop.

			if (!found_overlap && __DUMP_CRYPTANNOT_INDEL_MSGS__)
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
		fprintf(stderr, "Saving the insertion signal to %s\n", cur_chr_VCF_signal_op_fp);
		FILE* f_cur_chr_VCF_signal = open_f(cur_chr_VCF_signal_op_fp, "wb");

		// Save the max indel length.
		fwrite(insertion_genotype_signal, sizeof(unsigned int), cur_chr_EOI_covg, f_cur_chr_VCF_signal);

		close_f(f_cur_chr_VCF_signal, cur_chr_VCF_signal_op_fp);
		delete[] insertion_genotype_signal;

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
} // signalize_VCF_Indel_genotypes_per_EOI_regs

// Write the ternary signal on the EOI=gene_exons.
void signalize_VCF_Deletion_genotypes_per_EOI_regs(char* EOI_regs_BED_fp,
													char* per_chrom_VCF_dir,
													char* op_dir)
{
	fprintf(stderr, "Signalizing Deletes under %s over the target regions in %s and saving to %s\n", per_chrom_VCF_dir, EOI_regs_BED_fp, op_dir);

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

		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));

		// Allocat the per indel per l_indel selector.
		fprintf(stderr, "Allocating and initializing the %d long deletion VCF signal.\n", cur_chr_EOI_covg);
		unsigned int* deletion_genotype_signal = new unsigned int[cur_chr_EOI_covg + 2];
		memset(deletion_genotype_signal, 0, sizeof(unsigned int) * cur_chr_EOI_covg);

		// Assign the starting index to each EOI region.
		int cumul_signal_index = 1;
		for (int i_reg = 0; i_reg < cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;
			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "%s:%d-%d @ %d\n",
					cur_chr_EOI_regs->at(i_reg)->name,
					cur_chr_EOI_regs->at(i_reg)->start, cur_chr_EOI_regs->at(i_reg)->end,
					cumul_signal_index);
			}
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
			int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str()) + 1;

			// Exclude SNVs.
			int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());

			int eoi_reg_chr_i = t_string::get_i_str(restr_EOI_regs->chr_ids, var_chr);
			if (eoi_reg_chr_i == restr_EOI_regs->chr_ids->size() ||
				(l_alt_allele == 1 &&
				l_ref_allele == 1) ||
				(l_ref_allele < l_alt_allele)) // Insert check.
			{
				fprintf(stderr, "Skipping non-del variant: %s\n", cur_line);
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// We take the first allele in every SNP.
			int l_deletion = -1;
			if (l_ref_allele > l_alt_allele)
			{
				l_deletion = l_ref_allele - l_alt_allele;
			}
			else
			{
				fprintf(stderr, "We were not supposed to be here.\n");
				exit(0);
			}

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "l_indel: %d; %s\n", l_deletion, cur_line);
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
			int cur_EOI_reg_i = locate_posn_region_per_region_starts(var_posn, cur_chr_EOI_regs, 0, cur_chr_EOI_regs->size());
			while (cur_EOI_reg_i > 0 &&
				cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_end > var_posn)
			{
				cur_EOI_reg_i--;
			}

			// For every overlap, update the allelic signal.
			bool found_overlap = false;
			while (cur_EOI_reg_i < cur_chr_EOI_regs->size() &&
				cur_chr_EOI_regs->at(cur_EOI_reg_i)->sort_info->cumulative_sorted_start <= var_posn)
			{
				int overlap_start = MAX(cur_chr_EOI_regs->at(cur_EOI_reg_i)->start, var_posn);
				int overlap_end = MIN(cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn + l_deletion - 1);
							
				if (overlap_start <= overlap_end)
				{
					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding Deletion Overlap %s: %s:%d-%d:%s ; [%d-%d]\n", cur_line,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							cur_chr_EOI_regs->at(cur_EOI_reg_i)->name, 
							overlap_start, overlap_end);
					}

					// Set all the nucleotides between start and end to 1.
					for (int del_nuc_posn = overlap_start; 
						del_nuc_posn <= overlap_end; 
						del_nuc_posn++)
					{
						int track_posn = cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (del_nuc_posn - cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);
						deletion_genotype_signal[track_posn] = 1;

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "Setting Deletion Overlap @ %d (%d)\n", 
								del_nuc_posn, track_posn);
							getc(stdin);
						}
					} // del_nuc_i loop.

					found_overlap = true;
				} // overlap check.

				cur_EOI_reg_i++;
			} // region loop.

			if (!found_overlap && __DUMP_CRYPTANNOT_INDEL_MSGS__)
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
		fprintf(stderr, "Saving the deletion signal to %s\n", cur_chr_VCF_signal_op_fp);
		FILE* f_cur_chr_VCF_signal = open_f(cur_chr_VCF_signal_op_fp, "wb");

		// Save the max indel length.
		fwrite(deletion_genotype_signal, sizeof(unsigned int), cur_chr_EOI_covg, f_cur_chr_VCF_signal);

		close_f(f_cur_chr_VCF_signal, cur_chr_VCF_signal_op_fp);
		delete[] deletion_genotype_signal;

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
} // signalize_VCF_Indel_genotypes_per_EOI_regs

void multiply_Indel_variant_and_annotation_signals(char* annotation_signal_dir,
													char* variant_signal_dir,
													char* op_dir)
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
		sprintf(cur_chr_annotation_signal_fp, "%s/impact_signal_%s.bin.gz", annotation_signal_dir, var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_annotation_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_annotation_signal_fp, var_chr_ids->at(i_chr));
			continue;
		}

		char cur_chr_variant_signal_fp[1000];
		sprintf(cur_chr_variant_signal_fp, "%s/variant_signal_%s.bin.gz", variant_signal_dir, var_chr_ids->at(i_chr));
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
} // multiply_Indel_variant_and_annotation_signals

void translate_annotated_Deletions_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp)
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

	// Write the header.
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

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "%s:%d-%d @ %d\n",
					sorted_cur_chr_EOI_regs->at(i_reg)->name,
					sorted_cur_chr_EOI_regs->at(i_reg)->start, sorted_cur_chr_EOI_regs->at(i_reg)->end,
					cumul_signal_index);
			}
		} // i_reg loop.
		fprintf(stderr, "Loaded %d EOI regions on %s\n", sorted_cur_chr_EOI_regs->size(), anno_var_chr_ids->at(i_chr));

		// Load the impact signals.
		char cur_chr_anno_variant_signal_fp[1000];
		sprintf(cur_chr_anno_variant_signal_fp, "%s/annotated_variant_signal_%s.bin.gz", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_anno_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_anno_variant_signal_fp, anno_var_chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Loading annotated deletion signals.\n");
		unsigned long long* annotated_deletion_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		FILE* f_annotation_signal = open_f(cur_chr_anno_variant_signal_fp, "rb");
		if (fread(annotated_deletion_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_annotation_signal) != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Could not read annotated deletion signal from %s\n", cur_chr_anno_variant_signal_fp);
			exit(0);
		}
		close_f(f_annotation_signal, cur_chr_anno_variant_signal_fp);

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

		unsigned long long LLONE = 1;

		FILE* f_vcf = open_f(cur_chrom_vcf_fp, "r");
		int CHROM_col_i = 0;
		int POSITION_col_i = 1;
		int VAR_ID_col_i = 2;
		int REF_ALLELE_col_i = 3;
		int ALT_ALLELE_col_i = 4;
		int GENOTYPE_col_i = 9;

		int n_bits_per_nuc = 3;
		int n_neigh_nucs = 6;

		// This is the VCF reading loop, the idea is to read each deletion from the VCF, find the position on the vectorized annotations,
		// then interpret the list of deleted nucleotides into the final annotation.
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
			int var_posn = atoi(cur_line_toks->at(POSITION_col_i)->str()) + 1;
			char* var_name = cur_line_toks->at(VAR_ID_col_i)->str();

			// Exclude SNVs.
			int l_ref_allele = t_string::string_length(cur_line_toks->at(REF_ALLELE_col_i)->str());
			int l_alt_allele = t_string::string_length(cur_line_toks->at(ALT_ALLELE_col_i)->str());
			char* ref_allele_seq = cur_line_toks->at(REF_ALLELE_col_i)->str();
			char* alt_allele_seq = cur_line_toks->at(ALT_ALLELE_col_i)->str();

			int eoi_reg_chr_i = t_string::get_i_str(anno_var_chr_ids, var_chr);
			if (eoi_reg_chr_i == anno_var_chr_ids->size() ||
				(l_alt_allele == 1 &&
				l_ref_allele == 1) ||
				(l_ref_allele < l_alt_allele)) // Insert check.
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// We take the first allele in every SNP.
			int l_deletion = -1;
			if (l_ref_allele > l_alt_allele)
			{
				l_deletion = l_ref_allele - l_alt_allele;
			}
			else
			{
				fprintf(stderr, "We were not supposed to be here.\n");
				exit(0);
			}

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "l_indel: %d; %s\n", l_deletion, cur_line);
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
				int overlap_end = MIN(sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn + l_deletion - 1);

				// Make sure the deletion is engulfed.
				if (overlap_start <= overlap_end &&
					overlap_start <= var_posn &&
					overlap_end >= var_posn + l_deletion - 1)
				{
					fprintf(stderr, "----------------- Tracking Del: %s:%d-%d (%s) -----------------\n",
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding Deletion Overlap %s: %s:%d-%d:%s ; [%d-%d]\n", cur_line,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							overlap_start, overlap_end);
					}

					// For every overlap, update the allelic signal.
					vector<unsigned long long>* cur_var_sorted_annotation_signals = new vector<unsigned long long>();

					// Copy the impact signal values, set the first nucleotide's coding frame, and the junction sequence.
					char* junction_cdna_seq = new char[n_neigh_nucs * 2 + 3];
					//char* original_cdna_seq = new char[n_neigh_nucs * 2 + 3];
					memset(junction_cdna_seq, 0, (n_neigh_nucs * 2 + 3));
					//memset(original_cdna_seq, 0, (n_neigh_nucs * 2 + 3));
					vector<unsigned long long>* all_impacts_values = new vector<unsigned long long>();
					//int first_nuc_frame = 0;

					for (int del_nuc_posn = overlap_start;
						del_nuc_posn <= overlap_end;
						del_nuc_posn++)
					{
						// Copy the annotation signal value.
						int track_posn = sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (del_nuc_posn - sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "Setting Deletion Overlap @ %d (%d)\n",
								del_nuc_posn, track_posn);
							getc(stdin);
						}

						unsigned long long cur_del_annotation_signal_val = annotated_deletion_signal[track_posn];

						// Add the impact value here.
						cur_var_sorted_annotation_signals->push_back(cur_del_annotation_signal_val);

						// Unpack and translate the indel.
						char s_neigh_seq_buff[100];
						unsigned int s_coding_frame = 0;
						unsigned long long int s_impact_val = 0;
						unpack_impact_signal_values(cur_del_annotation_signal_val,
							s_coding_frame, s_impact_val, s_neigh_seq_buff,
							n_bits_per_nuc, n_neigh_nucs);

						// Add the impact value.
						all_impacts_values->push_back(s_impact_val);

						// Following blocks copy the original and post-del nucleotide sequences at the beginning of the deletion with respect to the strand of the feature.
						if(del_nuc_posn == overlap_start)
						{
							for (int nuc_i = 0; nuc_i < n_neigh_nucs / 2; nuc_i++)
							{
								junction_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
							} // nuc_i loop
						}

						if (del_nuc_posn == overlap_end)
						{
							for (int nuc_i = n_neigh_nucs / 2; nuc_i < n_neigh_nucs; nuc_i++)
							{
								junction_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
							} // nuc_i loop
						}
					} // del_nuc_i loop.

					int first_nuc_frame = 0;
					int first_tracked_nuc_frame = -1;
					char* deleted_CDS_nucs = track_original_del_neigh_sequence_per_overlap_posns(sorted_cur_chr_EOI_regs, cur_EOI_reg_i,
																								first_tracked_nuc_frame,
																								ref_allele_seq,
																								annotated_deletion_signal,
																								overlap_start, overlap_end,
																								n_bits_per_nuc, n_neigh_nucs,
																								vep_term_ctx);

					first_nuc_frame = first_tracked_nuc_frame;
					char* original_cdna_seq = t_string::copy_me_str(deleted_CDS_nucs);

					fprintf(stderr, "Extracted original del seq: %s (Frame: %d)\n", deleted_CDS_nucs, first_nuc_frame);

					// Concatenate the first and last sequences.
					// If the annotation is on the negative strand, reverse-complement the junction sequence.
					if (t_string::string_length(original_cdna_seq) > 0 &&
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand == '-')
					{
						reverse_complement_seq(junction_cdna_seq);
						reverse_complement_seq(original_cdna_seq);
					}

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						// Write the codon-by-codon junction sequence:
						char* new_codon = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);
						char* prev_codon = (char*)(&original_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);

						// Get the alternate codon.
						int first_nuc_frame_alt = ((first_nuc_frame + 1) % 3);
						char* new_codon_alt = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame_alt]);

						fprintf(stderr, "Original sequence frames:\n");
						int prev_i = 0;
						while (1)
						{
							if (prev_codon[prev_i] > 0 &&
								prev_codon[prev_i + 1] > 0 &&
								prev_codon[prev_i + 2] > 0)
							{
								fprintf(stderr, "|%c%c%c", prev_codon[prev_i], prev_codon[prev_i + 1], prev_codon[prev_i + 2]);
							}
							else
							{
								break;
							}
							prev_i += 3;
						} // prev_i loop.
						fprintf(stderr, "\n");

						fprintf(stderr, "Deleted sequence frames:\n");
						int new_i = 0;
						while (1)
						{
							if (new_codon[new_i] > 0 &&
								new_codon[new_i + 1] > 0 &&
								new_codon[new_i + 2] > 0)
							{
								fprintf(stderr, "|%c%c%c", new_codon[new_i], new_codon[new_i + 1], new_codon[new_i + 2]);
							}
							else
							{
								break;
							}
							new_i += 3;
						} // new_i loop.
						fprintf(stderr, "\n");
					}

					//////////////////////////////////////////////////////////////////////////////
					// We have setup the sequences and frames, interpret the impact.
					//////////////////////////////////////////////////////////////////////////////
					// Check for the stop gains: The first nucleotide on the 4th position takes the first nuc frame.
					bool start_loss = false;
					bool start_retain = false;
					bool stop_loss = false;
					bool stop_retained = false;
					bool stop_gain = false; // Stop gain is important only when there is an inframe insertion.

					// Set the AAs.
					char* prev_AA = NULL;
					char* new_AA = NULL;
					char* new_AA_alt = NULL;
					if (t_string::string_length(original_cdna_seq) > 0)
					{
						char* new_codon = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);
						char* prev_codon = (char*)(&original_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame]);

						// Get the alternate codon.
						int first_nuc_frame_alt = ((first_nuc_frame + 1) % 3);
						char* new_codon_alt = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_nuc_frame_alt]);

						new_AA = get_AA_code_per_codon(new_codon);
						prev_AA = get_AA_code_per_codon(prev_codon);
						new_AA_alt = get_AA_code_per_codon(new_codon_alt);
					}
					else
					{
						prev_AA = t_string::copy_me_str("-");
						new_AA = t_string::copy_me_str("-");
						new_AA_alt = t_string::copy_me_str("-");
					}

					fprintf(stderr, "%s:%d-%d: %s->%s (%s) :: %s->%s (%s) @ Frame: %d\n",
							anno_var_chr_ids->at(i_chr),
							var_posn,
							var_posn + l_deletion - 1,
							original_cdna_seq, junction_cdna_seq, deleted_CDS_nucs,
							prev_AA, new_AA, new_AA_alt,
							first_nuc_frame);

					// Check and set the stop losses.
					if (t_string::compare_strings(prev_AA, "STOP"))
					{
						if (t_string::compare_strings(new_AA, "STOP"))
						{
							stop_retained = true;
						}
						else
						{
							stop_loss = true;
						}
					}
					else
					{
						if (t_string::compare_strings(new_AA, "STOP"))
						{
							stop_gain = true;
						}
					}

					// Count the number of deleted coding nucleotides.
					int n_del_coding_nucs = 0;

					int n_deleted_coding_nucs = 0; int n_impacted_CDSs = 0;
					int n_del_splice_region_nucs = 0; int n_impacted_splice_regs = 0;
					int n_del_splice_acc_nucs = 0; int n_impacted_splice_acc_sites = 0;
					int n_del_splice_don_nucs = 0; int n_impacted_splice_don_sites = 0;
					int n_del_start_codon_nucs = 0; int n_impacted_start_codon_sites = 0;
					int n_del_stop_codon_nucs = 0; int n_impacted_stop_codon_sites = 0;
					int n_del_intron_nucs = 0; int n_impacted_introns = 0;
					int n_del_fp_utr_nucs = 0; int n_impacted_fp_utrs = 0;
					int n_del_tp_utr_nucs = 0; int n_impacted_tp_utrs = 0;

					// Count the number of nucleotides and elements.
					count_deleted_elements_left2right_impact_states(all_impacts_values,
						vep_term_ctx,
						n_deleted_coding_nucs, n_impacted_CDSs,
						n_del_splice_region_nucs, n_impacted_splice_regs,
						n_del_splice_acc_nucs, n_impacted_splice_acc_sites,
						n_del_splice_don_nucs, n_impacted_splice_don_sites,
						n_del_start_codon_nucs, n_impacted_start_codon_sites,
						n_del_stop_codon_nucs, n_impacted_stop_codon_sites,
						n_del_intron_nucs, n_impacted_introns,
						n_del_fp_utr_nucs, n_impacted_fp_utrs,
						n_del_tp_utr_nucs, n_impacted_tp_utrs);

					fprintf(stderr, "%s/%s: CDS: %d/%d\n\
Splice_Reg: %d/%d\n\
Splice_Acc: %d/%d\n\
Splice_Don: %d/%d\n\
Start_Codon: %d/%d\n\
Stop_Codon: %d/%d\n\
Intron: %d/%d\n\
5'-UTR_Codon: %d/%d\n\
3'-UTR_Codon: %d/%d\n", var_name, sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name, 
n_deleted_coding_nucs, n_impacted_CDSs,
n_del_splice_region_nucs, n_impacted_splice_regs,
n_del_splice_acc_nucs, n_impacted_splice_acc_sites,
n_del_splice_don_nucs, n_impacted_splice_don_sites,
n_del_start_codon_nucs, n_impacted_start_codon_sites,
n_del_stop_codon_nucs, n_impacted_stop_codon_sites,
n_del_intron_nucs, n_impacted_introns,
n_del_fp_utr_nucs, n_impacted_fp_utrs,
n_del_tp_utr_nucs, n_impacted_tp_utrs);

					int n_total_impacting_nucs = n_del_start_codon_nucs + n_del_stop_codon_nucs +
						n_del_splice_don_nucs + n_del_splice_acc_nucs + n_del_splice_region_nucs +
						n_del_intron_nucs + 
						n_deleted_coding_nucs +
						n_del_fp_utr_nucs + n_del_tp_utr_nucs;					

					// Depending on the number of nucs deleted, report the string.
					t_string* VEP_impact_str = new t_string();

					if (n_total_impacting_nucs == 0)
					{
						VEP_impact_str->copy("IMPACT_NOMINAL;");
						strcpy(prev_AA, "-");
						strcpy(new_AA, "-");
					}
					else
					{
						if (n_deleted_coding_nucs > 0)
						{
							if (n_deleted_coding_nucs % 3 == 0)
							{
								VEP_impact_str->concat_string("inframe_deletion;");
							}
							else
							{
								VEP_impact_str->concat_string("frameshift_variant;");
							}
						}

						if (n_del_splice_region_nucs > 0)
						{
							VEP_impact_str->concat_string("splice_region;");
						}

						if (n_del_splice_acc_nucs > 0)
						{
							VEP_impact_str->concat_string("splice_acceptor_variant;");
						}

						if (n_del_splice_don_nucs > 0)
						{
							VEP_impact_str->concat_string("splice_donor_variant;");
						}

						if (n_del_start_codon_nucs > 0)
						{
							VEP_impact_str->concat_string("start_lost;");

							// Check for start-retain: Check if this is an M variant.
							if (prev_AA == new_AA)
							{
								VEP_impact_str->concat_string("start_retained;");
							}
						}

						if (n_del_stop_codon_nucs > 0)
						{
							VEP_impact_str->concat_string("stop_lost;");

							if (t_string::compare_strings(prev_AA, new_AA) &&
								t_string::compare_strings(prev_AA, "STOP"))
							{
								VEP_impact_str->concat_string("stop_retained_variant;");
							}
						}

						if (n_del_intron_nucs > 0)
						{
							VEP_impact_str->concat_string("intron_variant;");
						}

						if (n_del_fp_utr_nucs > 0)
						{
							VEP_impact_str->concat_string("5_prime_UTR_variant;");
						}

						if (n_del_tp_utr_nucs > 0)
						{
							VEP_impact_str->concat_string("3_prime_UTR_variant;");
						}
					} // impact string building.

					// Stop gain needs to be explicitly checked.
					if (n_deleted_coding_nucs > 0 &&
						stop_gain)
					{
						VEP_impact_str->concat_string("stop_gained;");
					}

					// Sanity check.
					if(first_nuc_frame >= 3)
					{
						fprintf(stderr, "Sanity check failed: The nucleotide frame is out of range: %d\n", first_nuc_frame);
						exit(0);
					}

					// Write the translated annotation.
					fprintf(f_op, "%s\t%d\t%d\t%s\t%s\t%c\t%s\t%s\t%s\t%s\n",
							anno_var_chr_ids->at(i_chr),
							var_posn,
							var_posn + l_deletion - 1,
							var_name, 
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand,
							prev_AA, new_AA, new_AA_alt,
							VEP_impact_str->str());

					delete VEP_impact_str;

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Sorted Impacts @ %s:%d-%d on %s (%c):\n",
								anno_var_chr_ids->at(i_chr),
								var_posn, 
								var_posn + l_deletion - 1,
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand);

						for (int i_impact = 0; i_impact < cur_var_sorted_annotation_signals->size(); i_impact++)
						{
							// Unpack and translate the indel.
							char s_neigh_seq_buff[100];
							unsigned int s_coding_frame = 0;
							unsigned long long s_impact_val = 0;
							unpack_impact_signal_values(cur_var_sorted_annotation_signals->at(i_impact),
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
					delete cur_var_sorted_annotation_signals;

					found_overlap = true;

					fprintf(stderr, "----------------- End New Del. Tracking -----------------\n");
				} // overlap check.

				cur_EOI_reg_i++;
			} // forward searching loop.

			if (!found_overlap && __DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "Could not find an overlap for %s\n", cur_line);
			}

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			delete[] cur_line;
		} // file reading loop.
		close_f(f_vcf, cur_chrom_vcf_fp);

		fprintf(stderr, "\nFinished reading VCF file.\n");
	} // i_chr loop.

	  // Close output file.
	close_f(f_op, op_fp);
}

void translate_annotated_Insertions_from_annotated_signals(char* annotated_variant_signal_dir, char* per_chrom_VCF_dir, char* op_fp)
{
	fprintf(stderr, "Interpreting the annotated insertion signals from %s to %s.\n", annotated_variant_signal_dir, op_fp);

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

	// Write the header.
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

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "%s:%d-%d @ %d\n",
					sorted_cur_chr_EOI_regs->at(i_reg)->name,
					sorted_cur_chr_EOI_regs->at(i_reg)->start, sorted_cur_chr_EOI_regs->at(i_reg)->end,
					cumul_signal_index);
			}
		} // i_reg loop.
		fprintf(stderr, "Loaded %d EOI regions on %s\n", sorted_cur_chr_EOI_regs->size(), anno_var_chr_ids->at(i_chr));

		// Load the impact signals.
		char cur_chr_anno_variant_signal_fp[1000];
		sprintf(cur_chr_anno_variant_signal_fp, "%s/annotated_variant_signal_%s.bin.gz", annotated_variant_signal_dir, anno_var_chr_ids->at(i_chr));
		if (!check_file(cur_chr_anno_variant_signal_fp))
		{
			fprintf(stderr, "Could not find %s for %s, skipping.\n", cur_chr_anno_variant_signal_fp, anno_var_chr_ids->at(i_chr));
			exit(0);
		}

		fprintf(stderr, "Loading annotated insertion signals.\n");
		unsigned long long* annotated_insertion_signal = new unsigned long long[cur_chr_EOI_covg + 2];
		FILE* f_annotation_signal = open_f(cur_chr_anno_variant_signal_fp, "rb");
		if(fread(annotated_insertion_signal, sizeof(unsigned long long), cur_chr_EOI_covg, f_annotation_signal) != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Could not read annotated insertion signal from %s\n", cur_chr_anno_variant_signal_fp);
			exit(0);
		}

		close_f(f_annotation_signal, cur_chr_anno_variant_signal_fp);

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
			fprintf(stderr, "Reading inserts in %s\n", cur_chrom_vcf_fp);
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
				(l_alt_allele == 1 && l_ref_allele == 1) ||
					(l_ref_allele > l_alt_allele)) // Deletion check.
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				continue;
			}

			// We take the first allele in every SNP.
			int l_insertion = -1;
			if (l_ref_allele < l_alt_allele)
			{
				l_insertion = l_alt_allele - l_ref_allele;
			}
			else
			{
				fprintf(stderr, "We were not supposed to be here.\n");
				exit(0);
			}

			if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "l_indel: %d; %s\n", l_insertion, cur_line);
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
				int overlap_end = MIN(sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end, var_posn + 1);

				// Make sure the insertion is engulfed.
				if (overlap_start <= overlap_end &&
					overlap_start <= var_posn &&
					overlap_end >= var_posn + 1)
				{
					fprintf(stderr, "----------------- Tracking Ins.: %s:%d-%d (%s) -----------------\n",
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name);

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Adding Insertion Overlap %s: %s:%d-%d:%s ; [%d-%d]\n", cur_line,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->chrom,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->end,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							overlap_start, overlap_end);
					}

					// For every overlap, update the allelic signal.
					vector<unsigned long long>* cur_var_sorted_impacts = new vector<unsigned long long>();

					// Copy the impact signal values, set the first nucleotide's coding frame, and the junction sequence.
					char* junction_cdna_seq = new char[2*n_neigh_nucs + l_insertion + 3];
					char* original_cdna_seq = new char[2 * n_neigh_nucs * l_insertion + 3];
					memset(junction_cdna_seq, 0, (2 * n_neigh_nucs + l_insertion + 3));
					memset(original_cdna_seq, 0, (2 * n_neigh_nucs + l_insertion + 3));
					int first_ins_nuc_frame = 0;

					unsigned long long insertion_impact = 0;

					// Copy the annotation signal value.
					int track_posn = sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->score + (var_posn - sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->start);

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Setting Insertion Overlap @ %d (%d)\n",
							var_posn, track_posn);
						getc(stdin);
					}

					unsigned long long cur_ins_annotation_signal_val = annotated_insertion_signal[track_posn];

					// Unpack and translate the indel.
					char s_neigh_seq_buff[100];
					unsigned int s_coding_frame = 0;
					unsigned long long s_impact_val = 0;
					unpack_impact_signal_values(cur_ins_annotation_signal_val,
						s_coding_frame, s_impact_val, s_neigh_seq_buff,
						n_bits_per_nuc, n_neigh_nucs);

					// Set the insertion's impact.
					insertion_impact = s_impact_val;
					// Sequence copying: For the 
					{
						// Copy the neighborhood sequence around the deleted nuc, including the reference sequence.
						for (int nuc_i = 0; nuc_i < n_neigh_nucs / 2; nuc_i++)
						{
							original_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
						} // nuc_i loop

						for (int nuc_i = n_neigh_nucs / 2; nuc_i < n_neigh_nucs; nuc_i++)
						{
							original_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
						} // nuc_i loop

						// Concat the left neighbor sequence.
						for (int nuc_i = 0; nuc_i < n_neigh_nucs / 2; nuc_i++)
						{
							junction_cdna_seq[nuc_i] = s_neigh_seq_buff[nuc_i];
						} // nuc_i loop

						// Concat the inserted sequence: Must skip the first nucleotide, which is included.
						strcat(junction_cdna_seq, &alt_allele_seq[1]);
						int l_junc = strlen(junction_cdna_seq);

						// Concat the right neighbor sequence.
						for (int nuc_i = n_neigh_nucs / 2; nuc_i < n_neigh_nucs; nuc_i++)
						{
							junction_cdna_seq[l_junc + nuc_i - n_neigh_nucs / 2] = s_neigh_seq_buff[nuc_i];
						} // nuc_i loop
					}

					// This is the frame of the first inserted nucleotide among the inserted nucleotides, **taking the strand into account since VEP takes strand into account while reporting CDS position.**.
					first_ins_nuc_frame = (s_coding_frame + 1) % 3;

					// Concatenate the first and last sequences.
					// If the annotation is on the negative strand, reverse-complement the junction sequence.
					if (sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand == '-')
					{
						reverse_complement_seq(junction_cdna_seq);
						reverse_complement_seq(original_cdna_seq);
					}

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						// Write the codon-by-codon junction sequence:
						char* new_codon = (char*)(&junction_cdna_seq[n_neigh_nucs / 2 - first_ins_nuc_frame]);
						char* prev_codon = (char*)(&original_cdna_seq[n_neigh_nucs / 2 - first_ins_nuc_frame]);

						fprintf(stderr, "Original sequence frames:\n");
						int prev_i = 0;
						while (1)
						{
							if (prev_codon[prev_i] > 0 &&
								prev_codon[prev_i + 1] > 0 &&
								prev_codon[prev_i + 2] > 0)
							{
								fprintf(stderr, "|%c%c%c", prev_codon[prev_i], prev_codon[prev_i + 1], prev_codon[prev_i + 2]);
							}
							else
							{
								break;
							}
							prev_i += 3;
						} // prev_i loop.
						fprintf(stderr, "\n");

						fprintf(stderr, "Inserted sequence frames:\n");
						int new_i = 0;
						while (1)
						{
							if (new_codon[new_i] > 0 &&
								new_codon[new_i + 1] > 0 &&
								new_codon[new_i + 2] > 0)
							{
								fprintf(stderr, "|%c%c%c", new_codon[new_i], new_codon[new_i + 1], new_codon[new_i + 2]);
							}
							else
							{
								break;
							}
							new_i += 3;
						} // new_i loop.
						fprintf(stderr, "\n");
					}

					//////////////////////////////////////////////////////////////////////////////
					// We have setup the sequences and frames, interpret the impact.
					//////////////////////////////////////////////////////////////////////////////

					// Check for the stop gains: The first nucleotide on the 4th position takes the first nuc frame.					
					bool start_loss = false;
					bool start_retain = false;
					bool stop_loss = false;
					bool stop_retained = false;
					bool stop_gain = false; // Stop gain is important only when there is an inframe insertion.
					bool inframe_ins = false;

					if (l_insertion % 3 == 0)
					{
						// Set the in frame insertion.
						inframe_ins = true;
					}
					
					if (first_ins_nuc_frame >= 3)
					{
						fprintf(stderr, "Sanity check failed: The nucleotide frame is out of range: %d\n", first_ins_nuc_frame);
						exit(0);
					}

					// If the insertion is in-frame, check the stop gain.
					char* prev_codon = (char*)(&original_cdna_seq[n_neigh_nucs / 2 - first_ins_nuc_frame]);
					char* prev_AA = get_AA_code_per_codon(prev_codon);
					
					// Check for stop retain/gain.
					bool found_stop = false;
					if (t_string::compare_strings(prev_AA, "STOP"))
					{
						stop_loss = true;
					}

					// Stop-retainment/gain check block.
					int inframe_nuc_i = n_neigh_nucs / 2 - first_ins_nuc_frame;

					fprintf(stderr, "Checking stop-retention/gain for insertion for %s:%d :: junction_seq: %s on %s(first_ins_nuc_frame: %d)\n", var_name, var_posn, junction_cdna_seq, sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name, first_ins_nuc_frame);

					// Go over all the codons and 
					for (int i = inframe_nuc_i; 
						(i+2) < strlen(junction_cdna_seq); 
						i += 3)
					{
						char* cur_codon = (char*)(&junction_cdna_seq[i]);
						char* new_AA = get_AA_code_per_codon(cur_codon);

						fprintf(stderr, "%d: %c%c%c->%s\n", i, cur_codon[0], cur_codon[1], cur_codon[2], new_AA);

						if (t_string::compare_strings(new_AA, "STOP"))
						{
							if (t_string::compare_strings(prev_AA, "STOP"))
							{
								stop_retained = true;
							}
							else
							{
								stop_gain = true;
							}
							found_stop = true;
							break;
						}

						if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
						{
							fprintf(stderr, "%s::%s:%d-%d (%c): Inserted AA: %s (%d):%s <-> %s (%s)\n",
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
								anno_var_chr_ids->at(i_chr),
								var_posn, var_posn + 1,
								sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand,
								junction_cdna_seq, i,
								cur_codon, new_AA,
								prev_AA);
						}
					} // i loop.						

					fprintf(stderr, "%s::%s:%d-%d @ %d;; %llX\n",
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
						anno_var_chr_ids->at(i_chr),
						var_posn, var_posn + 1, track_posn, insertion_impact);

					//// Report using the VEP terms.
					//// Get the unique values of impacts.
					//vector<unsigned long long>* all_impacts_values = new vector<unsigned long long>();
					//all_impacts_values->push_back(insertion_impact);

					// Depending on where the nucleotide is inserted, report the impact.
					t_string* vep_impact_str = new t_string();
					if (check_element_context_per_impact_bitmap(vep_term_ctx, CDS_ctx, insertion_impact))
					{
						// This is a coding context insertion.
						if (l_insertion % 3 == 0)
						{
							vep_impact_str->concat_string("inframe_insertion;");
						}
						else
						{
							vep_impact_str->concat_string("frameshift_variant;");
						}

						// Stop gain is check under CDS context.
						if (stop_gain)
						{
							vep_impact_str->concat_string("stop_gained;");
						}
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, fp_UTR_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("5_prime_UTR_variant;");
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, intron_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("intron_variant;");
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, tp_UTR_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("3_prime_UTR_variant;");
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, start_codon_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("start_lost;");
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, stop_codon_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("stop_lost;");

						if (stop_retained)
						{
							vep_impact_str->concat_string("stop_retained_variant;");
						}
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_region_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("splice_region_variant;");
					}

					if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_acceptor_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("splice_acceptor_variant;");
					}
					
					if (check_element_context_per_impact_bitmap(vep_term_ctx, splice_donor_ctx, insertion_impact))
					{
						vep_impact_str->concat_string("splice_donor_variant;");
					}

					if (vep_impact_str->length() == 0)
					{
						vep_impact_str->concat_string("NOMINAL_IMPACT;");
					}

					/*if (!CDS_impact &&
						!start_codon_impact &&
						!stop_codon_impact)
					{
						prev_AA = t_string::copy_me_str("-");

						inframe_ins = false;
					}

					if (!stop_codon_impact)
					{
						stop_loss = false;
					}*/

					// Write the translated annotation.
					fprintf(f_op, "%s\t%d\t%d\t%s\t%s\t%c\t%s\t%s\t%s\t%s\n",
						anno_var_chr_ids->at(i_chr),
						var_posn,
						var_posn + l_insertion - 1,
						var_name,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
						sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand,
						prev_AA, ".", ".",
						vep_impact_str->str());

					delete vep_impact_str;

					if (__DUMP_CRYPTANNOT_INDEL_MSGS__)
					{
						fprintf(stderr, "Sorted Impacts @ %s:%d-%d on %s (%c):\n",
							anno_var_chr_ids->at(i_chr),
							var_posn,
							var_posn + l_insertion - 1,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->name,
							sorted_cur_chr_EOI_regs->at(cur_EOI_reg_i)->strand);

						for (int i_impact = 0; i_impact < cur_var_sorted_impacts->size(); i_impact++)
						{
							// Unpack and translate the indel.
							char s_neigh_seq_buff[100];
							unsigned int s_coding_frame = 0;
							unsigned long long s_impact_val = 0;
							unpack_impact_signal_values(cur_var_sorted_impacts->at(i_impact),
								s_coding_frame, s_impact_val, s_neigh_seq_buff,
								n_bits_per_nuc, n_neigh_nucs);

							fprintf(stderr, "%llX -- %s @ frame %d\n",
								s_impact_val,
								s_neigh_seq_buff,
								s_coding_frame);
						} // i_impact loop.

						fprintf(stderr, "Junction cDNA Seq: %s\n", junction_cdna_seq);

						getc(stdin);
					} // dumping check.
					delete cur_var_sorted_impacts;

					found_overlap = true;

					fprintf(stderr, "----------------- End New Ins. Tracking -----------------\n");
				} // overlap check.

				cur_EOI_reg_i++;
			} // region loop.

			if (!found_overlap && __DUMP_CRYPTANNOT_INDEL_MSGS__)
			{
				fprintf(stderr, "Could not find an overlap for %s\n", cur_line);
			}

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			delete[] cur_line;
		} // file reading loop.
		close_f(f_vcf, cur_chrom_vcf_fp);

		fprintf(stderr, "\nFinished reading VCF file.\n");
	} // i_chr loop.

	  // Close output file.
	close_f(f_op, op_fp);
}


