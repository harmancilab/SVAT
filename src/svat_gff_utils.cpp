#include <stdio.h>
#include <stdlib.h>
#include "svat_annot_region_tools.h"
#include <algorithm>
#include <string.h>
#include "file_utils.h"
#include "svat_nomenclature.h"
#include "svat_ansi_string.h"
#include "svat_genomics_coords.h"
#include "svat_genome_sequence_tools.h"
#include "svat_gff_utils.h"


void load_elements_set_annotation_information_per_GENCODE_GFF(char* gencode_gff_fp, char* genome_bin_seq_dir,
																vector<t_annot_region*>* gene_regs,
																vector<t_annot_region*>* transcript_regs,
																vector<t_annot_region*>* cds_regs,														
																vector<t_annot_region*>* UTR_5p_regs,
																vector<t_annot_region*>* UTR_3p_regs,
																vector<t_annot_region*>* intronic_regs,
																vector<t_annot_region*>* promotor_regs, int l_promotor)
{
	// Load the gff file: Load CDSs for each transcript.
	vector<t_annot_region*>* gff_regs = load_GFF(gencode_gff_fp);
	fprintf(stderr, "Loaded %d GFF regions from %s, parsing the elements.\n", (int)(gff_regs->size()), gencode_gff_fp);

	vector<t_annot_region*>* exon_regs = new vector<t_annot_region*>();
	for(int i_reg = 0; i_reg < gff_regs->size(); i_reg++)
	{
		t_gff_info* gff_info = (t_gff_info*)(gff_regs->at(i_reg)->data);

		// Parse the transcript name:
		t_string_tokens* cur_reg_toks = t_string::tokenize_by_chars(gff_info->group_str, ";=");
		char transcript_name[100];
		char gene_name[100];
		
		for(int i_t = 0; i_t < cur_reg_toks->size(); i_t++)
		{
			if(t_string::compare_strings(cur_reg_toks->at(i_t)->str(), "transcript_name"))
			{
				strcpy(transcript_name, cur_reg_toks->at(i_t+1)->str());
			}

			if(t_string::compare_strings(cur_reg_toks->at(i_t)->str(), "gene_name"))
			{
				strcpy(gene_name, cur_reg_toks->at(i_t+1)->str());
			}
		} // i_t loop.
		t_string::clean_tokens(cur_reg_toks);

		// If this is an exonic region, add to the list of exons.
		if(t_string::compare_strings(gff_info->feature_type, "exon"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			exon_regs->push_back(dup_reg);
		}

		if(t_string::compare_strings(gff_info->feature_type, "CDS"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(transcript_name);
			cds_regs->push_back(dup_reg);

			t_CDS_annotation_info* annot_info = new t_CDS_annotation_info();
			dup_reg->annotation_info = annot_info;
		} // CDS check
		else if(t_string::compare_strings(gff_info->feature_type, "transcript"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(transcript_name);
			transcript_regs->push_back(dup_reg);

			t_transcript_annotation_info* annot_info = new t_transcript_annotation_info();
			annot_info->gene_name = t_string::copy_me_str(gene_name);
			annot_info->CDSs = new vector<t_annot_region*>();
			annot_info->valid_3mer_frame = true;
			annot_info->l_transcript = 0;
			dup_reg->annotation_info = annot_info;

			// Add a promotor.
			t_annot_region* promotor_reg = duplicate_region(gff_regs->at(i_reg));
			promotor_reg->name = t_string::copy_me_str(transcript_name);
			if(promotor_reg->strand == '+')
			{
				promotor_reg->start = (gff_regs->at(i_reg)->start - l_promotor);
				promotor_reg->end = (gff_regs->at(i_reg)->start + l_promotor);
			}
			else
			{
				promotor_reg->start = (gff_regs->at(i_reg)->end - l_promotor);
				promotor_reg->end = (gff_regs->at(i_reg)->end + l_promotor);
			}
			promotor_regs->push_back(promotor_reg);
		} // transcript check.
		else if(t_string::compare_strings(gff_info->feature_type, "gene"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(gene_name);
			gene_regs->push_back(dup_reg);

			t_gene_annotation_info* annot_info = new t_gene_annotation_info();
			dup_reg->annotation_info = annot_info;
		} // gene chceck.
		else if (t_string::compare_strings(gff_info->feature_type, "UTR"))
		{
			// Get 5' UTR.
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(transcript_name);

			UTR_5p_regs->push_back(dup_reg);
		} // UTR check.
		else if(t_string::compare_strings(gff_info->feature_type, "five_prime_UTR"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(transcript_name);

			UTR_5p_regs->push_back(dup_reg);
		} // 5' UTR check.
		else if (t_string::compare_strings(gff_info->feature_type, "three_prime_UTR"))
		{
			t_annot_region* dup_reg = duplicate_region(gff_regs->at(i_reg));
			dup_reg->name = t_string::copy_me_str(transcript_name);

			UTR_3p_regs->push_back(dup_reg);
		} // 3' UTR check.
	} // i_reg loop.

	// Get the intronic regions: The intronic regions exclude the regions where we see an exon for ANY gene.
	fprintf(stderr, "Generating the intronic regions.\n");
	vector<t_annot_region*>* excluded_intronic_regions = subtract_annot_regions(gene_regs, exon_regs);
	delete_annot_regions(exon_regs);
	for(int i_int = 0; i_int < excluded_intronic_regions->size(); i_int++)
	{
		t_annot_region* cur_intronic_reg = duplicate_region(excluded_intronic_regions->at(i_int));
		t_intron_annotation_info* intron_info = new t_intron_annotation_info();
		intron_info->gene = (t_annot_region*)(excluded_intronic_regions->at(i_int)->data);
		intron_info->gene_name = (intron_info->gene->name);
		cur_intronic_reg->name = t_string::copy_me_str(intron_info->gene->name);
		cur_intronic_reg->annotation_info = intron_info;
		intronic_regs->push_back(cur_intronic_reg);
	} // i_int loop.
	delete excluded_intronic_regions;
	fprintf(stderr, "Extracted %d intronic regions.\n", intronic_regs->size());

	// Sort the transcript names.
	sort(transcript_regs->begin(), transcript_regs->end(), sort_transcripts_regions_per_name);

	vector<char*>* sorted_transcript_names = new vector<char*>();
	for(int i_r = 0; i_r < transcript_regs->size(); i_r++)
	{
		sorted_transcript_names->push_back(transcript_regs->at(i_r)->name);
	} // i_r loop.

	// Sort the genes per names.
	sort(gene_regs->begin(), gene_regs->end(), sort_transcripts_regions_per_name);

	vector<char*>* sorted_gene_names = new vector<char*>();
	for(int i_r = 0; i_r < gene_regs->size(); i_r++)
	{
		sorted_gene_names->push_back(gene_regs->at(i_r)->name);
	} // i_r loop.

	// Assign the CDSs to the transcripts.
	fprintf(stderr, "Assigning the CDS's to the transcripts.\n");
	for(int i_cds = 0; i_cds < cds_regs->size(); i_cds++)
	{
		t_CDS_annotation_info* cur_cds_info = (t_CDS_annotation_info*)(cds_regs->at(i_cds)->annotation_info);

		cur_cds_info->transcript = NULL;

		int trans_i = t_string::fast_search_string_per_prefix(cds_regs->at(i_cds)->name, sorted_transcript_names, 0, sorted_transcript_names->size()-1);

		while(trans_i > 0 &&
			trans_i < sorted_transcript_names->size() &&
			(t_string::sort_strings_per_prefix(cds_regs->at(i_cds)->name, sorted_transcript_names->at(trans_i)) ||
			t_string::compare_strings(cds_regs->at(i_cds)->name, sorted_transcript_names->at(trans_i))))
		{
			trans_i--;
		}

		while(trans_i >= 0 &&
			trans_i < sorted_transcript_names->size() &&
			(t_string::sort_strings_per_prefix(sorted_transcript_names->at(trans_i), cds_regs->at(i_cds)->name) ||
			t_string::compare_strings(cds_regs->at(i_cds)->name, sorted_transcript_names->at(trans_i))))
		{
			if(t_string::compare_strings(cds_regs->at(i_cds)->name, sorted_transcript_names->at(trans_i)))
			{
				vector<t_annot_region*>* transcript_cds_regs = ((t_transcript_annotation_info*)(transcript_regs->at(trans_i)->annotation_info))->CDSs;
				if(!t_string::compare_strings(cds_regs->at(i_cds)->name, transcript_regs->at(trans_i)->name))
				{
					fprintf(stderr, "Match does not hold: %s; %s, %s\n",  cds_regs->at(i_cds)->name, 
							sorted_transcript_names->at(trans_i),
							transcript_regs->at(trans_i)->name);
					exit(0);
				}

				cur_cds_info->transcript = transcript_regs->at(trans_i);

				transcript_cds_regs->push_back(cds_regs->at(i_cds));
				break;
			}

			trans_i++;
		}

		if(cur_cds_info->transcript == NULL)
		{
			fprintf(stderr, "Could not assign the transcript for %s (%d)\n", cds_regs->at(i_cds)->name, trans_i, sorted_transcript_names->at(trans_i));
			getc(stdin);
		}
	} // i_cds loop.

	fprintf(stderr, "Assigning the transcripts to the genes.\n");
	for(int i_tr = 0; i_tr < transcript_regs->size(); i_tr++)
	{
		t_transcript_annotation_info* cur_tr_annot_info = (t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info);

		int gene_i = t_string::fast_search_string_per_prefix(cur_tr_annot_info->gene_name, sorted_gene_names, 0, sorted_gene_names->size()-1);
		cur_tr_annot_info->gene = NULL;

		while(gene_i > 0 &&
			gene_i < sorted_gene_names->size() &&
			(t_string::sort_strings_per_prefix(cur_tr_annot_info->gene_name, sorted_gene_names->at(gene_i)) ||
			t_string::compare_strings(cur_tr_annot_info->gene_name, sorted_gene_names->at(gene_i))))
		{
			gene_i--;
		}

		while(gene_i > 0 &&
			gene_i < sorted_gene_names->size() &&
			(t_string::sort_strings_per_prefix(sorted_gene_names->at(gene_i), cur_tr_annot_info->gene_name) ||
			t_string::compare_strings(cur_tr_annot_info->gene_name, sorted_gene_names->at(gene_i))))
		{
			if(t_string::compare_strings(cur_tr_annot_info->gene_name, sorted_gene_names->at(gene_i)))
			{
				if(!t_string::compare_strings(cur_tr_annot_info->gene_name, gene_regs->at(gene_i)->name))
				{
					fprintf(stderr, "Match does not hold: %s; %s, %s\n", cur_tr_annot_info->gene_name, 
						sorted_gene_names->at(gene_i),
						gene_regs->at(gene_i)->name);
					exit(0);
				}

				cur_tr_annot_info->gene = gene_regs->at(gene_i);

				break;
			}

			gene_i++;
		}

		if(cur_tr_annot_info->gene == NULL)
		{
			fprintf(stderr, "Could not assign gene for %s\n", transcript_regs->at(i_tr)->name);
		}
	} // i_tr loop.

	fprintf(stderr, "Extracting the sequences.\n");
	vector<char*>* cds_chr_ids = get_chr_ids(cds_regs);
	for(int i_chr = 0; i_chr < cds_chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", cds_chr_ids->at(i_chr));
		int l_chr_seq = 0;
		char bin_seq_fp[1000];
		sprintf(bin_seq_fp, "%s/%s.bin", genome_bin_seq_dir, cds_chr_ids->at(i_chr));
		if (!check_file(bin_seq_fp))
		{
			sprintf(bin_seq_fp, "%s/%s.bin.gz", genome_bin_seq_dir, cds_chr_ids->at(i_chr));
			if (!check_file(bin_seq_fp))
			{
				fprintf(stderr, "Could not find the sequence file %s.\n", bin_seq_fp);
				exit(0);
			}
		}

		char* cur_chr_seq = load_binary_sequence_file(bin_seq_fp, l_chr_seq);
		fprintf(stderr, "Loaded %d nucleotide sequence.\n", l_chr_seq);

		vector<t_annot_region*>* cur_chr_cds_regs = get_regions_per_chromosome(cds_regs, cds_chr_ids->at(i_chr));
		for(int i_cds = 0; i_cds < cur_chr_cds_regs->size(); i_cds++)
		{
			char* cur_cds_reg_seq = new char[cur_chr_cds_regs->at(i_cds)->end - cur_chr_cds_regs->at(i_cds)->start + 5];
			int l_cds = cur_chr_cds_regs->at(i_cds)->end - cur_chr_cds_regs->at(i_cds)->start + 1;
			for(int i = 0; i < l_cds; i++)
			{
				cur_cds_reg_seq[i] = cur_chr_seq[cur_chr_cds_regs->at(i_cds)->start + i];
			} // i loop.
			cur_cds_reg_seq[l_cds] = 0;

			// Set the sequence of this cds.
			t_CDS_annotation_info* cur_cds_info = (t_CDS_annotation_info*)(cur_chr_cds_regs->at(i_cds)->annotation_info);
			cur_cds_info->sequence = cur_cds_reg_seq;
		} // i_cds loop.
		delete cur_chr_cds_regs;

		// Extract the sequences for introns.
		vector<t_annot_region*>* cur_chr_intron_regs = get_regions_per_chromosome(intronic_regs, cds_chr_ids->at(i_chr));
		for(int i_int = 0; i_int < cur_chr_intron_regs->size(); i_int++)
		{
			int l_intron = cur_chr_intron_regs->at(i_int)->end - cur_chr_intron_regs->at(i_int)->start + 1;
			char* cur_intron_reg_seq = new char[l_intron + 2];
			for(int i = 0; i < l_intron; i++)
			{
				cur_intron_reg_seq[i] = cur_chr_seq[cur_chr_intron_regs->at(i_int)->start + i];
			} // i loop.
			cur_intron_reg_seq[l_intron] = 0;

			t_intron_annotation_info* cur_int_info = (t_intron_annotation_info*)(cur_chr_intron_regs->at(i_int)->annotation_info);
			cur_int_info->sequence = cur_intron_reg_seq;
		} // i_int loop.

		delete [] cur_chr_seq;
	} // i_chr loop.

	// Go over all the transcripts and check sequence lengths.
	for(int i_tr = 0; i_tr < transcript_regs->size(); i_tr++)
	{
		vector<t_annot_region*>* transcript_cds_regs = ((t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info))->CDSs;
		sort(transcript_cds_regs->begin(), transcript_cds_regs->end(), sort_regions);
		t_string* cur_transcript_seq = new t_string();
		int l_transcript_so_far = 0;
		for(int i_cds = 0; i_cds < transcript_cds_regs->size(); i_cds++)
		{
			int l_cds = transcript_cds_regs->at(i_cds)->end - transcript_cds_regs->at(i_cds)->start + 1;

			t_CDS_annotation_info* cur_cds_info = (t_CDS_annotation_info*)(transcript_cds_regs->at(i_cds)->annotation_info);
			char* cur_cds_seq = cur_cds_info->sequence;
			cur_transcript_seq->concat_string(cur_cds_seq);

			// Set the first nucleotide's index in coding frame:
			cur_cds_info->coding_frame_i_per_first_nuc = l_transcript_so_far % 3;
			cur_cds_info->l_transcript_so_far = l_transcript_so_far;
			l_transcript_so_far += l_cds;
		} // i_cds loop.

		((t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info))->l_transcript = cur_transcript_seq->length();

		if(cur_transcript_seq->length() % 3 != 0)
		{
			/*
			fprintf(stderr, "%s @ %s:%d-%d, %d CDS's and %d nucleotides:\n", transcript_regs->at(i_tr)->name, 
												transcript_regs->at(i_tr)->chrom, transcript_regs->at(i_tr)->start, transcript_regs->at(i_tr)->end,
												transcript_cds_regs->size(), 
												cur_transcript_seq->length());

			fprintf(stderr, "%s\n", cur_transcript_seq->str());
			*/
			//getc(stdin);
			((t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info))->valid_3mer_frame = false;
		}
		else
		{
			//filtered_transcript_regs->push_back(transcript_regs->at(i_tr));
			((t_transcript_annotation_info*)(transcript_regs->at(i_tr)->annotation_info))->valid_3mer_frame = true;
		}
		delete cur_transcript_seq;
	} // i_tr loop.

	fprintf(stderr, "%d transcripts.\n", transcript_regs->size());
}

vector<t_annot_region*>* load_GFF_with_line_information(char* gff_fp)
{
	printf("Loading GFF file %s\n", gff_fp);

	vector<t_annot_region*>* gff_regions = new vector<t_annot_region*>();

	FILE* f_gff = open_f(gff_fp, "r");
	if(f_gff == NULL)
	{
		printf("Could not open GFF file @ %s\n", gff_fp);
		exit(0);
	}

	char* cur_line = getline(f_gff);
	while(cur_line != NULL)
	{
		int start;
		int end;
		char chrom[100];
		char strand;

		bool skip_line = check_line_skip(cur_line);

		if(!skip_line)
		{
			if(sscanf(cur_line, "%s %*s %*s %d %d %*s %c", chrom, &start, &end, &strand) != 4)
			{
				printf("Could not parse GFF file entry line:\n%s\n", cur_line);
				exit(0);
			}

			// Allocate and initialize the new region.
			t_annot_region* new_region = new t_annot_region();
			//new_region->start = start - GFF_START_BASE + CODEBASE_START_BASE;
			//new_region->end = end - GFF_END_BASE + CODEBASE_END_BASE;
			new_region->start = translate_coord(start, GFF_COORDS::start_base, CODEBASE_COORDS::start_base);
			new_region->end = translate_coord(end, GFF_COORDS::end_base, CODEBASE_COORDS::end_base);
			new_region->chrom = t_string::copy_me_str(chrom);
			normalize_chr_id(new_region->chrom);
			new_region->strand = strand;

			new_region->data = (void*)cur_line;

			// Add the new region to the list.
			gff_regions->push_back(new_region);
		}

		//delete [] cur_line;
		cur_line = getline(f_gff);
	} // narrowPeak file reading loop.

	fclose(f_gff);

	printf("Loaded %d regions from %s\n", gff_regions->size(), gff_fp);
	return(gff_regions);
}

/*
*/
char* get_gff_property_per_region(t_annot_region* entry, char* (*get_prop_per_group_str)(char*))
{
	t_gff_info* cur_reg_info = (t_gff_info*)(entry->data);
	char* res = (*get_prop_per_group_str)(cur_reg_info->group_str);

	return(res);
}

vector<t_annot_region*>* load_GFF(char* gff_fp)
{
	printf("Loading GFF file %s\n", gff_fp);

	vector<t_annot_region*>* gff_regions = new vector<t_annot_region*>();

	FILE* f_gff = NULL;
	if (t_string::compare_strings(gff_fp, "stdin"))
	{
		f_gff = stdin;
	}
	else
	{
		f_gff = open_f(gff_fp, "r");
	}

	char* cur_line = getline(f_gff);
	while(cur_line != NULL)
	{
		int start;
		int end;
		char chrom[100];
		char source[1000];
		char feature[1000];
		char strand;
		char score_str[1000];
		char* group_str;

		bool skip_line = check_line_skip(cur_line);

		if(!skip_line)
		{
			group_str = new char[strlen(cur_line) + 2];
			/*
0. seqname - The name of the sequence. Must be a chromosome or scaffold.
1. source - The program that generated this feature.
2. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
3. start - The starting position of the feature in the sequence. The first base is numbered 1.
4. end - The ending position of the feature (inclusive).
5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
6. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
7. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
8. group - All lines with the same group are linked together into a single item
			*/
			//Use t_string::get_first_n_tokens() to parse the loci and other information.
			int i_next_char;
			vector<t_string*>* line_tokens = t_string::get_first_n_tokens(cur_line, 9, "\t", i_next_char);
			//if(sscanf(cur_line, "%s %s %s %d %d %s %c %*s %s", chrom, source, feature, &start, &end, score_str, &strand, group_str) != 4)
			//{
			//	printf("Could not parse GFF file entry line:\n%s\n", cur_line);
			//	exit(0);
			//}	
			if(line_tokens->size() != 9)
			{
				fprintf(stderr, "Could not parse the line into 9 tokens: %s\n", cur_line);
				exit(0);
			}

			// Allocate and initialize the new region.
			t_annot_region* new_region = new t_annot_region();
			//new_region->start = start - GFF_START_BASE + CODEBASE_START_BASE;
			//new_region->end = end - GFF_END_BASE + CODEBASE_END_BASE;
			new_region->start = translate_coord(atoi(line_tokens->at(3)->str()), GFF_COORDS::start_base, CODEBASE_COORDS::start_base);
			new_region->end = translate_coord(atoi(line_tokens->at(4)->str()), GFF_COORDS::end_base, CODEBASE_COORDS::end_base);
			new_region->chrom = t_string::copy_me_str(line_tokens->at(0)->str());
			normalize_chr_id(new_region->chrom);
			new_region->strand = (line_tokens->at(6)->str())[0];
			new_region->score = (int)(atoi(line_tokens->at(5)->str()));
			new_region->name = NULL;

			// This is the pointer to the group string which contains the information.
			t_gff_info* cur_reg_info = new t_gff_info();
			cur_reg_info->feature_type = t_string::copy_me_str(line_tokens->at(2)->str());
			cur_reg_info->group_str = t_string::copy_me_str(line_tokens->at(8)->str());
			cur_reg_info->source_str = t_string::copy_me_str(line_tokens->at(1)->str());
			cur_reg_info->prop_ids = NULL;
			cur_reg_info->prop_vals = NULL;
			new_region->data = (void*)(cur_reg_info);
			new_region->intervals = NULL;	

			// Add the new region to the list.
			gff_regions->push_back(new_region);

			t_string::clean_tokens(line_tokens);
		}

		delete [] cur_line;
		cur_line = getline(f_gff);
	} // narrowPeak file reading loop.

	if (t_string::compare_strings(gff_fp, "stdin"))
	{
	}
	else
	{
		fclose(f_gff);
	}

	printf("Loaded %d regions from %s\n", gff_regions->size(), gff_fp);
	return(gff_regions);
}

vector<t_annot_region*>* get_gff_regions_per_source_string(vector<t_annot_region*>* all_regions,
															char* source_str)
{
	vector<t_annot_region*>* regions_of_source = new vector<t_annot_region*>();

	// Go over all the regions, extract the transcript id's for all the regions.	
	for(int i_reg = 0; i_reg < all_regions->size(); i_reg++)
	{
		// Extract the transcript id for the current transcript.
		t_gff_info* cur_reg_info = (t_gff_info*)(all_regions->at(i_reg)->data);
		char* cur_reg_src_str = cur_reg_info->source_str;

		// Is this a parent region?
		if(strcmp(cur_reg_src_str, source_str) == 0)
		{
			regions_of_source->push_back(all_regions->at(i_reg));
		}
	} // i_reg loop.

	return(regions_of_source);

}

void parse_group_strs(vector<t_annot_region*>* all_regions, void (*parse_grp_str)(t_annot_region*, vector<char*>*, vector<char*>*))
{
	for(int i_r = 0; i_r < all_regions->size(); i_r++)
	{
		vector<char*>* cur_reg_prop_ids = new vector<char*>();
		vector<char*>* cur_reg_prop_vals = new vector<char*>();

		// Sanity check.
		if(cur_reg_prop_ids->size() != cur_reg_prop_vals->size())
		{
			fprintf(stderr, "Group string parsing failed: %d, %d\n", cur_reg_prop_ids->size(), cur_reg_prop_vals->size());
			exit(0);
		}

		// Parse the group string for current entry.
		parse_grp_str(all_regions->at(i_r), cur_reg_prop_ids, cur_reg_prop_vals);
		
		t_gff_info* cur_gff_info = (t_gff_info*)(all_regions->at(i_r)->data);
		cur_gff_info->prop_ids = cur_reg_prop_ids;
		cur_gff_info->prop_vals = cur_reg_prop_vals;
	} // i_r loop.
}

vector<t_annot_region*>* get_gff_regions_per_feature_type(vector<t_annot_region*>* all_regions,
	char* feature_type)
{
	vector<t_annot_region*>* regions_of_type = new vector<t_annot_region*>();

	// Go over all the regions, extract the transcript id's for all the regions.	
	for(int i_reg = 0; i_reg < all_regions->size(); i_reg++)
	{
		// Extract the transcript id for the current transcript.
		t_gff_info* cur_reg_info = (t_gff_info*)(all_regions->at(i_reg)->data);
		char* cur_reg_feature_type = cur_reg_info->feature_type;

		// Is this a parent region?
		if(strcmp(cur_reg_feature_type, feature_type) == 0)
		{
			regions_of_type->push_back(all_regions->at(i_reg));
		}
	} // i_reg loop.

	return(regions_of_type);
}

char* get_prop_val_per_prop_id(t_annot_region* region, char* prop_id)
{
	if(region->data == NULL)
	{
		return(NULL);
	}

	vector<char*>* prop_ids = ((t_gff_info*)region->data)->prop_ids;
	vector<char*>* prop_vals = ((t_gff_info*)region->data)->prop_vals;

	for(int i_id = 0; i_id < prop_ids->size(); i_id++)
	{
		if(strcmp(prop_ids->at(i_id), prop_id) == 0)
		{
			return(prop_vals->at(i_id));
		}
	} // i_id loop.

	return(NULL);
}

// Form the list of regions with matching property value with id prop_id then get the regions containing each list.
void form_largest_regions_per_property_type(vector<t_annot_region*>* all_regions,
											vector<t_annot_region*>* merged_selected_regions,
											vector<char*>* sel_region_props, 
											char* prop_id)
{
	vector<vector<t_annot_region*>*>* sel_regions_lists = new vector<vector<t_annot_region*>*>();
	/*vector<char*>* sel_region_props = new vector<char*>();*/

	// Sort the regions with respect to their property values's, go over the regions and put the regions with the same proprty value together.
	fprintf(stderr, "Putting together the names of %d elements for merging features that belong to same element.\n", (int)all_regions->size());
	vector<char*>* all_region_props = new vector<char*>();
	for(int i_reg = 0; i_reg < all_regions->size(); i_reg++)
	{
		char* cur_reg_prop = get_prop_val_per_prop_id(all_regions->at(i_reg), prop_id);
		
		if(cur_reg_prop == NULL)
		{
			//fprintf(stderr, "Could not find property with id %s for region %d.\n", prop_id, i_reg);
			char* cur_noname_prop = new char[30];
			sprintf(cur_noname_prop, "NONAME_%d", i_reg);
			all_region_props->push_back(cur_noname_prop);
			all_regions->at(i_reg)->name = t_string::copy_me_str(cur_noname_prop);
		}
		else
		{
			all_region_props->push_back(cur_reg_prop);
			all_regions->at(i_reg)->name = t_string::copy_me_str(cur_reg_prop);
		}
	} // i_reg loop.
	fprintf(stderr, "Identified %d element names.\n", all_region_props->size());

	fprintf(stderr, "Generating the sorting index.\n");
	vector<int>* string_sort_idx = t_string::get_string_sorting_idx(all_region_props);
	//delete(all_region_props);

	char* cur_list_prop = all_region_props->at(string_sort_idx->at(0));
	sel_region_props->push_back(t_string::copy_me_str(cur_list_prop));
	sel_regions_lists->push_back(new vector<t_annot_region*>());
	sel_regions_lists->back()->push_back(all_regions->at(string_sort_idx->at(0)));
	for(int i_reg = 1; i_reg < string_sort_idx->size(); i_reg++)
	{
		if(i_reg % 10000 == 0)
		{
			fprintf(stderr, "%d. region.        \r", i_reg);
		}

		int cur_sorted_i = string_sort_idx->at(i_reg);

		// Does this entry have the same property as the previous entry?
		if(cur_list_prop != NULL &&
			strcmp(cur_list_prop, all_region_props->at(cur_sorted_i)) == 0)
		{
			//cur_sel_region_list->push_back(all_regions->at(cur_sorted_i));
			sel_regions_lists->back()->push_back(all_regions->at(cur_sorted_i));
		}
		else
		{
			// Copy the new entry.
			cur_list_prop = all_region_props->at(cur_sorted_i);

			sel_region_props->push_back(t_string::copy_me_str(cur_list_prop));

			// Allocate the new stuff, again.
			sel_regions_lists->push_back(new vector<t_annot_region*>());
			sel_regions_lists->back()->push_back(all_regions->at(cur_sorted_i));
		}
	} // i_reg loop.

	// Go over all the selected region lists, get the largest sets. Assign the sub regions to the intervals.
	//vector<t_annot_region*>* merged_selected_regions = new vector<t_annot_region*>();
	for(int i_reg = 0; i_reg < sel_regions_lists->size(); i_reg++)
	{
		if(i_reg % 1000 == 0)
		{
			fprintf(stderr, "Forming parent regions %d(%d)       \r", i_reg, sel_regions_lists->size());
		}

		// Sort the sub-regions for the current selected region.
		sort(sel_regions_lists->at(i_reg)->begin(), sel_regions_lists->at(i_reg)->end(), sort_regions);

		t_annot_region* new_merged_region = new t_annot_region();
		new_merged_region->chrom = t_string::copy_me_str(sel_regions_lists->at(i_reg)->at(0)->chrom);
		new_merged_region->start = sel_regions_lists->at(i_reg)->at(0)->start;
		new_merged_region->end = sel_regions_lists->at(i_reg)->back()->end;
		new_merged_region->strand = sel_regions_lists->at(i_reg)->at(0)->strand;
		new_merged_region->data = NULL;
		new_merged_region->name = t_string::copy_me_str(sel_regions_lists->at(i_reg)->at(0)->name);

		// Setup the intervals.
		// Sanity check: Do the selected regions have the same property value?
		for(int i_sg = 0; i_sg < sel_regions_lists->at(i_reg)->size(); i_sg++)
		{
			char* cur_region_prop = get_prop_val_per_prop_id(sel_regions_lists->at(i_reg)->at(i_sg), prop_id);
			if(cur_region_prop == NULL)
			{
				//if(strcmp(sel_region_props->at(i_reg), "NONAME") != 0)
				//{

				//}
			}
			else if(strcmp(sel_region_props->at(i_reg), cur_region_prop) != 0)
			{
				fprintf(stderr, "The property values are not matching: %s, %s\n", sel_region_props->at(i_sg), cur_region_prop);
				getc(stdin);
			}
		} // i_sg loop.

		new_merged_region->intervals = sel_regions_lists->at(i_reg);

		// Set the gff information for the new parent region.
		t_gff_info* new_gff_info = new t_gff_info();
		new_gff_info->feature_type = NULL;
		new_gff_info->group_str = NULL;
		new_gff_info->source_str = NULL;
		new_gff_info->prop_ids = new vector<char*>();
		new_gff_info->prop_ids->push_back(t_string::copy_me_str(prop_id));
		new_gff_info->prop_vals = new vector<char*>();
		new_gff_info->prop_vals->push_back(t_string::copy_me_str(sel_region_props->at(i_reg)));

		fprintf(stderr, "Setup the parent region %s: %d\n", sel_region_props->at(i_reg), new_merged_region->intervals->size());

		new_merged_region->data = new_gff_info;

		// Add the new merged region.
		merged_selected_regions->push_back(new_merged_region);
	}  // i_reg loop.

	//return(merged_selected_regions);
}
//
//vector<char*>* get_gff_string_properties_per_entry(vector<t_annot_region*>* all_regions,
//												char* (*get_property_string)(t_annot_region*))
//{
//	vector<char*>* props = new vector<char*>();
//
//	for(int i_t = 0; i_t < all_regions->size(); i_t++)
//	{
//		char* cur_prop_str = get_property_string(all_regions->at(i_t));
//
//		// The id is not found.
//		if(cur_prop_str== NULL)
//		{
//			fprintf(stderr, "Could not extract the property.\n");
//		}
//		else
//		{
//			props->push_back(cur_prop_str);
//		}
//	} // i_t loop.
//
//	return(props);
//}

// For a set of mixed regions, associate the parent regions with the children.
// For example, associate the exons with the transcripts,
// Associate the transcripts with the genes,
// Associate the start_codons with the transcripts, end_codons with transcripts,
// Associate custom annotations with each other.
void associate_gff_children_2_parents(vector<t_annot_region*>* parent_regions,
	vector<t_annot_region*>* children_regions,
	vector<char*>* (*get_parent_feature_id_per_child_group_str)(t_annot_region*), // Use this to extract the id for the parent using the group string for child.
	char* (*get_parent_feature_id_per_parent_group_str)(t_annot_region*)) // Use this to extract the id for the parent using the group string for parent.
{
	// Extract the parent id's and id's for all the entries and buffer them.
	vector<char*>* parent_feature_ids = new vector<char*>();
	for(int i_p = 0; i_p < parent_regions->size(); i_p++)
	{
		//t_gff_info* cur_parent_reg_info = (t_gff_info*)(parent_regions->at(i_p)->data);
		//char* cur_parent_reg_grp_str = cur_parent_reg_info->group_str;
		char* cur_parent_reg_feature_id = (*get_parent_feature_id_per_parent_group_str)(parent_regions->at(i_p));
		parent_feature_ids->push_back(t_string::copy_me_str(cur_parent_reg_feature_id));

		parent_regions->at(i_p)->intervals = new vector<t_annot_region*>();
	} // i_p loop.

	vector<vector<char*>*>* parent_feature_ids_per_children = new vector<vector<char*>*>();
	for(int i_c = 0; i_c < children_regions->size(); i_c++)
	{
		//t_gff_info* cur_child_reg_info = (t_gff_info*)(children_regions->at(i_c)->data);
		//char* cur_child_reg_grp_str = cur_child_reg_info->group_str;
		vector<char*>* cur_parent_reg_feature_ids = (*get_parent_feature_id_per_child_group_str)(children_regions->at(i_c));
		parent_feature_ids_per_children->push_back(cur_parent_reg_feature_ids);

		children_regions->at(i_c)->intervals = new vector<t_annot_region*>();
	} // i_c loop.

	// Go over the all the regions, associate the children regions with the parents.
	for(int i_c = 0; i_c < children_regions->size(); i_c++)
	{
		if(i_c % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. child. (%d)            \r", i_c, children_regions->size());
		}

		vector<char*>* cur_child_parent_ids = parent_feature_ids_per_children->at(i_c);

		// Go over all the parents of the current child and assign the child to the parents.
		for(int i_c_p = 0; i_c_p < cur_child_parent_ids->size(); i_c_p++)
		{
			bool found_parent = false;
			for(int i_p = 0;
				!found_parent && i_p < parent_regions->size();
				i_p++)
			{
				// Compare this parent id of the current child with the id of the current parent.
				if(strcmp(cur_child_parent_ids->at(i_c_p), parent_feature_ids->at(i_p)) == 0)
				{
					found_parent = true;
					parent_regions->at(i_p)->intervals->push_back(children_regions->at(i_c));
				}
			} // i_p loop.

			// Is the parent found?
			if(!found_parent)
			{
				fprintf(stderr, "Could not find the parent for child region: %s:%d-%d with parent region id %s\n", 
					children_regions->at(i_c)->chrom, children_regions->at(i_c)->start, children_regions->at(i_c)->end, 
					cur_child_parent_ids->at(i_c_p));

				getc(stdin);
			}
		} // i_c_p loop.
	} // i_reg loop.

	// Go over all the parents and determine the parents that do not have any children.
	fprintf(stderr, "Checking for parents with no children.\n");
	for(int i_p = 0;
		i_p < parent_regions->size();
		i_p++)
	{
		// Compare this parent id of the current child with the id of the current parent.
		if(parent_regions->at(i_p)->intervals->size() == 0)
		{
			fprintf(stderr, "There are no children for %s:%d-%d\n", parent_regions->at(i_p)->chrom, parent_regions->at(i_p)->start, parent_regions->at(i_p)->end);
			getc(stdin);
		}
	} // i_p loop.

	// All of the parents for all of the children are processed.
}

t_annot_region* copy_gff_entry(t_annot_region* gff_entry)
{
	t_annot_region* dup_entry = new t_annot_region();
	dup_entry->chrom = t_string::copy_me_str(gff_entry->chrom);
	dup_entry->start = gff_entry->start;
	dup_entry->end = gff_entry->end;
	dup_entry->strand = gff_entry->strand;
	dup_entry->intervals = NULL;

	// Copy gff_info.
	t_gff_info* orig_gff_info = (t_gff_info*)(gff_entry->data);
	t_gff_info* dup_gff_info = new t_gff_info();
	dup_gff_info->feature_type = t_string::copy_me_str(orig_gff_info->feature_type);
	dup_gff_info->group_str = t_string::copy_me_str(orig_gff_info->group_str);
	dup_gff_info->source_str = t_string::copy_me_str(orig_gff_info->source_str);
	dup_gff_info->prop_ids = new vector<char*>();
	dup_gff_info->prop_vals = new vector<char*>();
	
	for(int i_prop = 0; i_prop < orig_gff_info->prop_ids->size(); i_prop++)
	{
		dup_gff_info->prop_ids->push_back(t_string::copy_me_str(orig_gff_info->prop_ids->at(i_prop)));
		dup_gff_info->prop_vals->push_back(t_string::copy_me_str(orig_gff_info->prop_vals->at(i_prop)));
	} // i_prop loop.

	dup_entry->data = dup_gff_info;

	return(dup_entry);
}