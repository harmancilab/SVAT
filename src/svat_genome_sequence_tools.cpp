#include <stdio.h>
#include <stdlib.h>
#include "svat_genome_sequence_tools.h"
#include "file_utils.h"
#include "svat_ansi_string.h"
#include "svat_genomics_coords.h"
//#include "svat_alignment_tools.h"
#include "svat_annot_region_tools.h"
#include <math.h>
#include "svat_rng.h"
#include "svat_nomenclature.h"
#include "svat_seed_manager.h"
#include "svat_nucleotide.h"
#include <string.h>
#include <vector>
#include <algorithm>
#include <ctype.h>

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))


using namespace std;

bool _DUMP_GENOME_SEQUENCE_MESSAGES_ = false;

int map_nuc_to_sig(char nuc)
{
	switch(nuc)
	{
		case 'A':
			return(0);
			break;
		case 'C':
			return(1);
			break;
		case 'G':
			return(2);
			break;
		case 'U':
			return(0);
			break;
		case 'T':
			return(4);
			break;
		default:
			return(0);
			break;
	}
}

void get_kmer_signal_per_numerized_sequence_signal(char* seq_signal, int l_signal, double* kmer_signal, int k)
{
	// Make sure that k is odd.
	if(k % 2 != 1)
	{
		fprintf(stderr, "Please enter an odd k value.\n");
		exit(0);
	}

	int l_half_k = (k - 1) / 2;

	for(int i_nuc = 1; i_nuc <= l_signal; i_nuc++)
	{
		// Read the kmer around this nucleotide.
		int cur_val = 0;
		int cur_base = 1;

		//fprintf(stderr, "Current k-mer:\n");
		for(int j_nuc = MAX(i_nuc-l_half_k, 0); j_nuc <= MIN(i_nuc+l_half_k, l_signal-1); j_nuc++)
		{
			//fprintf(stderr, "%d ", numerized_seq_signal[j_nuc]);
			cur_val += nuc_2_num(seq_signal[j_nuc]) * cur_base;
			cur_base *= 4;
		} // j_nuc loop.
		//fprintf(stderr, "Value: %d\n", cur_val);
		//getc(stdin);

		kmer_signal[i_nuc] = ((double)cur_val);
	} // i_nuc loop.
}

void get_numerized_sequence_signals(char* sequence, int l_signal,
	double* walk_signal,
	double* nucleotide_signal,
	double* gc_content_signal)
{
	fprintf(stderr, "Generating the sequence signal for %d nucleotides.\n", l_signal);

	double gc_count = 0;
	int cur_cumul_walk_sig = 0;
	double cur_gc_content = 0;
	double l_sequenced = 0;
	for(int i_nuc = 0; i_nuc < l_signal; i_nuc++)
	{
		cur_gc_content = 0;
		if(toupper(sequence[i_nuc]) != 'N')
		{
			l_sequenced++;

			if(toupper(sequence[i_nuc]) == 'C' || 
						toupper(sequence[i_nuc]) == 'G')
			{
				cur_gc_content = 1;
				gc_count++;
			}
		}

		if(toupper(sequence[i_nuc]) == 'A' || toupper(sequence[i_nuc]) == 'G' || 
		toupper(sequence[i_nuc]) == 'C' || toupper(sequence[i_nuc]) == 'T' ||
		toupper(sequence[i_nuc]) == 'U')
		{
			if(toupper(sequence[i_nuc]) == 'A' || toupper(sequence[i_nuc]) == 'G')
			{
				cur_cumul_walk_sig--;
			}
			else
			{
				cur_cumul_walk_sig++;
			}
		}

		gc_content_signal[i_nuc+1] = cur_gc_content;
		walk_signal[i_nuc+1] = cur_cumul_walk_sig;
		nucleotide_signal[i_nuc+1] = nuc_2_num(toupper(sequence[i_nuc]));
	} // i_nuc loop.

	fprintf(stderr, "GC percentage: %lf\n", gc_count / l_sequenced);
}

// Given the 0-based seq_i indexing, returns the genome index with the assumed index convention.
bool get_genome_i_per_seq_i(char strand, 
	int genome_start, 
	int genome_end, 
	int& genome_i, int seq_i)
{
	int l_seq = genome_end - genome_start + 1;
	genome_i = seq_i + genome_start;

	// Revert the index if the mirna is on negative strand.
	if(strand == '-')
	{
		//seq_i = l_seq - seq_i - 1;
		genome_i = l_seq + genome_start - seq_i - 1;
	}

	return(true);
}

// This returns 0 based index for the sequence.
bool get_seq_i_per_genome_i(char strand, 
	int genome_start, 
	int genome_end, 
	int genome_i, int& seq_i)
{
	// genome_i must be within the element of interest.
	if(genome_i < genome_start ||
		genome_i > genome_end)
	{
		return(false);
	}

	// genome_i must be within the element of interest: If it is not, fix it.
	if(genome_i < genome_start)
	{
		genome_i = genome_start;
	}
	
	if(genome_i > genome_end)
	{
		genome_i = genome_end;
	}

	int l_seq = genome_end - genome_start + 1;
	seq_i = genome_i - genome_start;

	// Revert the index if the mirna is on negative strand.
	if(strand == '-')
	{
		seq_i = l_seq - seq_i - 1;
	}

	return(true);
}

//void binarize_buffered_fasta_file(char* fasta_fp, char* bin_dir)
void binarize_fasta_file(char* fasta_fp, char* bin_dir)
{
	printf("Binarizing %s.\n", fasta_fp);
	FILE* f_fasta = open_f(fasta_fp, "r");

	char* cur_line = NULL;
	char* cur_entry_buffer = new char[250 * 1000 * 1000];
	int cur_entry_i = 1; // This ensures that we are 1 based.
	char cur_entry_id[1000];

	char chr_ids_fp[1000];
	sprintf(chr_ids_fp, "%s/chr_ids.txt", bin_dir);
	FILE* f_chr_ids = open_f(chr_ids_fp, "a");

	while(1)
	{
		// Process the current buffer.
		cur_line = getline(f_fasta);

		if(cur_line == NULL)
		{
			// File ended, dump the last entry if there are values in it.
			if(cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin.gz", bin_dir, cur_entry_id);
				fprintf(f_chr_ids, "%s\n", cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				close_f(f_bin, cur_entry_bin_fp);
			}
			break;
		}
		else if(cur_line[0] == '>')
		{
			if(cur_entry_i > 1)
			{
				char cur_entry_bin_fp[1000];
				normalize_chr_id(cur_entry_id);
				sprintf(cur_entry_bin_fp, "%s/%s.bin.gz", bin_dir, cur_entry_id);
				fprintf(f_chr_ids, "%s\n", cur_entry_id);
				FILE* f_bin = open_f(cur_entry_bin_fp, "wb");

				// Open the new binary file.
				fprintf(stderr, "Dumping %s (%d)\n", cur_entry_id, cur_entry_i);

				// Dump the sequence buffer.
				fwrite(&cur_entry_i, sizeof(int), 1, f_bin);
				fwrite(cur_entry_buffer, sizeof(char), cur_entry_i, f_bin);

				// Close current file.
				close_f(f_bin, cur_entry_bin_fp);
			}

			// Update the id, reset the counter.
			strcpy(cur_entry_id, &(cur_line[1]));
			cur_entry_i = 1;
		}
		else
		{
			// Concatenate the current sequence line.
			int l_cur_line = t_string::string_length(cur_line);
			for(int i = 0; i < l_cur_line; i++)
			{
				cur_entry_buffer[cur_entry_i] = cur_line[i];
				cur_entry_i++;
			} // i loop.
		}

		delete [] cur_line;
	} // file reading loop.

	close_f(f_fasta, fasta_fp);
	fclose(f_chr_ids);

	delete [] cur_entry_buffer;
}

//char* numerize_sequence_signal(char* seq_signal, int l_seq)
//{
//	char* numerized_seq_buffer = new char[l_seq+1];
//
//	for(int i_nuc = 0; i_nuc < l_seq; i_nuc++)
//	{
//		numerized_seq_buffer[i_nuc] = nuc_2_num(seq_signal[i_nuc]);
//	} // i_nuc loop.
//
//	return(numerized_seq_buffer);
//}

char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq)
{
	FILE* f_bin_seq = open_f(bin_seq_fp, "rb");
	int l_bin_seq = 0;

	// Read the length.
	fread(&l_bin_seq, sizeof(int), 1, f_bin_seq);

	l_seq = l_bin_seq;

	// Read the sequence.
	char* bin_seq = new char[l_bin_seq+2];
	fread(bin_seq, sizeof(char), l_bin_seq+1, f_bin_seq);
	fclose(f_bin_seq);

	return(bin_seq);
}

//void binarize_fasta_file(char* fasta_fp, char* bin_dir)
//{
//	printf("Binarizing %s.\n", fasta_fp);
//	FILE* f_fasta = open_f(fasta_fp, "r");
//
//	// Read the first line, which is the identifier line.
//	char bin_fp[1000];
//	FILE* f_bin = NULL;
//
//	char* cur_line = NULL;
//
//	while(1)
//	{
//		cur_line = getline(f_fasta);
//
//		if(cur_line == NULL)
//		{
//			break;
//		}
//
//		if(cur_line[0] == '>')
//		{
//			// Close the current binary file.
//			if(f_bin != NULL)
//			{
//				fclose(f_bin);
//			}
//
//			// Open the new binary file.
//			fprintf(stderr, "Loading %s\n", &(cur_line[1]));
//			char chrom[1000];
//			strcpy(chrom, &(cur_line[1]));
//			normalize_chr_id(chrom);
//			fprintf(stderr, "Normalized chromosome name: %s\n", chrom);
//			sprintf(bin_fp, "%s/%s.bin", bin_dir, chrom);
//			f_bin = open_f(bin_fp, "wb");
//		}
//		else
//		{
//			// Concatenate the sequence information.
//			if(f_bin != NULL)
//			{
//				fwrite(cur_line, strlen(cur_line), 1, f_bin);
//			}
//			else
//			{
//				printf("The file pointer is currently null while dumping:\n%s\n", cur_line);
//				exit(0);
//			}
//		}
//	} // file reading loop.
//
//	fclose(f_fasta);
//	fclose(f_bin);
//}

void binarize_genome_per_plain(char* fp, char* op_dir)
{
}

vector<t_annot_region*>* load_BED_with_ancestral_derived_seqs(char* bed_w_seq_fp)
{
	vector<t_annot_region*>* regs_w_seqs = new vector<t_annot_region*>();

	FILE* f_bed_w_seq = open_f(bed_w_seq_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_bed_w_seq);
		if(cur_line == NULL)
		{
			break;
		}

		char chrom[100];
		int start;
		int end;
		char strand;
		char* name = new char[t_string::string_length(cur_line) + 2];
		char** cur_anc_der_seqs = new char*[2];
		char* cur_ancestral_seq = new char[t_string::string_length(cur_line) + 2];
		char* cur_derived_seq = new char[t_string::string_length(cur_line) + 2];
		cur_anc_der_seqs[0] = cur_ancestral_seq;
		cur_anc_der_seqs[1] = cur_derived_seq;

		if(sscanf(cur_line, "%s %d %d %s %*s %c %s %s", chrom, &start, &end, name, &strand, 
			cur_ancestral_seq, cur_derived_seq) != 7)
		{
			fprintf(stderr, "Could not parse the sequence bed file line:\n%s\n", cur_line);
			exit(0);
		}

		t_annot_region* cur_region = get_empty_region();
		cur_region->chrom = t_string::copy_me_str(chrom);
		cur_region->start = start;
		cur_region->end = end;
		cur_region->name = t_string::copy_me_str(name);
		cur_region->strand = strand;

		// The sequence is stored in the data entry.
		cur_region->data = (void*)cur_anc_der_seqs;

		regs_w_seqs->push_back(cur_region);

		delete [] name;
	} // file reading loop.

	fclose(f_bed_w_seq);

	return(regs_w_seqs);
}

// Load the bed file with sequence information that is dumped by sequence extraction.
vector<t_annot_region*>* load_BED_with_sequences(char* bed_w_seq_fp)
{
	vector<t_annot_region*>* regs_w_seqs = new vector<t_annot_region*>();

	FILE* f_bed_w_seq = open_f(bed_w_seq_fp, "r");

	while(1)
	{
		char* cur_line = getline(f_bed_w_seq);
		if(cur_line == NULL)
		{
			break;
		}

		char chrom[100];
		int start;
		int end;
		char strand;
		char name[100];
		char* cur_seq = new char[t_string::string_length(cur_line) + 2];		

		if(sscanf(cur_line, "%s %d %d %s %*s %c %s", chrom, 
			&start, &end, name, &strand, cur_seq) != 6)
		{
			fprintf(stderr, "Could not parse the sequence bed file line:\n%s\n", cur_line);
			exit(0);
		}

		t_annot_region* cur_region = get_empty_region();
		cur_region->chrom = t_string::copy_me_str(chrom);
		cur_region->start = translate_coord(start, BED_COORDS::start_base, CODEBASE_COORDS::start_base);
		cur_region->end = translate_coord(end, BED_COORDS::end_base, CODEBASE_COORDS::end_base);
		cur_region->strand = strand;

		// The sequence is stored in the data entry.
		cur_region->data = (void*)cur_seq;

		regs_w_seqs->push_back(cur_region);
	} // file reading loop.

	fclose(f_bed_w_seq);

	return(regs_w_seqs);
}

void batch_extract_interval_sequences(char* chr_data_dir, vector<t_annot_region*>* regions)
{
	t_restr_annot_region_list* retstr_regions = restructure_annot_regions(regions);

	for(int i_chr = 0; i_chr < (int)retstr_regions->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing %s\n", retstr_regions->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_regions = retstr_regions->regions_per_chrom[i_chr];
		for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
		{
			// Process the current transcript: Merge the exons
			vector<t_annot_region*>* merged_exons = merge_annot_regions(cur_chr_regions->at(i_reg)->intervals, 0, true);

			// Sort the merged regions.
			sort(merged_exons->begin(), merged_exons->end(), sort_regions);

			t_string* cur_transcript_seq = new t_string();

			// Extract the sequence for all the exons.
			for(int i_ex = 0; i_ex < (int)merged_exons->size(); i_ex++)
			{
				if(merged_exons->at(i_ex)->strand == '-')
				{
					merged_exons->at(i_ex)->strand = '+';
				}

				char* cur_exon_seq = extract_sequence(chr_data_dir, retstr_regions->chr_ids->at(i_chr), merged_exons->at(i_ex)->start, merged_exons->at(i_ex)->end, merged_exons->at(i_ex)->strand);
				cur_transcript_seq->concat_string(cur_exon_seq);
				delete [] cur_exon_seq;
			} // i_ex loop.

			// If the transcript is on negative strand, revert the sequence.
			char* transcript_seq = t_string::copy_me_str(cur_transcript_seq->str());
			if(cur_chr_regions->at(i_reg)->strand == '-')
			{
				//cur_transcript_seq->revert();
				delete [] transcript_seq;
				transcript_seq = get_reverse_complement(cur_transcript_seq->str(), 0, t_string::string_length(cur_transcript_seq->str())-1);
			}

			// Copy the data.
			cur_chr_regions->at(i_reg)->data = transcript_seq;

			delete(cur_transcript_seq);
			// Memory for merged exons must be freeed.
		} // i_reg loop.
	} // i_chr loop.
}

/*
Extract the sequences for multiple regions at one time.
Start and end are 1 based.
*/
//vector<t_annot_region*>* batch_extract_region_sequences(char* chr_data_dir, char* regions_bed_fp)
void batch_extract_region_sequences(char* chr_data_dir, vector<t_annot_region*>* regions)
{
	//vector<t_annot_region*>* regions = load_BED(regions_bed_fp);

	vector<char*>* chr_ids = get_chr_ids(regions);

	for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_regions = get_regions_per_chromosome(regions, chr_ids->at(i_chr));

		fprintf(stderr, "Extracting sequences in %s\n", chr_ids->at(i_chr));

		// Sort the regions.
		sort(cur_chr_regions->begin(), cur_chr_regions->end(), sort_regions);

		char chr_seq_fp[1000];
		sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, chr_ids->at(i_chr));
		if(check_file(chr_seq_fp))
		{
			FILE* f_chr = open_f(chr_seq_fp, "rb");

			// Go over all the regions.
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				fprintf(stderr, "%d (%d)             \r", i_reg, cur_chr_regions->size());

				int l_reg = cur_chr_regions->at(i_reg)->end - cur_chr_regions->at(i_reg)->start + 1;
				int f_start = cur_chr_regions->at(i_reg)->start - 1 + sizeof(int);
				int f_end = f_start + l_reg - 1;

				int length = f_end - f_start + 1;

				char* seq_buffer = new char[length + 2];
				memset(seq_buffer, 0, sizeof(char) * (length + 1));

				// Read the data.
				fseek(f_chr, f_start, SEEK_SET);

				// Fill the buffer.
				fread(seq_buffer, length, 1, f_chr);

				// If the requested strand is on negative strand, do reverse complementing.
				if(cur_chr_regions->at(i_reg)->strand == '-')
				{
					// Replace the sequence buffer with its reverse complement.
					char* reverse_complement_seq = get_reverse_complement(seq_buffer, 0, t_string::string_length(seq_buffer) - 1);
					delete [] seq_buffer;
					seq_buffer = reverse_complement_seq;
				}

				//fprintf(stderr, "%s\t%d\t%d\t.\t.\t%c\t%s\n", cur_chr_regions->at(i_reg)->chrom, 
				//	cur_chr_regions->at(i_reg)->start, cur_chr_regions->at(i_reg)->end, 
				//	cur_chr_regions->at(i_reg)->strand, seq_buffer);
			
				// Save the sequence to data entry of region.
				cur_chr_regions->at(i_reg)->data = (void*)seq_buffer;
			} // i_reg.

			fprintf(stderr, "\n");

			fclose(f_chr);
		}
		else
		{
			fprintf(stderr, "Could not find %s\n", chr_seq_fp);
			for(int i_reg = 0; i_reg < (int)cur_chr_regions->size(); i_reg++)
			{
				cur_chr_regions->at(i_reg)->data = NULL;
			} // i_reg loop.
		}
	} // i_chr loop.
}

char* load_whole_chromosome_sequence_per_bin_file(char* chr_fp, int& n_nucs)
{
	FILE* f_chr = open_f(chr_fp, "rb");
	
	// Get the file size.
	fseek (f_chr, 0 , SEEK_END);
	n_nucs = (int)(ftell(f_chr));
	rewind(f_chr);

	char* nucs = new char[n_nucs+2]; 

	fread(nucs, sizeof(char), n_nucs, f_chr);

	fclose(f_chr);

	return(nucs);
}

/*
Start and end are 1 based.
*/
char* extract_sequence(char* chr_data_dir, char* raw_chr_name, int start, int end, char strand, bool upper_case)
{
	char chr_seq_fp[1000];
	sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, raw_chr_name);
	FILE* f_chr = fopen(chr_seq_fp, "rb");

	if(f_chr == NULL)
	{
		// Try to fix the chromosome file path: Add 'chr' to beginning.
		sprintf(chr_seq_fp, "%s/%s.bin", chr_data_dir, raw_chr_name);
		f_chr = fopen(chr_seq_fp, "rb");

		if(f_chr == NULL)
		{
			fprintf(stderr, "Could not identify the conforming chromosome name from id %s, returning NULL.\n", raw_chr_name);
			return(NULL);
		}
	}

	// Read the binary sequence length.
	int l_bin_seq = 0;
	fread(&l_bin_seq, sizeof(int), 1, f_chr);

	int f_start = start + sizeof(int);
	int f_end = f_start + (end-start);

	char* seq_buffer = NULL;

	if(f_start < f_end)
	{
		int length = (end - start + 1);

		seq_buffer = new char[length + 2];
		memset(seq_buffer, 0, sizeof(char) * (length + 1));

		// Go to the starting position: Note that the binary file is 1 based; i.e., there is an extra byte at the beginning
		// of the file that corresponds to 0th character, which is not meaningful in codebase coordinates. Therefore
		// following puts us to the right position in the file.
		fseek(f_chr, f_start, SEEK_SET);

		// Fill the buffer.
		fread(seq_buffer, sizeof(char), length, f_chr);

		seq_buffer[length] = 0;
	} // start, end coparison.
	else
	{
		// Get the file size.
		fseek (f_chr, 0 , SEEK_END);
		int length = (int)(ftell(f_chr));
		
		seq_buffer = new char[length+2]; 
		memset(seq_buffer, 0, sizeof(char) * (length + 1));

		rewind(f_chr);

		fread(seq_buffer, sizeof(char), length, f_chr);

		seq_buffer[length] = 0;
	} // whole chromosome reading check.

	// If the requested strand is on negative strand, do reverse complementing.
	if(strand == '-')
	{
		// Replace the sequence buffer with its reverse complement.
		char* reverse_complement_seq = get_reverse_complement(seq_buffer, 0, t_string::string_length(seq_buffer) - 1);
		delete [] seq_buffer;
		seq_buffer = reverse_complement_seq;
	}

	fclose(f_chr);

	if(upper_case)
	{
		t_string::to_upper(seq_buffer);
	}

	return(seq_buffer);
}

/*
For reverse strand, the transcribed RNA is exactly same as forward strand with the exception that T(t)'s
are converted to U(u)'s.
*/
char* get_reverse_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* reverse_RNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(reverse_RNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	int i_rna_nuc = 0;
	for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	{
		if(cur_aln_line[i_nuc] == 'T')
		{
			reverse_RNA[i_rna_nuc] = 'U';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			reverse_RNA[i_rna_nuc] = 'u';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'A' || cur_aln_line[i_nuc] == 'a' ||
			cur_aln_line[i_nuc] == 'C' || cur_aln_line[i_nuc] == 'c' ||
			cur_aln_line[i_nuc] == 'G' || cur_aln_line[i_nuc] == 'g' ||
			cur_aln_line[i_nuc] == '-')
		{
			reverse_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
		else
		{
			printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			reverse_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
	} // i_nuc loop.

	return(reverse_RNA);
}

void get_reverse_complement(char* cur_aln_line, char* rev_comp_DNA, int i_nuc_min, int i_nuc_max)
{
	int i_dna_nuc = 0;

	// The forward strand RNA is transcribed from 3' end.
	for(int i_nuc = i_nuc_max; i_nuc >= i_nuc_min; i_nuc--)
	{
		if(cur_aln_line[i_nuc] == 'A')
		{
			rev_comp_DNA[i_dna_nuc] = 'T';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'T')
		{
			rev_comp_DNA[i_dna_nuc] = 'A';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'G')
		{
			rev_comp_DNA[i_dna_nuc] = 'C';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'C')
		{
			rev_comp_DNA[i_dna_nuc] = 'G';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'a')
		{
			rev_comp_DNA[i_dna_nuc] = 't';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			rev_comp_DNA[i_dna_nuc] = 'a';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'g')
		{
			rev_comp_DNA[i_dna_nuc] = 'c';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'c')
		{
			rev_comp_DNA[i_dna_nuc] = 'g';
			i_dna_nuc++;
		}
		//else if(cur_aln_line[i_nuc] == '-')
		//{
		//	rev_comp_DNA[i_dna_nuc] = '-';
		//	i_dna_nuc++;
		//}
		else
		{
			//printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			rev_comp_DNA[i_dna_nuc] = cur_aln_line[i_nuc];
			i_dna_nuc++;
		}
	} // i_nuc loop.
}

char* get_forward_complement(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* fore_comp_DNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(fore_comp_DNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	int i_dna_nuc = 0;

	// Go straight.
	for (int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	{
		if (cur_aln_line[i_nuc] == 'A')
		{
			fore_comp_DNA[i_dna_nuc] = 'T';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'T')
		{
			fore_comp_DNA[i_dna_nuc] = 'A';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'G')
		{
			fore_comp_DNA[i_dna_nuc] = 'C';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'C')
		{
			fore_comp_DNA[i_dna_nuc] = 'G';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'a')
		{
			fore_comp_DNA[i_dna_nuc] = 't';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 't')
		{
			fore_comp_DNA[i_dna_nuc] = 'a';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'g')
		{
			fore_comp_DNA[i_dna_nuc] = 'c';
			i_dna_nuc++;
		}
		else if (cur_aln_line[i_nuc] == 'c')
		{
			fore_comp_DNA[i_dna_nuc] = 'g';
			i_dna_nuc++;
		}
		//else if(cur_aln_line[i_nuc] == '-')
		//{
		//	rev_comp_DNA[i_dna_nuc] = '-';
		//	i_dna_nuc++;
		//}
		else
		{
			//printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			fore_comp_DNA[i_dna_nuc] = cur_aln_line[i_nuc];
			i_dna_nuc++;
		}
	} // i_nuc loop.

	  //printf("%s", rev_comp_DNA);
	  //getc(stdin);

	return(fore_comp_DNA);
}

char* get_reverse_complement(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* rev_comp_DNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(rev_comp_DNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	//printf("Forward strand RNA for strand:\n");
	//for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	//{
	//	printf("%c", cur_aln_line[i_nuc]);
	//}
	//printf("\n");

	int i_dna_nuc = 0;

	// The forward strand RNA is transcribed from 3' end.
	for(int i_nuc = i_nuc_max; i_nuc >= i_nuc_min; i_nuc--)
	{
		if(cur_aln_line[i_nuc] == 'A')
		{
			rev_comp_DNA[i_dna_nuc] = 'T';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'T')
		{
			rev_comp_DNA[i_dna_nuc] = 'A';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'G')
		{
			rev_comp_DNA[i_dna_nuc] = 'C';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'C')
		{
			rev_comp_DNA[i_dna_nuc] = 'G';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'a')
		{
			rev_comp_DNA[i_dna_nuc] = 't';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			rev_comp_DNA[i_dna_nuc] = 'a';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'g')
		{
			rev_comp_DNA[i_dna_nuc] = 'c';
			i_dna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'c')
		{
			rev_comp_DNA[i_dna_nuc] = 'g';
			i_dna_nuc++;
		}
		//else if(cur_aln_line[i_nuc] == '-')
		//{
		//	rev_comp_DNA[i_dna_nuc] = '-';
		//	i_dna_nuc++;
		//}
		else
		{
			//printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			rev_comp_DNA[i_dna_nuc] = cur_aln_line[i_nuc];
			i_dna_nuc++;
		}
	} // i_nuc loop.

	//printf("%s", rev_comp_DNA);
	//getc(stdin);
	
	return(rev_comp_DNA);
}

char* get_forward_strand_RNA(char* cur_aln_line, int i_nuc_min, int i_nuc_max)
{
	char* forward_RNA = new char[i_nuc_max - i_nuc_min + 2];
	memset(forward_RNA, 0, sizeof(char) * (i_nuc_max - i_nuc_min + 2));

	//printf("Forward strand RNA for strand:\n");
	//for(int i_nuc = i_nuc_min; i_nuc <= i_nuc_max; i_nuc++)
	//{
	//	printf("%c", cur_aln_line[i_nuc]);
	//}
	//printf("\n");

	int i_rna_nuc = 0;

	// The forward strand RNA is transcribed from 3' end.
	for(int i_nuc = i_nuc_max; i_nuc >= i_nuc_min; i_nuc--)
	{
		if(cur_aln_line[i_nuc] == 'A')
		{
			forward_RNA[i_rna_nuc] = 'U';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'T')
		{
			forward_RNA[i_rna_nuc] = 'A';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'G')
		{
			forward_RNA[i_rna_nuc] = 'C';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'C')
		{
			forward_RNA[i_rna_nuc] = 'G';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'a')
		{
			forward_RNA[i_rna_nuc] = 'u';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 't')
		{
			forward_RNA[i_rna_nuc] = 'a';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'g')
		{
			forward_RNA[i_rna_nuc] = 'c';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == 'c')
		{
			forward_RNA[i_rna_nuc] = 'g';
			i_rna_nuc++;
		}
		else if(cur_aln_line[i_nuc] == '-')
		{
			forward_RNA[i_rna_nuc] = '-';
			i_rna_nuc++;
		}
		else
		{
			printf("Found %c in the alignment line @ %s(%d):\n%s\n",  cur_aln_line[i_nuc], __FILE__, __LINE__, cur_aln_line);
			//exit(0);
			forward_RNA[i_rna_nuc] = cur_aln_line[i_nuc];
			i_rna_nuc++;
		}
	} // i_nuc loop.

	//printf("%s", forward_RNA);
	//getc(stdin);
	
	return(forward_RNA);
}

//void binarize_genome_per_fasta(char* fasta_fp, char* op_dir)
//{
//	// Load fasta file.
//	t_chr_sequence* chr_seq = load_genome_fasta(fasta_fp);
//
//	// Dump the sequence data to the directory.
//	//for(int i_chr = 0; i_chr < chr_seqs->size(); i_chr++)
//	//{		
//	char seq_fp[1000];
//	sprintf(seq_fp, "%s/%s.bin", op_dir, chr_seq->chr_name);
//	int n_cur_chr_nucs = strlen(chr_seq->nucs);
//	fprintf(stderr, "Dumping %d nucs to %s.\n", n_cur_chr_nucs, seq_fp);
//
//	FILE* f_seq = open_f(seq_fp, "wb");
//	fwrite(chr_seq->nucs, n_cur_chr_nucs, 1, f_seq);
//	fclose(f_seq);
//
//	fprintf(stderr, "Done.\n", n_cur_chr_nucs, seq_fp);
//	//} // i_chr loop.
//}
//
//void separate_multi_fasta_file(char* multi_fasta_fp)
//{
//	FILE* f_multi_fasta = open_f(fasta_fp, "r");
//
//	// Read the first line, which is the identifier line.
//	char* cur_line = getline(f_fasta);
//	char bin_fp[1000];
//	sprintf(bin_fp, "%s/%s.bin", bin_dir, &(cur_line[1]));
//	FILE* f_bin = open_f(bin_fp, "wb");
//
//	while(1)
//	{
//		cur_line = getline(f_fasta);
//
//		if(cur_line == NULL)
//		{
//			break;
//		}
//
//		if(cur_line[0] == '>')
//		{
//			printf("Loading %s\n", &(cur_line[1]));
//		}
//		else
//		{
//			// Concatenate the sequence information.
//			fwrite(cur_line, strlen(cur_line), 1, f_bin);
//		}
//	} // file reading loop.
//
//	fclose(f_multi_fasta);
//}

//t_chr_sequence* load_genome_fasta(char* fasta_fp)
//{
//	FILE* f_fasta = open_f(fasta_fp, "r");
//
//	t_chr_sequence* cur_chr_seq = new t_chr_sequence();
//
//    int file_size = 0;
//    fseek(f_fasta, 0, SEEK_END);
//    file_size = (int)ftell(f_fasta);
//    fseek(f_fasta, 0, SEEK_SET);
//
//    printf("Size of %s is around %.2f megabytes\n", fasta_fp, (double)file_size/ (1000.0 * 1000.0));
//
//	cur_chr_seq->nucs = new char[file_size + 3];
//    cur_chr_seq->nucs[0] = 0;
//
//	while(1)
//	{
//		char* cur_line = getline(f_fasta);
//
//		if(cur_line == NULL)
//		{
//			break;
//		}
//
//		if(cur_line[0] == '>')
//		{
//			// Set the chromosome name.
//			cur_chr_seq->chr_name = new char[strlen(cur_line) + 2];
//			strcpy(cur_chr_seq->chr_name, &(cur_line[1]));
//
//			printf("Loading %s\n", cur_chr_seq->chr_name);
//		}
//		else
//		{
//			// Concatenate the sequence information.
//			strcat(cur_chr_seq->nucs, cur_line);
//		}
//	} // file reading loop.
//
//	fclose(f_fasta);
//
//	return(cur_chr_seq);
//}

void assemble_reads_iterative(char* fastq_fp)
{
}

// This function takes a raw sequence list file path, not a formatted file. The formatting should be converted to a sequence list file before calling this function.
//void get_kmer_probability_profiles(char* seq_list_fp, int min_k, int max_k, int max_n_iter)
void get_kmer_probability_profiles(vector<char*>* positive_sequences, vector<char*>* negative_sequences,	// These are the list of sequences to scan.
									vector<int*>* positive_sequence_per_nuc_depths,							// This is the per nucleotide score for the positive sequences.
									int min_k, int max_k, int n_inits, int n_iters_per_init,
									int max_matches_per_sequence,
									char* op_dir) // This is the number of matches to be traced per sequence.
{
	if(positive_sequences == NULL)
	{
		fprintf(stderr, "Must have positive sequences for doing k-mer enrichment.\n");
		exit(0);
	}

	if(negative_sequences != NULL)
	{
		fprintf(stderr, "Using %d positive sequences and %d negative sequences.\n", 
			positive_sequences->size(),
			negative_sequences->size());
	}
	else
	{
		fprintf(stderr, "Using %d positive sequences and no negative sequences.\n", 
			positive_sequences->size());
	}

	// Store the numerized nucleotide sequences for each sequence.
	fprintf(stderr, "Numerizing the positive sequences and setting per nucleotide RD sequence.\n");
	vector<int*>* numerized_pos_sequences = new vector<int*>();
	vector<int>* numerized_pos_seq_lengths = new vector<int>();
	vector<int*>* pos_seq_per_nuc_depths = new vector<int*>();
	int n_pos_skipped_seqs = 0;
	for(int i_seq = 0; i_seq < positive_sequences->size(); i_seq++)
	{
		int* cur_seq_num_seq = new int[t_string::string_length(positive_sequences->at(i_seq))+2];	
		bool is_valid = true;
		for(int i = 0; i < t_string::string_length(positive_sequences->at(i_seq)); i++)
		{
			cur_seq_num_seq[i] = nuc_2_num(positive_sequences->at(i_seq)[i]);
			if(cur_seq_num_seq[i] <= 3 && cur_seq_num_seq[i] >= 0)
			{
			}
			else
			{
				is_valid = false;
				break;
			}
		} // i loop

		if(is_valid)
		{
			numerized_pos_sequences->push_back(cur_seq_num_seq);
			numerized_pos_seq_lengths->push_back(t_string::string_length(positive_sequences->at(i_seq)));

			// Add the per nucleotide sequencing depth if it is supplied.
			if(positive_sequence_per_nuc_depths != NULL)
			{
				pos_seq_per_nuc_depths->push_back(positive_sequence_per_nuc_depths->at(i_seq));
			}
			else
			{
				// Add the sequence depth for the negative set.
				int l_cur_seq = t_string::string_length(positive_sequences->at(i_seq));
				int* cur_seq_depths = new int[l_cur_seq + 2];
				for(int i_nuc = 0; i_nuc <= l_cur_seq; i_nuc++)
				{
					cur_seq_depths[i_nuc] = 1;
				} // i_nuc loop.
				pos_seq_per_nuc_depths->push_back(cur_seq_depths);
			}
		}
		else
		{
			delete [] cur_seq_num_seq;
			fprintf(stderr, "Skipping:\n%s\n", positive_sequences->at(i_seq));
			n_pos_skipped_seqs++;
		}
	} // i_seq loop.

	fprintf(stderr, "Skipped %d positive sequences (%d)\n", n_pos_skipped_seqs, positive_sequences->size());

	// Numerize and check the negative sequences.
	vector<int*>* numerized_neg_sequences = new vector<int*>();
	vector<int>* numerized_neg_seq_lengths = new vector<int>();
	vector<int*>* neg_seq_per_nuc_depths = new vector<int*>();
	int n_skipped_neg_seqs = 0;
	if(negative_sequences == NULL)
	{
		fprintf(stderr, "Skipping the negative sequences.\n");
	}
	else
	{
		for(int i_seq = 0; 
			i_seq < negative_sequences->size(); 
			i_seq++)
		{
			int* cur_seq_num_seq = new int[t_string::string_length(negative_sequences->at(i_seq))+2];	
			bool is_valid = true;
			for(int i = 0; i < t_string::string_length(negative_sequences->at(i_seq)); i++)
			{
				cur_seq_num_seq[i] = nuc_2_num(negative_sequences->at(i_seq)[i]);
				if(cur_seq_num_seq[i] <= 3 && cur_seq_num_seq[i] >= 0)
				{
				}
				else
				{
					is_valid = false;
					break;
				}
			} // i loop

			if(is_valid)
			{
				int l_cur_seq = t_string::string_length(negative_sequences->at(i_seq));
				numerized_neg_sequences->push_back(cur_seq_num_seq);
				numerized_neg_seq_lengths->push_back(t_string::string_length(negative_sequences->at(i_seq)));

				// Add the sequence depth for the negative set.
				int* cur_seq_depths = new int[l_cur_seq + 2];
				for(int i_nuc = 0; i_nuc <= l_cur_seq; i_nuc++)
				{
					cur_seq_depths[i_nuc] = 1;
				} // i_nuc loop.
				neg_seq_per_nuc_depths->push_back(cur_seq_depths);
			}
			else
			{
				delete [] cur_seq_num_seq;
				fprintf(stderr, "Skipping:\n%s\n", negative_sequences->at(i_seq));
				n_skipped_neg_seqs++;
			}
		} // i_seq loop.

		fprintf(stderr, "Skipped %d sequences (%d)\n", n_skipped_neg_seqs, negative_sequences->size());
	} // existence check for the negative sequences.

	// k-mer discovery loop.
	// Initialize the k-mer probability profiles.
	double*** kmer_emission_log_probs = new double**[max_k-min_k+1];
	double*** kmer_observation_log_counts = new double**[max_k-min_k+1];
	kmer_emission_log_probs -= min_k;
	kmer_observation_log_counts -= min_k;
	for(int k = min_k; k <= max_k; k++)
	{
		kmer_emission_log_probs[k] = new double*[k+2];
		kmer_observation_log_counts[k] = new double*[k+2];

		// pos = 0 is the beginning position.
		for(int pos = 0; pos <= k; pos++)
		{
			// Allocate a number for each nucleotide.
			kmer_emission_log_probs[k][pos] = new double[6];
			kmer_observation_log_counts[k][pos] = new double[6];

			//kmer_emission_log_probs[k][pos][0] = xlog(0.0);
			//kmer_emission_log_probs[k][pos][1] = xlog(0.0);
			//kmer_emission_log_probs[k][pos][2] = xlog(0.0);
			//kmer_emission_log_probs[k][pos][3] = xlog(0.0);

			//kmer_observation_log_counts[k][pos][0] = xlog(0.0);
			//kmer_observation_log_counts[k][pos][1] = xlog(0.0);
			//kmer_observation_log_counts[k][pos][2] = xlog(0.0);
			//kmer_observation_log_counts[k][pos][3] = xlog(0.0);

			kmer_emission_log_probs[k][pos][0] = (0.0);
			kmer_emission_log_probs[k][pos][1] = (0.0);
			kmer_emission_log_probs[k][pos][2] = (0.0);
			kmer_emission_log_probs[k][pos][3] = (0.0);
			kmer_emission_log_probs[k][pos][4] = (0.0);

			kmer_observation_log_counts[k][pos][0] = (0.0);
			kmer_observation_log_counts[k][pos][1] = (0.0);
			kmer_observation_log_counts[k][pos][2] = (0.0);
			kmer_observation_log_counts[k][pos][3] = (0.0);
			kmer_observation_log_counts[k][pos][4] = (0.0);
		} // pos loop.
	} // k loop.

	// 7 Different probability profiles for a position.
	const int n_init_prob_sets = 6;
	double init_probs[n_init_prob_sets][4] = {.23, .27, .24, .26,		// Almost uniform.
												.25, .25, .25, .25,		// Totally uniform.
												.50, .25, .125, .125,	// One highest, one follows, two same and low.
												.10, .20, .30, .40,		// Increasing/Decreasing.
												.05, .225, .325, .40,	// Steeper Increasing/Decreasing.
												.40, .40, .10, .10};	// Two high, two low.

	// Start the inits.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for(int i_init = 0; i_init < n_inits; i_init++)
	{
		// For all the k-mers, initialize the current emission probabilities.
		char cur_init_probs_fp[1000];
		sprintf(cur_init_probs_fp, "%s/init_%d_probs.list", op_dir, i_init);
		FILE* f_cur_init_probs = open_f(cur_init_probs_fp, "w");
		for(int k = min_k; k <= max_k; k++)
		{
			fprintf(f_cur_init_probs, "k\t%d\n", k);
			for(int pos = 0; pos <= k; pos++)
			{
				// Select a random probability set from the predefined sets.
				int cur_prob_set_i = (int)(floor(rng->random_double_ran3() * n_init_prob_sets)) % n_init_prob_sets;

				// Choose a shuffle for the current probability set.
				vector<int>* cur_shuff = rng->permute_indices(4, 4);

				// Randomly set the probabilities.
				kmer_emission_log_probs[k][pos][0] = init_probs[cur_prob_set_i][cur_shuff->at(0)];
				kmer_emission_log_probs[k][pos][1] = init_probs[cur_prob_set_i][cur_shuff->at(1)];
				kmer_emission_log_probs[k][pos][2] = init_probs[cur_prob_set_i][cur_shuff->at(2)];
				kmer_emission_log_probs[k][pos][3] = init_probs[cur_prob_set_i][cur_shuff->at(3)];

				fprintf(f_cur_init_probs, "%lf\t%lf\t%lf\t%lf\n", kmer_emission_log_probs[k][pos][0], 
																	kmer_emission_log_probs[k][pos][1], 
																	kmer_emission_log_probs[k][pos][2], 
																	kmer_emission_log_probs[k][pos][3]);

				delete cur_shuff;
			} // pos loop.
		} // k loop.
		fclose(f_cur_init_probs);

		// Start the training iterations for the current initialization.
		for(int i_iter = 0; i_iter < n_iters_per_init; i_iter++)
		{
			fprintf(stderr, "Started %d. iteration.\n", i_iter);
			//getc(stdin);

			// Re-initialize the current counts to 0 over the profiles.
			for(int k = min_k; k <= max_k; k++)
			{
				for(int pos = 0; pos < k; pos++)
				{
					//kmer_observation_log_counts[k][pos][0] = xlog(0);
					//kmer_observation_log_counts[k][pos][1] = xlog(0);
					//kmer_observation_log_counts[k][pos][2] = xlog(0);
					//kmer_observation_log_counts[k][pos][3] = xlog(0);
					kmer_observation_log_counts[k][pos][0] = (0);
					kmer_observation_log_counts[k][pos][1] = (0);
					kmer_observation_log_counts[k][pos][2] = (0);
					kmer_observation_log_counts[k][pos][3] = (0);
				} // pos loop.
			} // k loop.

			// Update the motif counts with positive and negative sequences.
			fprintf(stderr, "Updating counts with positive sequences.\n");
			update_motif_counts_per_sequence_list(numerized_pos_sequences, 
				numerized_pos_seq_lengths,
				pos_seq_per_nuc_depths,
				min_k, max_k, 
				kmer_emission_log_probs,
				kmer_observation_log_counts,
				true);

			// Update the motif counts with negative sequences.
			if(numerized_neg_sequences->size() > 0)
			{
				fprintf(stderr, "Updating counts with negative sequences.\n");
				update_motif_counts_per_sequence_list(numerized_neg_sequences, 
					numerized_neg_seq_lengths,
					neg_seq_per_nuc_depths,
					min_k, max_k, 
					kmer_emission_log_probs,
					kmer_observation_log_counts,
					false);
			}
			else
			{
				fprintf(stderr, "Skipping negative sequence based updating of the counts.\n");
			}

			fprintf(stderr, "Updating the emission probabilities.\n");
			//getc(stdin);
			char cur_init_profiles_fp[1000];
			sprintf(cur_init_profiles_fp, "%s/motifs_per_iterations_init_%d.txt", op_dir, i_init);
			FILE* f_motifs_per_iteration = open_f(cur_init_profiles_fp, "a");
			for(int k = min_k; k <= max_k; k++)
			{
				// For each position, compute the frequencies.
				fprintf(f_motifs_per_iteration, "Iteration\t%d\n", i_iter);
				fprintf(f_motifs_per_iteration, "k: %d\n", k);
				for(int pos = 0; pos < k; pos++)
				{
					// Sanity check.
					//double total_posn_log_counts = xlog_sum(xlog_sum(xlog_sum(kmer_observation_log_counts[k][pos][0], kmer_observation_log_counts[k][pos][1]), kmer_observation_log_counts[k][pos][2]), kmer_observation_log_counts[k][pos][3]);
					double total_posn_log_counts = (kmer_observation_log_counts[k][pos][0] +
													kmer_observation_log_counts[k][pos][1] +
													kmer_observation_log_counts[k][pos][2] + 
													kmer_observation_log_counts[k][pos][3]);

					//// Check if one of the counts is 0, update it to the 1/10th of the smallest count.
					//int non_zero_min_count = 1000*1000;
					//for(int i_nuc = 0; i_nuc < 4; i_nuc++)
					//{
					//	if(non_zero_min_count > kmer_observation_log_counts[k][pos][i_nuc] &&
					//		kmer_observation_log_counts[k][pos][i_nuc] > 0)
					//	{
					//		non_zero_min_count = kmer_observation_log_counts[k][pos][i_nuc];
					//	}
					//} // i loop.

					//// Sanity check on the value of next min.
					//if(non_zero_min_count == 0)
					//{
					//	fprintf(stderr, "Next min value is 0.\n");
					//	exit(0);
					//}

					//// Go over all the counts.
					//int pseudocount = (non_zero_min_count / 10 > 1)?(non_zero_min_count / 10):(1);
					//for(int i_nuc = 0; i_nuc < 4; i_nuc++)
					//{
					//	// Is this value 0?
					//	if(kmer_observation_log_counts[k][pos][i_nuc] == 0)
					//	{
					//		kmer_observation_log_counts[k][pos][i_nuc] = pseudocount;
					//	}
					//} // i loop.

					////if(!xlog_comp(xlog(sequences->size()), total_posn_log_counts))
					//if((numerized_pos_sequences->size()+numerized_neg_sequences->size()) != total_posn_log_counts)
					//{
					//	fprintf(stderr, "The counts do not add upto the sequence number: %d, %lf\n", 
					//		(numerized_pos_sequences->size()+numerized_neg_sequences->size()), 
					//		total_posn_log_counts);
					//}

					// Update the probabilities: These are log probabilities.
					//kmer_emission_log_probs[k][pos][0] = xlog_div(kmer_observation_log_counts[k][pos][0], total_posn_log_counts);
					//kmer_emission_log_probs[k][pos][1] = xlog_div(kmer_observation_log_counts[k][pos][1], total_posn_log_counts);
					//kmer_emission_log_probs[k][pos][2] = xlog_div(kmer_observation_log_counts[k][pos][2], total_posn_log_counts);
					//kmer_emission_log_probs[k][pos][3] = xlog_div(kmer_observation_log_counts[k][pos][3], total_posn_log_counts);
					kmer_emission_log_probs[k][pos][0] = (kmer_observation_log_counts[k][pos][0] / total_posn_log_counts);
					kmer_emission_log_probs[k][pos][1] = (kmer_observation_log_counts[k][pos][1] / total_posn_log_counts);
					kmer_emission_log_probs[k][pos][2] = (kmer_observation_log_counts[k][pos][2] / total_posn_log_counts);
					kmer_emission_log_probs[k][pos][3] = (kmer_observation_log_counts[k][pos][3] / total_posn_log_counts);

					// Sanity check on the number of entries at 4th position, which is not an existing value.
					if(kmer_observation_log_counts[k][pos][4] > 0)
					{
						fprintf(stderr, "Sanity check failed: 4th entry has non-zero values.");
						exit(0);
					}


					// Dump the profile.
					fprintf(f_motifs_per_iteration, "%.2f\t%.2f\t%.2f\t%.2f\n", kmer_emission_log_probs[k][pos][0], kmer_emission_log_probs[k][pos][1], kmer_emission_log_probs[k][pos][2], kmer_emission_log_probs[k][pos][3]);
				} // pos loop.
			} // k loop.
			fclose(f_motifs_per_iteration);

			// Must add a maximum iteration number and a convergence check for the probabilities.
		} // i_iter loop.

		// Dump the final probabilities for each k-mer.
		for(int k = min_k; k <= max_k; k++)
		{
			char cur_k_final_freqs_fp[1000];
			sprintf(cur_k_final_freqs_fp, "%s/k_%d_init_%d_freqs.list", op_dir, k, i_init);
			FILE* f_cur_k_final_freqs = open_f(cur_k_final_freqs_fp, "w");
			for(int pos = 0; pos < k; pos++)
			{
				fprintf(f_cur_k_final_freqs, "%lf\t%lf\t%lf\t%lf\n", kmer_emission_log_probs[k][pos][0], 
																	kmer_emission_log_probs[k][pos][1], 
																	kmer_emission_log_probs[k][pos][2], 
																	kmer_emission_log_probs[k][pos][3]);

			} // pos loop.
			fclose(f_cur_k_final_freqs);
		} // k loop.
	} // i_init loop.

	delete rng;

	// Must report a match score between the identified motifs and the given positive and negative sequences like a likelihood. This is used to sort the motifs.
}

// Generate random k-mers and build the score distribution for random sequences.
void get_shuffled_motif_matching_statistics(char* sequence, double** kmer_emission_probs, int k, double& mean, double& std_dev, int n_randomizations)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	int l_seq = t_string::string_length(sequence);

	char* shuffled_seq = new char[l_seq + 2];
	vector<double>* all_scores = new vector<double>();
	double* cur_seq_negtv_match_scores = new double[l_seq+2];
	double* cur_seq_postv_match_scores = new double[l_seq+2];
	for(int i_rand = 0; i_rand < n_randomizations; i_rand++)
	{
		// Shuffle the sequence.
		vector<int>* shuffled_seq_i = rng->fast_permute_indices(0, l_seq);
		for(int i_nuc = 0; i_nuc < l_seq; i_nuc++)
		{
			shuffled_seq[i_nuc] = sequence[shuffled_seq_i->at(i_nuc)];
		} // i_nuc loop.
		shuffled_seq[l_seq] = 0;

		// For the current shuffled sequence, get the score.
		scan_motif_in_sequence(shuffled_seq, kmer_emission_probs, k, cur_seq_postv_match_scores, cur_seq_negtv_match_scores);

		// Find the maximum matching score over two strands.
		double maximum_score = 0.0;
		for(int i_nuc = 0; i_nuc < l_seq-k; i_nuc++)
		{
			if(maximum_score < cur_seq_postv_match_scores[i_nuc])
			{
				maximum_score = cur_seq_postv_match_scores[i_nuc];
			}
		} // i_nuc loop.

		for(int i_nuc = 0; i_nuc < l_seq-k; i_nuc++)
		{
			if(maximum_score < cur_seq_negtv_match_scores[i_nuc])
			{
				maximum_score = cur_seq_negtv_match_scores[i_nuc];
			}
		} // i_nuc loop.

		// Add the maximum score.
		all_scores->push_back(maximum_score);

		// Free memory.
		delete shuffled_seq_i;
	} // i_rand loop.

	if(n_randomizations != all_scores->size())
	{
		fprintf(stderr, "Shuffled scores contain %d elements.\n", all_scores->size());
		exit(0);
	}

	delete [] shuffled_seq;
	delete rng;

	// Compute the mean and standard deviation.
	mean = 0.0;
	std_dev = 0.0;
	for(int i_score = 0; i_score < all_scores->size(); i_score++)
	{
		mean += all_scores->at(i_score);
	} // i_score loop.
	mean /= all_scores->size();

	for(int i_score = 0; i_score < all_scores->size(); i_score++)
	{
		std_dev += (all_scores->at(i_score) - mean) * (all_scores->at(i_score) - mean);
	} // i_score loop.
	std_dev /= all_scores->size();

	std_dev = pow(std_dev, .5);

	delete [] cur_seq_postv_match_scores;
	delete [] cur_seq_negtv_match_scores;
	delete all_scores;
}

void scan_motif_in_sequence(char* sequence, double** kmer_emission_probs, int k, double* postv_match_score, double* negtv_match_score)
{
	int l_seq = t_string::string_length(sequence);
	double* match_score_per_posn = new double[l_seq + 2];
	for(int profile_posn = 0; profile_posn < l_seq; profile_posn++)
	{
		postv_match_score[profile_posn] = 0.0;
		negtv_match_score[profile_posn] = 0.0;
	} // profile_posn loop.

	// Slide over all the positions in this sequence.
	double optimal_log_score = 0.0;
	int optimizing_scoring_profile_posn = 0;
	for(int strand = 0; strand <= 1; strand++)
	{
		for(int profile_posn = 0; profile_posn < l_seq-k; profile_posn++)
		{
			//double cur_profile_posn_log_score = xlog(1.0);
			double cur_profile_posn_log_score = (1.0);

			// Compute the current matching score.
			if(strand == 1)
			{
				// This is the negative strand. profile_posn holds the leftmost position of the profile in sequence.
				// posn holds the leftmost position of the nucleotide in the profile.
				// For negative strand, the nucleotide needs to be complemented and the profile position is in the opposite direction.
				for(int posn = 0; posn < k; posn++)
				{
					int cur_nuc_num = 3 - nuc_2_num(sequence[profile_posn+posn]);
					//int cur_nuc_num = numerized_sequences->at(i_seq)[profile_posn+posn];
					//cur_profile_posn_log_score = xlog_mul(cur_profile_posn_log_score, kmer_emission_log_probs[k][posn][cur_nuc_num]);
					cur_profile_posn_log_score = (cur_profile_posn_log_score * kmer_emission_probs[k-posn-1][cur_nuc_num]);
				} // posn loop.

				//match_score_per_posn[profile_posn] = cur_profile_posn_log_score;
				negtv_match_score[profile_posn] = cur_profile_posn_log_score;
			}
			else
			{
				for(int posn = 0; posn < k; posn++)
				{
					int cur_nuc_num = nuc_2_num(sequence[profile_posn+posn]);
					//int cur_nuc_num = numerized_sequences->at(i_seq)[profile_posn+posn];
					//cur_profile_posn_log_score = xlog_mul(cur_profile_posn_log_score, kmer_emission_log_probs[k][posn][cur_nuc_num]);
					cur_profile_posn_log_score = (cur_profile_posn_log_score * kmer_emission_probs[posn][cur_nuc_num]);
				} // posn loop.

				//match_score_per_posn[profile_posn] = cur_profile_posn_log_score;
				postv_match_score[profile_posn] = cur_profile_posn_log_score;
			}
		} // profile_posn loop.
	} // strand loop.
	//return(match_score_per_posn);
}

bool sort_kmers_per_count(t_kmer_freq_info* info1, t_kmer_freq_info* info2)
{
	return(info1->count > info2->count);
}

// Get k-mer frequencies and do an exhaustive enumeration of the most 
void get_kmer_frequency_enrichments(vector<char*>* positive_sequences, vector<char*>* negative_sequences, 
								vector<int*>* positive_sequence_per_nuc_depths,
								int min_k, int max_k)
{
	if(positive_sequences == NULL)
	{
		fprintf(stderr, "Must have positive sequences for doing k-mer enrichment.\n");
		exit(0);
	}

	if(negative_sequences != NULL)
	{
		fprintf(stderr, "Using %d positive sequences and %d negative sequences.\n", 
			positive_sequences->size(),
			negative_sequences->size());
	}
	else
	{
		fprintf(stderr, "Using %d positive sequences and no negative sequences.\n", 
			positive_sequences->size());
	}

	// Store the numerized nucleotide sequences for each sequence.
	fprintf(stderr, "Numerizing the positive sequences and setting per nucleotide RD sequence.\n");
	vector<int*>* numerized_pos_sequences = new vector<int*>();
	vector<int>* numerized_pos_seq_lengths = new vector<int>();
	vector<int*>* pos_seq_per_nuc_depths = new vector<int*>();
	int n_pos_skipped_seqs = 0;
	for(int i_seq = 0; i_seq < positive_sequences->size(); i_seq++)
	{
		int* cur_seq_num_seq = new int[t_string::string_length(positive_sequences->at(i_seq))+2];	
		bool is_valid = true;
		for(int i = 0; i < t_string::string_length(positive_sequences->at(i_seq)); i++)
		{
			cur_seq_num_seq[i] = nuc_2_num(positive_sequences->at(i_seq)[i]);
			if(cur_seq_num_seq[i] <= 3 && cur_seq_num_seq[i] >= 0)
			{
			}
			else
			{
				is_valid = false;
				break;
			}
		} // i loop

		if(is_valid)
		{
			numerized_pos_sequences->push_back(cur_seq_num_seq);
			numerized_pos_seq_lengths->push_back(t_string::string_length(positive_sequences->at(i_seq)));

			// Add the per nucleotide sequencing depth if it is supplied.
			if(positive_sequence_per_nuc_depths != NULL)
			{
				pos_seq_per_nuc_depths->push_back(positive_sequence_per_nuc_depths->at(i_seq));
			}
			else
			{
				// Add the sequence depth for the negative set.
				int l_cur_seq = t_string::string_length(positive_sequences->at(i_seq));
				int* cur_seq_depths = new int[l_cur_seq + 2];
				for(int i_nuc = 0; i_nuc <= l_cur_seq; i_nuc++)
				{
					cur_seq_depths[i_nuc] = 1;
				} // i_nuc loop.
				pos_seq_per_nuc_depths->push_back(cur_seq_depths);
			}
		}
		else
		{
			delete [] cur_seq_num_seq;
			fprintf(stderr, "Skipping:\n%s\n", positive_sequences->at(i_seq));
			n_pos_skipped_seqs++;
		}
	} // i_seq loop.

	fprintf(stderr, "Skipped %d positive sequences (%d)\n", n_pos_skipped_seqs, positive_sequences->size());

	// Numerize and check the negative sequences.
	vector<int*>* numerized_neg_sequences = new vector<int*>();
	vector<int>* numerized_neg_seq_lengths = new vector<int>();
	vector<int*>* neg_seq_per_nuc_depths = new vector<int*>();
	int n_skipped_neg_seqs = 0;
	if(negative_sequences == NULL)
	{
		fprintf(stderr, "Skipping the negative sequences.\n");
	}
	else
	{
		for(int i_seq = 0; 
			i_seq < negative_sequences->size(); 
			i_seq++)
		{
			int* cur_seq_num_seq = new int[t_string::string_length(negative_sequences->at(i_seq))+2];	
			bool is_valid = true;
			for(int i = 0; i < t_string::string_length(negative_sequences->at(i_seq)); i++)
			{
				cur_seq_num_seq[i] = nuc_2_num(negative_sequences->at(i_seq)[i]);
				if(cur_seq_num_seq[i] <= 3 && cur_seq_num_seq[i] >= 0)
				{
				}
				else
				{
					is_valid = false;
					break;
				}
			} // i loop

			if(is_valid)
			{
				int l_cur_seq = t_string::string_length(negative_sequences->at(i_seq));
				numerized_neg_sequences->push_back(cur_seq_num_seq);
				numerized_neg_seq_lengths->push_back(t_string::string_length(negative_sequences->at(i_seq)));

				// Add the sequence depth for the negative set.
				int* cur_seq_depths = new int[l_cur_seq + 2];
				for(int i_nuc = 0; i_nuc <= l_cur_seq; i_nuc++)
				{
					cur_seq_depths[i_nuc] = 1;
				} // i_nuc loop.
				neg_seq_per_nuc_depths->push_back(cur_seq_depths);
			}
			else
			{
				delete [] cur_seq_num_seq;
				fprintf(stderr, "Skipping:\n%s\n", negative_sequences->at(i_seq));
				n_skipped_neg_seqs++;
			}
		} // i_seq loop.

		fprintf(stderr, "Skipped %d sequences (%d)\n", n_skipped_neg_seqs, negative_sequences->size());
	} // existence check for the negative sequences.

	vector<t_kmer_freq_info*>* valid_raw_enriched_motifs = new vector<t_kmer_freq_info*>();
	vector<t_kmer_freq_info*>* cur_kmer_valid_raw_enriched_motifs = new vector<t_kmer_freq_info*>();
	for(int k = min_k; k <= max_k; k++)
	{
		t_kmer_count_tree_entry* main_entry = build_kmer_frequency_tree(numerized_pos_sequences, numerized_pos_seq_lengths, k);

		// Remove the entries that are exact reverse complement of each other: These are redundant.

		// Do we have negative sequences?
		if(numerized_neg_sequences == NULL || numerized_neg_sequences->size() == 0)
		{
			// Get the nucleotide frequencies.
			double nuc_freqs[6];
			memset(nuc_freqs, 0, sizeof(double)*6);
			for(int i_seq = 0; i_seq < numerized_pos_sequences->size(); i_seq++)
			{
				for(int i_nuc = 0; i_nuc < numerized_pos_seq_lengths->at(i_seq); i_nuc++)
				{
					int cur_nuc = numerized_pos_sequences->at(i_seq)[i_nuc];
					nuc_freqs[cur_nuc]++;
					nuc_freqs[3-cur_nuc]++;
				} // i_nuc loop.
			} // i_seq loop.

			double total_n_nucs = (nuc_freqs[0] + nuc_freqs[1] + nuc_freqs[2] + nuc_freqs[3]);

			// Get the nucleotide frequencies to build a simple null model which doesnt work well.
			nuc_freqs[0] = nuc_freqs[0] / total_n_nucs;
			nuc_freqs[1] = nuc_freqs[1] / total_n_nucs;
			nuc_freqs[2] = nuc_freqs[2] / total_n_nucs;
			nuc_freqs[3] = nuc_freqs[3] / total_n_nucs;

			// Sort the entries with respect to frequency, then align the most frequent motifs to identify motifs with highest likelihood.
			vector<t_kmer_freq_info*>* kmer_list = new vector<t_kmer_freq_info*>();
			char cur_kmer[100];
			memset(cur_kmer, 0, sizeof(char) * 100);
			get_kmer_count_array(kmer_list, main_entry, cur_kmer);

			// Sort the kmers with respect to the k-mer count.
			sort(kmer_list->begin(), kmer_list->end(), sort_kmers_per_count);

			// Compute the number of times we expect to see the kmers.
			FILE* f_sorted_kmer_list = open_f("sorted_kmer_list.txt", "w");
			fprintf(stderr, "Dumping %d kmers.\n", kmer_list->size());
			for(int i_kmer = 0; i_kmer < kmer_list->size(); i_kmer++)
			{
				char* cur_kmer = kmer_list->at(i_kmer)->kmer;
				double null_prob = 1.0;
				for(int i = 0; i < k; i++)
				{
					null_prob *= nuc_freqs[nuc_2_num(cur_kmer[i])];
				} // i loop.

				double n_expected = null_prob * (total_n_nucs - k);

				if(kmer_list->at(i_kmer)->count > n_expected)
					fprintf(f_sorted_kmer_list, "%s\t%d\t%lf\n", cur_kmer, kmer_list->at(i_kmer)->count, n_expected);
			} // i_kmer loop.
			fclose(f_sorted_kmer_list);

			fprintf(stderr, "Done!\n");
		}
		else if(numerized_neg_sequences->size() > 0)
		{
			fprintf(stderr, "Building the negative list frequency tree.\n");

			// Build the count tree for the negative sequences.
			t_kmer_count_tree_entry* neg_main_entry = build_kmer_frequency_tree(numerized_neg_sequences, numerized_neg_seq_lengths, k);

			int n_total_neg_kmers = get_total_kmer_count(neg_main_entry, k);
			int n_total_pos_kmers = get_total_kmer_count(main_entry, k);

			fprintf(stderr, "%d positive set kmers and %d negative set kmers.\n", n_total_pos_kmers, n_total_neg_kmers);

			// Set the background counts for the positive kmers.
			get_bckgrnd_counts(main_entry, neg_main_entry);

			// Now we can delete the negative k-mer frequency tree.
			delete_kmer_freq_tree(neg_main_entry);

			vector<t_kmer_freq_info*>* kmer_list = new vector<t_kmer_freq_info*>();
			char cur_kmer[100];
			memset(cur_kmer, 0, sizeof(char) * 100);
			get_kmer_count_array(kmer_list, main_entry, cur_kmer); 

			// Sort the negative kmer list.
			sort(kmer_list->begin(), kmer_list->end(), sort_kmers_per_count);

			//FILE* f_sorted_kmer_list = open_f("sorted_kmer_list.txt", "w");
			char cur_sorted_kmer_list_fp[1000];
			sprintf(cur_sorted_kmer_list_fp, "sorted_%dmer_list.txt", k);
			FILE* f_sorted_kmer_list = open_f(cur_sorted_kmer_list_fp, "w");
			fprintf(stderr, "Computing the significance for %d kmers.\n", kmer_list->size());
			for(int i_kmer = 0; i_kmer < kmer_list->size(); i_kmer++)
			{
				char* cur_kmer = kmer_list->at(i_kmer)->kmer;

				double pos_fraction = ((double)kmer_list->at(i_kmer)->count / (double)n_total_pos_kmers);
				double neg_fraction = ((double)kmer_list->at(i_kmer)->bckgrnd_count / (double)n_total_neg_kmers);
			
				int cur_kmer_numseq[1000];
				for(int i = 0; i < t_string::string_length(kmer_list->at(i_kmer)->kmer); i++)
				{
					cur_kmer_numseq[i] = nuc_2_num(kmer_list->at(i_kmer)->kmer[i]);
				}

				if(pos_fraction > neg_fraction &&
					check_kmer_composition(cur_kmer_numseq, 0, k))
				{
					fprintf(f_sorted_kmer_list, "%s\t%d\t%lf\t%lf\n", cur_kmer, kmer_list->at(i_kmer)->count, pos_fraction, neg_fraction);
					cur_kmer_valid_raw_enriched_motifs->push_back(kmer_list->at(i_kmer));
				}
				else
				{
					// Free the memory for the current kmer.
					delete kmer_list->at(i_kmer);
				}
			} // i_kmer loop.

			fclose(f_sorted_kmer_list);
		} // negative sequence list existence check.

		// Add the valid kmers to the list of valid kmers.
		const int MAX_KMER_PER_KMER_LENGTH = 2000;
		for(int i_kmer = 0; 
			i_kmer < MAX_KMER_PER_KMER_LENGTH && i_kmer < cur_kmer_valid_raw_enriched_motifs->size(); 
			i_kmer++)
		{
			valid_raw_enriched_motifs->push_back(cur_kmer_valid_raw_enriched_motifs->at(i_kmer));
		} // i_kmer loop.
		
		fprintf(stderr, "Computing the distance between the valid k-mers.\n");
		double distance_buffer[100][100];
		//double score_buffer[100][100];
		double** distances = new double*[valid_raw_enriched_motifs->size() + 2];
		for(int i_kmer = 0; i_kmer <= valid_raw_enriched_motifs->size(); i_kmer++)
		{
			distances[i_kmer] = new double[valid_raw_enriched_motifs->size()+2];
		} // i_kmer loop.

		// Loop over all the kmers.
		const double GAP_DISTANCE = 1;
		const double MISMATCH_DISTANCE = 2;
		const double MATCH_DISTANCE = 0;
		fprintf(stderr, "Computing edit distances between the k-mers with mismatch distance of %lf and gap distance of %lf\n", MISMATCH_DISTANCE, GAP_DISTANCE);
		FILE* f_dists = open_f("distances.list", "w");
		FILE* f_avg_dists = open_f("per_kmer_avg_dists.txt", "w");
		for(int i_kmer = 0; i_kmer < valid_raw_enriched_motifs->size(); i_kmer++)
		{
			//distances[i_kmer] = new double[valid_raw_enriched_motifs->size() + 2];
			char* kmer1 = valid_raw_enriched_motifs->at(i_kmer)->kmer;
			for(int j_kmer = 0; j_kmer < valid_raw_enriched_motifs->size(); j_kmer++)
			{
				char* kmer2 = valid_raw_enriched_motifs->at(j_kmer)->kmer;

				// Go over all the values and find the max value.
				int l_kmer1 = t_string::string_length(kmer1);
				int l_kmer2 = t_string::string_length(kmer2);

				//score_buffer[0][0] = 0;
				//score_buffer[0][0] += (kmer1[0]==kmer2[0]);
				//for(int i = 0; i < l_kmer1; i++)
				//{
				//	score_buffer[i][0] = (kmer1[i]==kmer2[0]);
				//} //i loop.

				//for(int j = 0; j < l_kmer2; j++)
				//{
				//	score_buffer[0][j] = (kmer1[0]==kmer2[j]);
				//} // j loop

				//// Go over all the values and find the max value.
				//int cur_kmers_match_score = 0;
				//int max_match_i = 0;
				//int max_match_j = 0;
				//for(int i = 1; i < l_kmer1; i++)
				//{
				//	for(int j = 1; j < l_kmer2; j++)
				//	{
				//		score_buffer[i][j] = score_buffer[i-1][j-1] + (kmer1[i]==kmer2[j]);

				//		// Update the score.
				//		if(cur_kmers_match_score < score_buffer[i][j])
				//		{
				//			cur_kmers_match_score = score_buffer[i][j];
				//			max_match_i = i;
				//			max_match_j = j;
				//		}
				//	} // j loop.
				//} // i loop.


				distance_buffer[0][0] = 0;

				// Initialize all the insertion starts.
				for(int i = 0; i <= k; i++)
				{
					distance_buffer[i][0] = i * GAP_DISTANCE;
				} //i loop.

				for(int j = 0; j <= k; j++)
				{
					distance_buffer[0][j] = j * GAP_DISTANCE;
				} // j loop
				for(int i = 1; i <= k; i++)
				{
					for(int j = 1; j <= k; j++)
					{
						int seq_i = i-1;
						int seq_j = j-1;

						// Find the minimum distance between the strings.
						distance_buffer[i][j] = distance_buffer[i-1][j-1] + MATCH_DISTANCE * (kmer1[seq_i]==kmer2[seq_j]) + MISMATCH_DISTANCE * (kmer1[seq_i]!=kmer2[seq_j]);

						// Add gaps at the end of the sequences.
						if(j == k)
						{
							distance_buffer[i][j] = MIN(distance_buffer[i][j], distance_buffer[i-1][j] + GAP_DISTANCE);
						}

						if(i == k)
						{
							distance_buffer[i][j] = MIN(distance_buffer[i][j], distance_buffer[i][j-1] + GAP_DISTANCE);
						}

						//// Update the score.
						//if(cur_kmers_distace > distance_buffer[i][j])
						//{
						//	cur_kmers_distace = distance_buffer[i][j];
						//	max_match_i = i;
						//	max_match_j = j;
						//}
					} // j loop.
				} // i loop.

if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
				fprintf(stderr, "%s\t%s\t%lf\n", valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer, distance_buffer[k][k]);
				getc(stdin);
}

				//// Update the distance.
				//distances[i_kmer+1][j_kmer+1] = 1;
				//if(cur_kmers_match_score > 7)
				//{
				//	distances[i_kmer+1][j_kmer+1] = 0;
				//}

				distances[i_kmer+1][j_kmer+1] = distance_buffer[k][k];
				fprintf(f_dists, "%lf\t", distance_buffer[k][k]);
			} // j_kmer loop.

			fprintf(f_dists, "\n");

			double cur_avg_dist = 0.0;			
			for(int i = 0; i < valid_raw_enriched_motifs->size(); i++)
			{
				cur_avg_dist += distances[i_kmer+1][i+1];
			} // i loop/
			fprintf(f_avg_dists, "%s\t%lf\n", valid_raw_enriched_motifs->at(i_kmer)->kmer, cur_avg_dist/valid_raw_enriched_motifs->size());
		} // i_kmer loop.
		fclose(f_avg_dists);

		fclose(f_dists);

		exit(0);

		//// Cluster the k-mers for this length.
		//t_matrix* cur_kmers_distance_matrix = new t_matrix(distances, valid_raw_enriched_motifs->size(), valid_raw_enriched_motifs->size(), false);
		//t_obj_set* cur_kmers_obj_list = new t_obj_set(valid_raw_enriched_motifs->size(), cur_kmers_distance_matrix);
		//cur_kmers_obj_list->diana();

		//// Write the 10 clustering scheme indices.
		//int n_clusters = cur_kmers_obj_list->cluster_hierarchy->at(10)->size();

		//for(int i_clust = 0; i_clust < n_clusters; i_clust++)
		//{
		//	//cur_kmers_obj_list->cluster_hierarchy->at(10)->at(0)->elements->at(0)->dist_matrix_i
		//}

		// Select the height for the clustering, then dump the clusterings.

		// Clear the list of valid kmers.
		cur_kmer_valid_raw_enriched_motifs->clear();

		// Free memory for the current frequency tree for the positive sequences.
		delete_kmer_freq_tree(main_entry);

		// Free distance memory.
		for(int i_kmer = 0; i_kmer <= valid_raw_enriched_motifs->size(); i_kmer++)
		{
			delete [] distances[i_kmer];
		} // i_kmer loop.
	} // k loop.

	fprintf(stderr, "%d valid raw enriched kmers are left, computing the distance between the valid k-mers.\n", valid_raw_enriched_motifs->size());

	exit(0);

	// Before starting the matches, go over each k-mer and find the position of it over the sequences so that we can refer back to the sequences while building the PWMs.
	fprintf(stderr, "Setting the matches for each k-mer.\n");
	for(int i_kmer = 0; i_kmer < valid_raw_enriched_motifs->size(); i_kmer++)
	{
		int cur_k = t_string::string_length(valid_raw_enriched_motifs->at(i_kmer)->kmer);
		valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns = new vector<t_k_mer_match_posn_info*>();

		// For the current k-mer, set the position information on the positive sequences.
		for(int i_seq = 0; i_seq < numerized_pos_sequences->size(); i_seq++)
		{
			int l_cur_seq = numerized_pos_seq_lengths->at(i_seq);
			int* cur_seq_nums = numerized_pos_sequences->at(i_seq);

			for(int strand = 0; strand <= 1; strand++)
			{
				// Check if the current k-mer matches this sequence.
				for(int left_kmer_posn = 0; left_kmer_posn <= l_cur_seq-cur_k; left_kmer_posn++)
				{
					bool cur_kmer_match = true;
					for(int kmer_i = 0; kmer_i < cur_k; kmer_i++)
					{
						int seq_kmer_i = kmer_i;
						if(strand == 1)
						{
							seq_kmer_i = cur_k - seq_kmer_i - 1;
						}

						// Get the sequence number.
						char cur_seq_num = cur_seq_nums[left_kmer_posn+seq_kmer_i];
						if(strand == 1)
						{
							cur_seq_num = 3 - cur_seq_num;
						}

						// Compare the kmer nuclotide with the nucleotide from the sequence.
						if(nuc_2_num(valid_raw_enriched_motifs->at(i_kmer)->kmer[kmer_i]) != cur_seq_num)
						{
							cur_kmer_match = false;
							break;
						}
					} // i loop.

					// If this kmer matches, set it.
					if(cur_kmer_match)
					{
						t_k_mer_match_posn_info* new_match_posn = new t_k_mer_match_posn_info();
						new_match_posn->seq_i = i_seq;
						new_match_posn->left_kmer_posn = left_kmer_posn;
						new_match_posn->strand = strand;
						valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->push_back(new_match_posn);
					}
				} // left_kmer_posn loop.
			} // strand loop.
		} // i_seq loop.

		// For this kmer, compare the counts from the matches with the counts from tree.
		if(valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->size() != valid_raw_enriched_motifs->at(i_kmer)->count)
		{
			fprintf(stderr, "The number of matches does not match between the tree count and sequence search: %s: %d, %d\n", 
				valid_raw_enriched_motifs->at(i_kmer)->kmer,
				valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->size(), 
				valid_raw_enriched_motifs->at(i_kmer)->count);

			for(int i = 0; i < valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->size(); i++)
			{
				fprintf(stderr, "Seq %d (%d): %d (%d)\n", valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i)->seq_i, 
					numerized_pos_seq_lengths->at(valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i)->seq_i),
					valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i)->left_kmer_posn,
					(int)(valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i)->strand));
			}

			exit(0);
		}
	} // i_kmer loop for filling the sequence positions.

	// Go over all the kmers and compute the distance between the kmers, if their lengths are the same.
	int score_buffer[100][100];

	const double match_score = 1;
	const double mismatch_score = 0;
	const double gap_score = 0;
	FILE* f_per_kmer_n_matches = open_f("per_motif_match_counts.list", "w");
	FILE* f_motif_matches = open_f("motif_match_scores.list", "w");
	for(int i_kmer = 0; i_kmer < valid_raw_enriched_motifs->size(); i_kmer++)
	{
		if(i_kmer % 10000 == 0)
		{
			fprintf(stderr, "Processing %d. kmer                  \r", i_kmer);
		}

		int l_kmer1 = t_string::string_length(valid_raw_enriched_motifs->at(i_kmer)->kmer);

		int n_kmer_matches_per_cur_kmer = 0;

		// Build the pwm for the current k-mer:
		double** cur_kmer_pwm = new double*[l_kmer1*3+2];
		for(int i = 0; i < l_kmer1*3+2; i++)
		{
			cur_kmer_pwm[i] = new double[4];
			memset(cur_kmer_pwm[i], 0, sizeof(double)*4);
		} // i loop.

		char cur_kmer_vicinity_string[1000];
		memset(cur_kmer_vicinity_string, 0,1000);

		char cur_kmer_pwm_string_fp[1000];
		sprintf(cur_kmer_pwm_string_fp, "kmer_%d_pwm_string.list", i_kmer);
		FILE* f_cur_kmer_pwm_string = open_f(cur_kmer_pwm_string_fp, "w");	

		// Update the k-mer counts: Extract the sequence around the k-mer and add to the counts: Using the matches of this k-mer over the sequences, update the counts 
		// of 3*l_kmer1 long k-mers to build the pwm.
		for(int i_m = 0; i_m < valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->size(); i_m++)
		{
			int i_match_seq = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->seq_i;
			int match_strand = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->strand;
			int match_left_kmer_posn = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn;
			int* matching_seq_numseq = numerized_pos_sequences->at(i_match_seq);

			if(valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn < l_kmer1 ||
				valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn + 2 * l_kmer1 >= numerized_pos_seq_lengths->at(i_match_seq))
			{
				// This is not a valid case to process.
			}
			else
			{	
				int cur_match_start = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn - l_kmer1;
				int cur_match_end = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn + 2 * l_kmer1;

				// Based on the strand, update the counts in the pwm for the current k-mer.
				if(match_strand == 0)
				{
					int pwm_i = 0;
					for(int i = cur_match_start; i < cur_match_end; i++)
					{
						cur_kmer_pwm[pwm_i][matching_seq_numseq[i]]++;
						cur_kmer_vicinity_string[pwm_i] = num_2_nuc(matching_seq_numseq[i]);
						pwm_i++;
					} // i loop.
				}
				else
				{
					int pwm_i = 0;
					for(int i = cur_match_end-1; i >= cur_match_start; i--)
					{
						cur_kmer_pwm[pwm_i][3-matching_seq_numseq[i]]++;
						cur_kmer_vicinity_string[pwm_i] = num_2_nuc(3-matching_seq_numseq[i]);
						pwm_i++;
					} // i loop.
				}
			}

			fprintf(f_cur_kmer_pwm_string, "%s\n", cur_kmer_vicinity_string);
		} // i_m loop.
				
		// Go over all remaining kmers.
		for(int j_kmer = i_kmer+1; j_kmer < valid_raw_enriched_motifs->size(); j_kmer++)
		{
			//int i_match;
			if(strcmp(valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer) == 0)
			{
				fprintf(stderr, "The compared k-mers are the same: %s, %s\n",valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer);
				exit(0);
			}

			if(t_string::string_length(valid_raw_enriched_motifs->at(i_kmer)->kmer) != t_string::string_length(valid_raw_enriched_motifs->at(j_kmer)->kmer))
			{
				// If the k-mers have different lengths, do not compare them.
			}
			else
			{
				int l_kmer2 = t_string::string_length(valid_raw_enriched_motifs->at(j_kmer)->kmer);
				char* kmer1 = valid_raw_enriched_motifs->at(i_kmer)->kmer;
				char* kmer2 = valid_raw_enriched_motifs->at(j_kmer)->kmer;

				score_buffer[0][0] = 0;
				score_buffer[0][0] += (kmer1[0]==kmer2[0]);
				for(int i = 0; i < l_kmer1; i++)
				{
					score_buffer[i][0] = (kmer1[i]==kmer2[0]);
				} //i loop.

				for(int j = 0; j < l_kmer2; j++)
				{
					score_buffer[0][j] = (kmer1[0]==kmer2[j]);
				} // j loop

				// Go over all the values and find the max value.
				int cur_kmers_match_score = 0;
				int max_match_i = 0;
				int max_match_j = 0;
				for(int i = 1; i < l_kmer1; i++)
				{
					for(int j = 1; j < l_kmer2; j++)
					{
						score_buffer[i][j] = score_buffer[i-1][j-1] + (kmer1[i]==kmer2[j]);

						// Update the score.
						if(cur_kmers_match_score < score_buffer[i][j])
						{
							cur_kmers_match_score = score_buffer[i][j];
							max_match_i = i;
							max_match_j = j;
						}
					} // j loop.
				} // i loop.

				// Is this a significant match?
				int l_min_kmer = (l_kmer2>l_kmer1)?(l_kmer1):(l_kmer2);
			
				// If the match score is greater than or equal to the .75 of the length of the smaller k-mer, process this pair.
				const double ALIGNMENT_SCORE_CUTOFF = .75;
				if(cur_kmers_match_score >= l_min_kmer*ALIGNMENT_SCORE_CUTOFF)
				{
					fprintf(f_motif_matches, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n", 
						valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(i_kmer)->count, valid_raw_enriched_motifs->at(i_kmer)->bckgrnd_count,
						valid_raw_enriched_motifs->at(j_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->count, valid_raw_enriched_motifs->at(j_kmer)->bckgrnd_count, 
						cur_kmers_match_score);

					// Update the number of kmers that this kmer matches to.
					n_kmer_matches_per_cur_kmer++;

					// Update the k-mer PWM for this match:
					for(int i_m = 0; i_m < valid_raw_enriched_motifs->at(j_kmer)->k_mer_match_posns->size(); i_m++)
					{
						int i_match_seq = valid_raw_enriched_motifs->at(j_kmer)->k_mer_match_posns->at(i_m)->seq_i;
						int match_strand = valid_raw_enriched_motifs->at(j_kmer)->k_mer_match_posns->at(i_m)->strand;
						int match_left_kmer_posn = valid_raw_enriched_motifs->at(j_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn;
						int* matching_seq_numseq = numerized_pos_sequences->at(i_match_seq);
	
						//int cur_match_start = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn - l_kmer1;
						//int cur_match_end = valid_raw_enriched_motifs->at(i_kmer)->k_mer_match_posns->at(i_m)->left_kmer_posn + 2 * l_kmer1;
						int cur_match_end = match_left_kmer_posn + 2 * l_kmer1 - max_match_i + max_match_j;

						// The match end must be changed for negative strand because 
						if(match_strand == 1)
						{
							cur_match_end = match_left_kmer_posn + 2 * l_kmer1 - max_match_j + max_match_i;
						}

						int cur_match_start = cur_match_end - 3 * l_kmer1;

						if(cur_match_end > 0 && 
							cur_match_start > 0 &&
							cur_match_end < numerized_pos_seq_lengths->at(i_match_seq))
						{
							char cur_match_kmer_vicinity_string[1000];
							memset(cur_match_kmer_vicinity_string, 0,1000);

							// Based on the strand, update the counts in the pwm for the current k-mer.
							if(match_strand == 0)
							{
								int pwm_i = 0;
								for(int i = cur_match_start; i < cur_match_end; i++)
								{
									cur_kmer_pwm[pwm_i][matching_seq_numseq[i]]++;
									cur_match_kmer_vicinity_string[pwm_i] = num_2_nuc(matching_seq_numseq[i]);
									pwm_i++;
								} // i loop.
							}
							else
							{
								int pwm_i = 0;
								for(int i = cur_match_end-1; i >= cur_match_start; i--)
								{
									cur_kmer_pwm[pwm_i][3-matching_seq_numseq[i]]++;
									cur_match_kmer_vicinity_string[pwm_i] = num_2_nuc(3 - matching_seq_numseq[i]);
									pwm_i++;
								} // i loop.
							}

							fprintf(f_cur_kmer_pwm_string, "%s\n", cur_match_kmer_vicinity_string);

if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
							// Dump the sequences.
							fprintf(stderr, "Current matching sequences (%d, %d, %s, %s, %d):\n%s\n%s\n", max_match_i, max_match_j, valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer, match_strand, cur_kmer_vicinity_string, cur_match_kmer_vicinity_string);
							fprintf(stderr, "The sequence info: j_kmer: %d, matching_seq_i: %d, left_kmer_posn: %d\n", j_kmer, i_match_seq, match_left_kmer_posn);
							getc(stdin);
}
						} // The check for validity of the coordinates.
					} // i_m loop.

if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
					fprintf(stderr, "Pass: %s\t%s\t%d\n", valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer, cur_kmers_match_score);
				}
				else
				{
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
					fprintf(stderr, "Not Pass: %s\t%s\t%d\n", valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(j_kmer)->kmer, cur_kmers_match_score);
				}
				//getc(stdin);
			}
		} // j_kmer loop.

		fprintf(f_per_kmer_n_matches, "%s\t%d\t%d\n", 
			valid_raw_enriched_motifs->at(i_kmer)->kmer, valid_raw_enriched_motifs->at(i_kmer)->count, n_kmer_matches_per_cur_kmer);

		// Close the file with pwm strings.
		fclose(f_cur_kmer_pwm_string);

		// Dump the pwm for the current kmer.
		char cur_kmer_pwm_fp[1000];
		sprintf(cur_kmer_pwm_fp, "kmer_%d_%d_pwm.list", t_string::string_length(valid_raw_enriched_motifs->at(i_kmer)->kmer), i_kmer);
		FILE* f_cur_kmer_pwm = open_f(cur_kmer_pwm_fp, "w");
		//double** cur_kmer_pwm = new double*[l_kmer1*3+2];
		for(int i = 0; i < 3*l_kmer1; i++)
		{
			//cur_kmer_pwm[i] = new double[4];
			double cur_total = cur_kmer_pwm[i][0] + cur_kmer_pwm[i][1] + cur_kmer_pwm[i][2] + cur_kmer_pwm[i][3];
			fprintf(f_cur_kmer_pwm, "%lf\t%lf\t%lf\t%lf\n", cur_kmer_pwm[i][0]/cur_total, 
				cur_kmer_pwm[i][1]/cur_total,
				cur_kmer_pwm[i][2]/cur_total,
				cur_kmer_pwm[i][3]/cur_total);

			delete [] cur_kmer_pwm[i];
		} // i loop.
		delete [] cur_kmer_pwm;
		fclose(f_cur_kmer_pwm);

		//exit(0);
		getc(stdin);
	} // i_kmer loop
	fclose(f_motif_matches);
	fclose(f_per_kmer_n_matches);
	
	// Cluster the motifs with respect to the match scores: Use the match score as a similarity measure and do merging of the motifs.
}

void update_kmer_frequency_tree(t_kmer_count_tree_entry* main_entry, 
	int* new_seq, 
	int l_cur_seq, 
	int k)
{
	// Process both strands.
	for(int strand = 0; strand <= 1; strand++)
	{
		for(int win_pos_left_i = 0; win_pos_left_i <= (l_cur_seq - k); win_pos_left_i++)
		{
			// Build the kmer frequency tree from the list of positive sequences.
			int* cur_num_seq = &(new_seq[win_pos_left_i]);
			t_kmer_count_tree_entry* cur_entry = main_entry;

			//// Skip the kmers that have weird compositions.
			//if(!check_kmer_composition(positive_num_seqs->at(i_seq), win_pos_left_i, k))
			//{
			//	continue;
			//}

			// Go over all the nucleotides.
			for(int nuc_i = 0; nuc_i < k; nuc_i++)
			{
				// Move to the next entry: We are currently at the entry for the previous nucleotide (main_entry if we are at the first nuc), 
				// move to the entry for the current nucleotide. Note that the next nuc entries may not have been allocated for this k-mer, yet.
				if(cur_entry->next_nuc_entries == NULL)
				{
					cur_entry->next_nuc_entries = new t_kmer_count_tree_entry*[4];
					for(int i = 0; i < 4; i++)
					{
						cur_entry->next_nuc_entries[i] = new t_kmer_count_tree_entry();
						cur_entry->next_nuc_entries[i]->k = cur_entry->k + 1;
						cur_entry->next_nuc_entries[i]->count = 0;
						cur_entry->next_nuc_entries[i]->bckgrnd_count = 0;
						cur_entry->next_nuc_entries[i]->next_nuc_entries = NULL;
					} // i loop.
				}

				// Move to the next entry.
				cur_entry = cur_entry->next_nuc_entries[cur_num_seq[nuc_i]];

				// Update the count.
				cur_entry->count++;
			} // nuc_i loop.
		} // win_pos_left_i loop.

		// Reverse complement the string.
		for(int i_nuc = 0; i_nuc < l_cur_seq/2; i_nuc++)
		{
			int cur_nuc = 3 - new_seq[i_nuc];
			int cur_nuc_other_half = 3 - new_seq[l_cur_seq - i_nuc - 1];

			// Switch the nucleotides.
			new_seq[i_nuc] = cur_nuc_other_half;
			new_seq[l_cur_seq - i_nuc - 1] = cur_nuc;
		} // i_nuc loop.

		// If the sequence length is odd, complement the nucleotide in the middle.
		if(l_cur_seq % 2 == 1)
		{
			new_seq[l_cur_seq/2] = 3 - new_seq[l_cur_seq/2];
		}
	} // strand loop.
}

t_kmer_count_tree_entry* init_kmer_frequency_tree(int k)
{
	// Allocate an entry: Do not allocate all the next nucleotide entries. Those are going to be allocated 
	t_kmer_count_tree_entry* main_entry = new t_kmer_count_tree_entry();
	main_entry->k = 0;
	main_entry->count = 1;
	main_entry->bckgrnd_count = 1;
	main_entry->next_nuc_entries = NULL;

	return(main_entry);
}

t_kmer_count_tree_entry* build_kmer_frequency_tree(vector<int*>* positive_num_seqs, vector<int>* positive_num_seq_lengths, int k)
{
	// Allocate an entry: Do not allocate all the next nucleotide entries. Those are going to be allocated 
	t_kmer_count_tree_entry* main_entry = new t_kmer_count_tree_entry();
	main_entry->k = 0;
	main_entry->count = 1;
	main_entry->bckgrnd_count = 1;
	main_entry->next_nuc_entries = NULL;

	// Go over all the sequences.
	for(int i_seq = 0; i_seq < positive_num_seqs->size(); i_seq++)
	{
		int l_cur_seq = positive_num_seq_lengths->at(i_seq);

		// Process both strands.
		for(int strand = 0; strand <= 1; strand++)
		{
			for(int win_pos_left_i = 0; win_pos_left_i <= (l_cur_seq - k); win_pos_left_i++)
			{
				// Build the kmer frequency tree from the list of positive sequences.
				int* cur_num_seq = &(positive_num_seqs->at(i_seq)[win_pos_left_i]);
				t_kmer_count_tree_entry* cur_entry = main_entry;

				//// Skip the kmers that have weird compositions.
				//if(!check_kmer_composition(positive_num_seqs->at(i_seq), win_pos_left_i, k))
				//{
				//	continue;
				//}

				// Go over all the nucleotides.
				for(int nuc_i = 0; nuc_i < k; nuc_i++)
				{
					// Move to the next entry: We are currently at the entry for the previous nucleotide (main_entry if we are at the first nuc), 
					// move to the entry for the current nucleotide. Note that the next nuc entries may not have been allocated for this k-mer, yet.
					if(cur_entry->next_nuc_entries == NULL)
					{
						cur_entry->next_nuc_entries = new t_kmer_count_tree_entry*[4];
						for(int i = 0; i < 4; i++)
						{
							cur_entry->next_nuc_entries[i] = new t_kmer_count_tree_entry();
							cur_entry->next_nuc_entries[i]->k = cur_entry->k + 1;
							cur_entry->next_nuc_entries[i]->count = 0;
							cur_entry->next_nuc_entries[i]->bckgrnd_count = 0;
							cur_entry->next_nuc_entries[i]->next_nuc_entries = NULL;
						} // i loop.
					}

					// Move to the next entry.
					cur_entry = cur_entry->next_nuc_entries[cur_num_seq[nuc_i]];

					// Update the count.
					cur_entry->count++;
				} // nuc_i loop.
			} // win_pos_left_i loop.

			// Reverse complement the string.
			for(int i_nuc = 0; i_nuc < l_cur_seq/2; i_nuc++)
			{
				int cur_nuc = 3 - positive_num_seqs->at(i_seq)[i_nuc];
				int cur_nuc_other_half = 3 - positive_num_seqs->at(i_seq)[l_cur_seq - i_nuc - 1];

				// Switch the nucleotides.
				positive_num_seqs->at(i_seq)[i_nuc] = cur_nuc_other_half;
				positive_num_seqs->at(i_seq)[l_cur_seq - i_nuc - 1] = cur_nuc;
			} // i_nuc loop.

			// If the sequence length is odd, complement the nucleotide in the middle.
			if(l_cur_seq % 2 == 1)
			{
				positive_num_seqs->at(i_seq)[l_cur_seq/2] = 3 - positive_num_seqs->at(i_seq)[l_cur_seq/2];
			}
		} // strand loop.
	} // i_seq loop.

	return(main_entry);
}

// Recursively dump the kmers and the frequencies.
void dump_seqs_freqs(FILE* f_op, t_kmer_count_tree_entry* cur_entry, char* cur_kmer)
{
	if(cur_entry->next_nuc_entries == NULL)
	{
		// End the string.
		cur_kmer[cur_entry->k] = 0;

		if(cur_entry->count > 0)
		{
			fprintf(f_op, "%s\t%d\n", cur_kmer, cur_entry->count);
		}
	}
	else
	{
		for(int i = 0; i < 4; i++)
		{	
			if(cur_entry->next_nuc_entries[i] != NULL)
			{
				cur_kmer[cur_entry->k] = num_2_nuc(i);
				dump_seqs_freqs(f_op, cur_entry->next_nuc_entries[i], cur_kmer);
			}
		} // i loop.
	}
}

// Recursively count the total number of k-mers in the tree.
void delete_kmer_freq_tree(t_kmer_count_tree_entry* cur_tree_entry_2_delete)
{
	// Are there next entries next to the current node?
	if(cur_tree_entry_2_delete->next_nuc_entries != NULL)
	{
		int total_counts = 0;
		for(int i = 0; i < 4; i++)
		{
			delete_kmer_freq_tree(cur_tree_entry_2_delete->next_nuc_entries[i]);
		} // i loop.

		delete [] cur_tree_entry_2_delete->next_nuc_entries;
	}
	
	delete cur_tree_entry_2_delete;
}

// Recursively count the total number of k-mers in the tree.
int get_total_kmer_count(t_kmer_count_tree_entry* cur_tree_entry, int k_2_count)
{
	if(cur_tree_entry->k == k_2_count)
	{
		return(cur_tree_entry->count);
	}
	else if(cur_tree_entry->k < k_2_count)
	{
		// Are there next entries next to the current node?
		if(cur_tree_entry->next_nuc_entries != NULL)
		{
			int total_counts = 0;
			for(int i = 0; i < 4; i++)
			{
				total_counts += (get_total_kmer_count(cur_tree_entry->next_nuc_entries[i], k_2_count));
			}

			return(total_counts);
		}
		else
		{
			return(0);
		}
	}
	else
	{
		// We are past the k we want to process, return 0 without further recursing.
		return(0);
	}
}

// Recursively set the negative counts for kmers from the negative kmer frequency tree.
void get_bckgrnd_counts(t_kmer_count_tree_entry* cur_tree_entry, t_kmer_count_tree_entry* cur_neg_tree_entry)
{
	if(cur_tree_entry->next_nuc_entries == NULL)
	{
		// This is a leaf node, set the background count.
		cur_tree_entry->bckgrnd_count = cur_neg_tree_entry->count;
	}
	else
	{
		if(cur_tree_entry->next_nuc_entries != NULL &&
			cur_neg_tree_entry->next_nuc_entries != NULL)
		{
		}
		else
		{
			return;
		}

		// Recurse to all the branches.
		for(int i = 0; i < 4; i++)
		{	
			// Make sure that both entries have a branch here.
			if(cur_tree_entry->next_nuc_entries[i] != NULL &&
				cur_neg_tree_entry->next_nuc_entries[i] != NULL)
			{
				get_bckgrnd_counts(cur_tree_entry->next_nuc_entries[i], cur_neg_tree_entry->next_nuc_entries[i]);
			}
		} // i loop.
	}
}

// Recursively dump the kmers and the frequencies.
void get_kmer_count_array(vector<t_kmer_freq_info*>* cur_kmer_list, t_kmer_count_tree_entry* cur_entry, char* cur_kmer)
{
	if(cur_entry->next_nuc_entries == NULL)
	{
		// End the string.
		cur_kmer[cur_entry->k] = 0;

		if(cur_entry->count > 0)
		{
			//fprintf(stderr, "%s\t%d\n", cur_kmer, cur_entry->count);
			t_kmer_freq_info* new_kmer_freq_info = new t_kmer_freq_info();
			new_kmer_freq_info->count = cur_entry->count;
			new_kmer_freq_info->bckgrnd_count = cur_entry->bckgrnd_count;
			new_kmer_freq_info->kmer = t_string::copy_me_str(cur_kmer);
			cur_kmer_list->push_back(new_kmer_freq_info);
		}
	}
	else
	{
		for(int i = 0; i < 4; i++)
		{	
			if(cur_entry->next_nuc_entries[i] != NULL)
			{
				cur_kmer[cur_entry->k] = num_2_nuc(i);
				get_kmer_count_array(cur_kmer_list, cur_entry->next_nuc_entries[i], cur_kmer);
			}
		} // i loop.
	}
}

// Must add support for precomputing the positions that have weird compositions so as to decrease run time.
// Must add support for known weird motifs: Tandem repeats, etc, i.e., AGAGAG, CGCGCG...
bool check_kmer_composition(int* num_seq, int left_win_pos, int k)
{
	int nuc_nums[8];
	memset(nuc_nums, 0, sizeof(int) * 8);
	int consecutive_nuc_num = 0;

	int max_valid_nuc_count = k / 2;

	// If there is a consecutive run of longer than 2 nucleotides, return.
	//int max_valid_consecutive_nuc_count = k / 3;

	// Dinucleotide frequencies of the nucleotides.
	int homonuc_cnts = 1;
	int dinuc_freqs[20];
	memset(dinuc_freqs, 0, sizeof(int) * 20);
	for(int i_nuc = left_win_pos; i_nuc < (left_win_pos+k); i_nuc++)
	{
		// Dinucleotide frequency check on the validity.
		if(i_nuc > left_win_pos)
		{
			dinuc_freqs[num_seq[i_nuc]*4 + num_seq[i_nuc-1]]++;

			// If the more than 60% of the nucleotides are in dinuc repeats, exclude the k-mer.
			if((dinuc_freqs[num_seq[i_nuc]*4 + num_seq[i_nuc-1]]*2) >= k*5/10)
			{
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
				fprintf(stderr, "Skipping wrt dinuc frequency criteria:\n");
				for(int i = left_win_pos; i < left_win_pos+k; i++)
				{
					fprintf(stderr, "%c", num_2_nuc(num_seq[i]));
				} // i loop.
				fprintf(stderr, "\n");
				getc(stdin);
}
				return(false);
			}
		}

		// Do homonucleotide count check.
		if(i_nuc > left_win_pos)
		{
			if(num_seq[i_nuc] == num_seq[i_nuc - 1])
			{
				homonuc_cnts++;

				// If more than half of the nucleotides are in a homonucleotide sequence, exclude the k-mer.
				if(homonuc_cnts >= k/2)
				{
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
					fprintf(stderr, "Skipping wrt homonucleotide frequency criteria:\n");
					for(int i = left_win_pos; i < left_win_pos+k; i++)
					{
						fprintf(stderr, "%c", num_2_nuc(num_seq[i]));
					} // i loop.
					fprintf(stderr, "\n");
					getc(stdin);
}

					return(false);
				}
			}
			else
			{
				homonuc_cnts = 1;
			}
		}
		else
		{
			homonuc_cnts = 1;
		}
		
		nuc_nums[num_seq[i_nuc]]++;

		// Did one of the nucleotides become too frequent in the window?
		if(nuc_nums[num_seq[i_nuc]] > max_valid_nuc_count)
		{
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
			fprintf(stderr, "Skipping wrt nucleotide number criteria:\n");
			for(int i = left_win_pos; i < left_win_pos+k; i++)
			{
				fprintf(stderr, "%c", num_2_nuc(num_seq[i]));
			} // i loop.
			fprintf(stderr, "\n");
			getc(stdin);
}
			return(false);
		}
	} // i_nuc loop.

	// Go over the dinucleotide frequencies and break if there are too many dinucleotide's repeated.

	return(true);
}

void update_motif_counts_per_sequence_list(vector<int*>* numerized_sequences,
	vector<int>* numerized_seq_lengths,
	vector<int*>* per_nuc_depths,
	int min_k, int max_k, 
	double*** kmer_emission_log_probs,
	double*** kmer_observation_log_counts,
	bool positive_seq)
{
	// Go over all the sequences and align the profiles.
	for(int i_seq = 0; i_seq < (int)(numerized_sequences->size()); i_seq++)
	{
		if((i_seq % 1000) == 0)
		{
			fprintf(stderr, "Processing %d. sequence.                \r", i_seq);
		}

		// Scan the positive strand, then scan the reverse complement.

		int l_seq = numerized_seq_lengths->at(i_seq);

		// Loop over k values.
		for(int k = min_k; k <= max_k; k++)
		{
			// Process both strand for the current sequence.
			// Align the current sequence to the current k-mer profile: Move the current profile over the sequence and find the 
			int optimizing_scoring_profile_posn = l_seq;
			double optimal_log_score = (0.0);
			// If this set is the negative sequences, we will find the position that minimizes the score.
			if(!positive_seq)
			{
				optimal_log_score = 1000*1000;
			}

			int optimal_strand = 0;
			for(int strand = 0; strand <= 1; strand++)
			{
				// Slide over all the positions in this sequence.
				for(int left_seq_posn = 0; left_seq_posn < l_seq-k; left_seq_posn++)
				{
					// First check if this position has a weird composition that should be excluded from motif identification: Low complexity regions.
					if(!check_kmer_composition(numerized_sequences->at(i_seq), left_seq_posn, k))
					{
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
						fprintf(stderr, "Skipping:\n");
						for(int i_nuc = left_seq_posn; i_nuc < (left_seq_posn+k); i_nuc++)
						{
							fprintf(stderr, "%c", num_2_nuc(numerized_sequences->at(i_seq)[i_nuc]));
						} // i_nuc loop.
						fprintf(stderr, "\n");
						getc(stdin);
}
						continue;
					}

					//double cur_profile_posn_log_score = xlog(1.0);
					double cur_profile_posn_log_score = (1.0);

					// Compute the current matching score.
					if(strand == 0)
					{
						for(int left_nuc_posn = 0; left_nuc_posn < k; left_nuc_posn++)
						{
							int profile_posn = left_nuc_posn;
							int cur_nuc_num = numerized_sequences->at(i_seq)[left_seq_posn+left_nuc_posn];

							//int cur_nuc_num = nuc_2_num(sequences->at(i_seq)[profile_posn+posn]);
							//int cur_nuc_num = numerized_sequences->at(i_seq)[left_seq_posn+posn];
							//cur_profile_posn_log_score = xlog_mul(cur_profile_posn_log_score, kmer_emission_log_probs[k][posn][cur_nuc_num]);
					
							cur_profile_posn_log_score = (cur_profile_posn_log_score * kmer_emission_log_probs[k][profile_posn][cur_nuc_num]);
						} // posn loop.
					}
					else
					{
						// Go in the negative direction.
						for(int left_nuc_posn = 0; left_nuc_posn < k; left_nuc_posn++)
						{
							int profile_posn = k - left_nuc_posn - 1;
							int cur_nuc_num = 3 - numerized_sequences->at(i_seq)[left_seq_posn+left_nuc_posn];

							//int cur_nuc_num = nuc_2_num(sequences->at(i_seq)[profile_posn+posn]);
							//int cur_nuc_num = numerized_sequences->at(i_seq)[left_seq_posn+posn];
							//cur_profile_posn_log_score = xlog_mul(cur_profile_posn_log_score, kmer_emission_log_probs[k][posn][cur_nuc_num]);
					
							cur_profile_posn_log_score = (cur_profile_posn_log_score * kmer_emission_log_probs[k][profile_posn][cur_nuc_num]);
						} // posn loop.
					}

if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
					fprintf(stderr, "profile_posn: %d, log_score: %lf\n", left_seq_posn, cur_profile_posn_log_score);

					// If this is the positive sequence list, find the maximizing position.
					if(positive_seq)
					{
						if(cur_profile_posn_log_score > optimal_log_score)
						{
							optimizing_scoring_profile_posn = left_seq_posn;
							optimal_log_score = cur_profile_posn_log_score;
							optimal_strand = strand;
						}
					}
					else
					{
						// If this is the negative sequence list, find the minimizing position.
						if(cur_profile_posn_log_score < optimal_log_score)
						{
							optimizing_scoring_profile_posn = left_seq_posn;
							optimal_log_score = cur_profile_posn_log_score;
							optimal_strand = strand;
						}
					}
				} // left_seq_posn loop.
			} // strand loop.

if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
{
			fprintf(stderr, "k: %d; Highest scoring profile position: %d, Highest score: %lf, Highest scoring strand: %d\n", k, optimizing_scoring_profile_posn, optimal_log_score, optimal_strand);
			for(int left_nuc_posn = 0; left_nuc_posn < k; left_nuc_posn++)
			{
				int profile_posn = left_nuc_posn;

				char cur_nuc = num_2_nuc(numerized_sequences->at(i_seq)[optimizing_scoring_profile_posn+left_nuc_posn]);
				fprintf(stderr, "%c", cur_nuc);
			} // left_nuc_posn loop.
			getc(stdin);
			fprintf(stderr, "\n");
}

			// Based on the highest scoring profile posn, and strand, update the counts. There may be a filter here to filter based on the scores.
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
			fprintf(stderr, "Matching sequence:\n");
			if(optimal_strand == 0)
			{
				for(int left_nuc_posn = 0; left_nuc_posn < k; left_nuc_posn++)
				{
					int profile_posn = left_nuc_posn;

					//int cur_nuc_num = nuc_2_num(sequences->at(i_seq)[highest_scoring_profile_posn+posn]);
					int cur_nuc_num = numerized_sequences->at(i_seq)[optimizing_scoring_profile_posn+left_nuc_posn];

//if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
//					fprintf(stderr, "%c", sequences->at(i_seq)[highest_scoring_profile_posn+posn]);
					
					// Add the score of the current alignment to the counts.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = xlog_sum(kmer_observation_log_counts[k][posn][cur_nuc_num], highest_log_score);

					// Add one to the counts.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = xlog_sum(kmer_observation_log_counts[k][posn][cur_nuc_num], 0.0);				

					// Update the total count using the read depth as a bias.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = (kmer_observation_log_counts[k][posn][cur_nuc_num] + 1.0);
					kmer_observation_log_counts[k][profile_posn][cur_nuc_num] = (kmer_observation_log_counts[k][profile_posn][cur_nuc_num] + per_nuc_depths->at(i_seq)[optimizing_scoring_profile_posn+left_nuc_posn]);
				} // posn loop.
			}
			else
			{
				for(int left_nuc_posn = 0; left_nuc_posn < k; left_nuc_posn++)
				{
					// Position on the profile.
					int profile_posn = k - left_nuc_posn - 1;

					//int cur_nuc_num = nuc_2_num(sequences->at(i_seq)[highest_scoring_profile_posn+posn]);
					int cur_nuc_num = 3 - numerized_sequences->at(i_seq)[optimizing_scoring_profile_posn+left_nuc_posn];

//if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
//					fprintf(stderr, "%c", sequences->at(i_seq)[highest_scoring_profile_posn+posn]);
					
					// Add the score of the current alignment to the counts.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = xlog_sum(kmer_observation_log_counts[k][posn][cur_nuc_num], highest_log_score);

					// Add one to the counts.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = xlog_sum(kmer_observation_log_counts[k][posn][cur_nuc_num], 0.0);				

					// Update the total count using the read depth as a bias.
					//kmer_observation_log_counts[k][posn][cur_nuc_num] = (kmer_observation_log_counts[k][posn][cur_nuc_num] + 1.0);
					kmer_observation_log_counts[k][profile_posn][cur_nuc_num] = (kmer_observation_log_counts[k][profile_posn][cur_nuc_num] + per_nuc_depths->at(i_seq)[optimizing_scoring_profile_posn+left_nuc_posn]);
				} // posn loop.
			}
if(_DUMP_GENOME_SEQUENCE_MESSAGES_)
			fprintf(stderr, "\n");
		} // k loop.
	} // i_seq loop.
}

void dump_mutated_embedded_sequences_w_multiple_kmers(vector<char*>* kmers, vector<double>* kmer_fractions, int l_seqs, int n_seqs, double mutation_prob, char* op_fp)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	if(kmer_fractions->size() != kmers->size())
	{
		fprintf(stderr, "The fractions must be of the same length as the kmers: %d, %d\n", kmer_fractions->size(), kmers->size());
		exit(0);
	}

	// Compute the number of sequences.
	vector<int>* nseqs_per_kmers = new vector<int>();
	for(int k_i = 0; k_i < (int)kmers->size(); k_i++)
	{	
		nseqs_per_kmers->push_back((int)(n_seqs * kmer_fractions->at(k_i)));
		fprintf(stderr, "%s: %d sequences.\n", kmers->at(k_i), nseqs_per_kmers->at(k_i));
	} // k loop.

	// Generate the sequences for each kmer for the number of sequences for the corresponding kmer (with mutations.)
	vector<char*>* sequences = new vector<char*>();

	// For each kmer, dump the sequences.
	for(int k_i = 0; k_i < kmers->size(); k_i++)
	{
		fprintf(stderr, "Generating sequences for %s\n", kmers->at(k_i));

		char* cur_kmer = kmers->at(k_i);

		int l_flank = l_seqs - t_string::string_length(kmers->at(k_i));

		// Generate sequences for the current kmer.
		for(int i_seq = 0; i_seq < nseqs_per_kmers->at(k_i); i_seq++)
		{
			char* cur_seq = new char[l_seqs+3];
			memset(cur_seq, 0, l_seqs+2);
			int i_nuc = 0;

			fprintf(stderr, "Dumping %d. sequence               \r", i_seq);

			// Select the position of the k-mer randomly.
			int cur_rand_kmer_posn = (int)(floor(rng->random_double_ran3() * l_flank));
			if(cur_rand_kmer_posn == l_flank)
			{
				cur_rand_kmer_posn--;
			}

			// Dump the first set of random nucleotides.
			for(int i = 0; i < cur_rand_kmer_posn; i++)
			{
				// Dump.
				char nucs[] = "ACGU";
				char cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
				if(cur_nuc_i == 4)
				{
					cur_nuc_i--;
				}

				//fprintf(f_op, "%c", nucs[cur_nuc_i]);
				cur_seq[i_nuc] = nucs[cur_nuc_i];
				i_nuc++;
			} // i loop.

			// Dump the k-mer: Add mutations where necessary.
			for(int i = 0; i < t_string::string_length(kmers->at(k_i)); i++)
			{
				char cur_nuc = 0;
				if(rng->random_double_ran3() < mutation_prob)
				{
					char nucs[] = "ACGU";
					int cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
					if(cur_nuc_i == 4)
					{
						cur_nuc_i--;
					}

					cur_nuc = nucs[cur_nuc_i];
				}
				else
				{
					cur_nuc = kmers->at(k_i)[i];
				}
				cur_seq[i_nuc] = cur_nuc;
				i_nuc++;
				//fprintf(f_op, "%c", cur_nuc);
			}

			// Dump the remaining random nucleotides.
			for(int i = cur_rand_kmer_posn+t_string::string_length(kmers->at(k_i)); i < l_seqs; i++)
			{
				// Dump.
				char nucs[] = "ACGU";
				char cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
				if(cur_nuc_i == 4)
				{
					cur_nuc_i--;
				}

				//fprintf(f_op, "%c", nucs[cur_nuc_i]);
				cur_seq[i_nuc] = nucs[cur_nuc_i];
				i_nuc++;
			} // i loop.

			//fprintf(f_op, "\n");

			// Add the current sequence to the list of sequences.
			sequences->push_back(cur_seq);
		} // i_seq loop.
	} // k_i loop.

	fprintf(stderr, "Generated %d sequences.\n", sequences->size());

	// Permute the sequences.
	vector<int>* permuted_indices = rng->permute_indices(sequences->size(), sequences->size());

	// Dump the sequences.
	FILE* f_op = open_f(op_fp, "w");
	for(int i_seq = 0; i_seq < sequences->size(); i_seq++)
	{
		fprintf(f_op, "%s\n", sequences->at(permuted_indices->at(i_seq)));
	} // i_seq loop.
	fclose(f_op);
}

void dump_embedded_sequences_w_kmer(char* kmer, int l_seqs, int n_seqs, double mutation_prob, char* op_fp)
{
	t_rng* rng = new t_rng(t_seed_manager::seed_me());

	int l_flank = l_seqs - t_string::string_length(kmer);

	FILE* f_op = open_f(op_fp, "w");

	for(int i_seq = 0; i_seq < n_seqs; i_seq++)
	{
		fprintf(stderr, "Dumping %d. sequence               \r", i_seq);

		// Select the position of the k-mer randomly.
		int cur_rand_kmer_posn = (int)(floor(rng->random_double_ran3() * l_flank));
		if(cur_rand_kmer_posn == l_flank)
		{
			cur_rand_kmer_posn--;
		}

		// Dump the first set of random nucleotides.
		for(int i = 0; i < cur_rand_kmer_posn; i++)
		{
			// Dump.
			char nucs[] = "ACGU";
			char cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
			if(cur_nuc_i == 4)
			{
				cur_nuc_i--;
			}

			fprintf(f_op, "%c", nucs[cur_nuc_i]);
		} // i loop.

		// Dump the k-mer: Add mutations where necessary.
		for(int i = 0; i < t_string::string_length(kmer); i++)
		{
			char cur_nuc = 0;
			if(rng->random_double_ran3() < mutation_prob)
			{
				char nucs[] = "ACGU";
				int cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
				if(cur_nuc_i == 4)
				{
					cur_nuc_i--;
				}

				cur_nuc = nucs[cur_nuc_i];
			}
			else
			{
				cur_nuc = kmer[i];
			}
			fprintf(f_op, "%c", cur_nuc);
		}

		// Dump the remaining random nucleotides.
		for(int i = cur_rand_kmer_posn+t_string::string_length(kmer); i < l_seqs; i++)
		{
			// Dump.
			char nucs[] = "ACGU";
			char cur_nuc_i = (int)(floor(rng->random_double_ran3() * 4));
			if(cur_nuc_i == 4)
			{
				cur_nuc_i--;
			}

			fprintf(f_op, "%c", nucs[cur_nuc_i]);
		} // i loop.

		fprintf(f_op, "\n");
	} // i_seq loop.

	fprintf(stderr, "\n");

	fclose(f_op);
}
