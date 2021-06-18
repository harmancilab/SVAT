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

bool __DUMP_SECURE_CRYPTAGGREGATE_INDEL_MSGS__ = false;

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Following is for encryption of the SNV genotype matrix. Please feel free to change this interface.
void encrypt_encoded_SNV_genotypes_2_matrix(char* encoded_vectorized_snvs_dir,
											char* encrypted_snv_genotype_matrix_dir,
											char* EOI_regs_BED_fp,
											char* sample_ids_list_fp, // We need this only to know the number of inidividuals in the dataset, nothing else.
											int genotype_encoding_type)
{
	fprintf(stderr, "Encrypting the encoded SNV genotypes.\n");

	// EOI stands for elements of interest. These are the regions that we use to focus on. These are the regions that are interesting and much much smaller than the whole genome.
	vector<t_annot_region*>* EOI_regs = load_BED(EOI_regs_BED_fp);
	fprintf(stderr, "Loaded %d EOI regions.\n", EOI_regs->size());

	// Note that EOIs are public knowledge and they are not meant to be encrypted. The vector coordinates are also not sensitive. They can be plainly stored.

	// We will process each chromosome separately. Get the chromosome id's.
	vector<char*>* chr_ids = get_chr_ids(EOI_regs);

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
		// Get the regions on this chromosome.
		vector<t_annot_region*>* cur_chr_EOI_regs = get_regions_per_chromosome(EOI_regs, chr_ids->at(i_chr));
		fprintf(stderr, "%d EOI regions on %s\n", cur_chr_EOI_regs->size(), chr_ids->at(i_chr));

		// Following is a special sorting that is used throughout to set the "vectorized coordinate system". (Figure 1c -- slide 4).
		// It first sorts with respect to start position, then wrt name of the element.
		sort_set_sorting_info(cur_chr_EOI_regs, sort_regions_per_start_end_name);

		// Setup the target regions on this chromosome.
		// This is necessary to setup the mapping between "genomic coordinates" and the "vectorized coordinates".
		int cur_chr_EOI_covg = (int)(coverage(cur_chr_EOI_regs));
		int cumul_signal_index = 1;
		int* per_track_posn_genomic_posn = new int[cur_chr_EOI_covg + 1];
		for (int i_reg = 0; i_reg < (int)cur_chr_EOI_regs->size(); i_reg++)
		{
			cur_chr_EOI_regs->at(i_reg)->score = cumul_signal_index;

			// Following sets mapping between vector-coordinates to the genomic-coordinates. Note that I use "track" and "vector" interchangeably.
			for (int cur_track_posn = cumul_signal_index;
				cur_track_posn < cumul_signal_index + cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start;
				cur_track_posn++)
			{
				// This array holds the genomic position for a given vector (or track) position 
				per_track_posn_genomic_posn[cur_track_posn] = cur_chr_EOI_regs->at(i_reg)->start + (cur_track_posn - cumul_signal_index);
			} // cur_track_posn loop.

			cumul_signal_index += (cur_chr_EOI_regs->at(i_reg)->end - cur_chr_EOI_regs->at(i_reg)->start + 1);
		} // i_reg loop.

		// Sanity check on the coverage; we should not have passed the total coverage of elements-of-interest.
		if (cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score != cur_chr_EOI_covg)
		{
			fprintf(stderr, "Score does not meet up coverage: %d, %d\n",
				cur_chr_EOI_regs->back()->end - cur_chr_EOI_regs->back()->start + cur_chr_EOI_regs->back()->score,
				cur_chr_EOI_covg);

			exit(0);
		}

		// Set the encoded genotype matrix name for this chromosome. This is the plaintext encoded genotype matrix.
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

				if (__DUMP_SECURE_CRYPTAGGREGATE_INDEL_MSGS__)
				{
					fprintf(stderr, "Read the encoded genotypes @ %d\n", cur_track_posn);
				}

				per_allele_encoded_posn_genotype_data[allele_i][cur_track_posn] = cur_posn_snv_genotypes;
			}
		} // file reading loop.
		close_f(f_geno_matrix, cur_chr_encoded_geno_matrix_fp);

		// We have loaded the genotype matrix. We need to encrypt it and save it below.

		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// MK TO FILL HERE
		// We can encrypt the genotype matrix position-by-position where position is plaintext. This is because 
		// the vector coordinates are standard and are not sensitive. They do not reveal any information since elements are public knowledge.
		// Only the genotype values are sensitive. 
		//////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// Feel free to change following.
		char cur_chr_encrypted_matrix_fp[1000];
		sprintf(cur_chr_encrypted_matrix_fp, "%s/%s_encrypted.bin.gz", encrypted_snv_genotype_matrix_dir, chr_ids->at(i_chr));

		// ......
		// ......
	} // i_chr loop.
}

// I am writing this function declaration for MK to complete.
void re_encrypt_genotype_matrix(char* encrypted_genotype_matrix_directory_list,
								char* reencrypted_snv_genotype_matrix_dir,
								char* EOI_regs_BED_fp,
								char* sample_ids_list_fp,
								int genotype_encoding_type,
								char* public_keys_file_path)
{

}

// Following is the plain version that we need to convert into HE.
unsigned int* secure_aggregate_encrypted_SNV_genotype_matrix(char* encrypted_vectorized_snv_genotypes_dir,
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
	fprintf(stderr, "Loaded %d samples.\n", (int)sample_ids->size());

	// This is the # of bytes for each variant that necessary to store the genotypes of all individuals; We need this to allocate/load the genotype matrix later.
	int n_bytes_per_var = (int)(ceil(((double)n_bits_per_sample * (sample_ids->size() + 1)) / 8));
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
		sprintf(cur_chr_encoded_geno_matrix_fp, "%s/%s_matrix.bin.gz", encrypted_vectorized_snv_genotypes_dir, chr_ids->at(i_chr));
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

				if (__DUMP_SECURE_CRYPTAGGREGATE_INDEL_MSGS__)
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

			fprintf(stderr, "Found genotypes for SNV @ %d, counting the alleles over %d samples\n", cur_SNV_posn, (int)sample_ids->size());

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
					if (__DUMP_SECURE_CRYPTAGGREGATE_INDEL_MSGS__)
					{
						fprintf(stderr, "Sample %d: %d\n", sample_i, cur_sample_snv_geno);
					}

					int allele_count = cur_sample_snv_geno;
					cur_var_count += allele_count;
				}
			} // sample_i loop.

			fprintf(stderr, "%s:%d (%s) -- %d/%d over positions.\n",
				cur_chr_VOI_regs->at(i_voi)->chrom, cur_SNV_posn, cur_chr_VOI_regs->at(i_voi)->name,
				cur_var_count, (int)sample_ids->size());
		} // i_voi loop.
	} // i_chr loop.

	return(NULL);
}
