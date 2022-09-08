#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svat_utils.h"
#include "svat_aggregation_utils.h"
#include "svat_genome_sequence_tools.h"
#include "svat_ansi_string.h"
#include "file_utils.h"

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
Variant Simulation:\n\
	-generate_simulated_InDels_on_EOI\n\
	-generate_simulated_SNVs_on_EOI\n\
VEP Track Generation:\n\
	-get_EOIs_BED_per_exonic_regions\n\
	-get_frame_gene_unique_EOIs_BED_per_CDS_regions\n\
	-generate_EOI_enumerating_SNV_VCF\n\
	-generate_EOI_enumerating_Indel_VCF\n\
Vectorization:\n\
	VEP Annotated Variants:\n\
		-separate_VEP_output_per_chromosome\n\
		-signalize_VEP_annotated_SNVs\n\
		-signalize_VEP_annotated_Inserts\n\
		-signalize_VEP_annotated_Deletes\n\
		-signalize_VEP_annotated_Deletes_w_CDS_vicinity\n\
	VCF Variants:\n\
		-separate_VCF_per_chromosome\n\
		-signalize_non_coding_variants_per_VCF\n\
		-signalize_SNV_genotypes_per_VCF\n\
		-signalize_Insertion_genotypes_per_VCF\n\
		-signalize_Deletion_genotypes_per_VCF\n\
Annotation Signal Processing:\n\
	Plaintext:\n\
		-multiply_signalized_SNVs_and_annotations\n\
		-multiply_signalized_Indels_and_annotations\n\
Variant Genotype Matrix Aggregation:\n\
	Plaintext Variant Encoding:\n\
		-encode_SNV_genotypes_2_matrix_per_VCF\n\
		-encode_DEL_genotypes_2_matrix_per_VCF\n\
	Plaintext Aggregation:\n\
		-plain_aggregate_encoded_Del_genotype_matrix\n\
		-plain_aggregate_encoded_SNV_genotype_matrix\n\
Plain text Impact Analysis:\n\
	-compare_VEP_annotations_with_SVAT_annotations\n\
	-compare_ANNOVAR_vs_SVAT_variant_wise\n\
	-compare_ANNOVAR_annotations_with_SVAT_annotations\n\
	-extract_annotated_SNVs_from_signal\n\
	-translate_annotated_Deletions_from_annotated_signals\n\
	-translate_annotated_Insertions_from_annotated_signals\n\
	-translate_annotated_SNVs_from_annotated_signals\n\
----------------------------------------------------------------\n\n", argv[0]);

		exit(0);
	}
	
	if (t_string::compare_strings(argv[1], "-convert_Deletion_annotations_XCVATR_2_VEP"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [XCVATR annotated Deletion path] [XCVATR 2 VEP term mapping list] [Output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* XCVATR_deletion_annotation_fp = argv[2];
		char* XCVATR_2_VEP_mapping_list_fp = argv[3];
		char* vep_op_fp = argv[4];

		// Load the mapper.
		fprintf(stderr, "Loading XCVATR->VEP term mapping from %s\n", XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* term_mapper_lines = buffer_file(XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* xcvatr_terms = new vector<char*>();
		vector<char*>* vep_terms = new vector<char*>();
		for (int i_l = 0; i_l < term_mapper_lines->size(); i_l++)
		{
			fprintf(stderr, "Parsing: %s\n", term_mapper_lines->at(i_l));

			t_string_tokens* cur_toks = t_string::tokenize_by_chars(term_mapper_lines->at(i_l), "\t");
			if (cur_toks->size() < 2)
			{
				fprintf(stderr, "Could not correctly parse: %d tokens.\n", cur_toks->size());
				exit(0);
			}

			xcvatr_terms->push_back(t_string::copy_me_str(cur_toks->at(0)->str()));
			vep_terms->push_back(t_string::copy_me_str(cur_toks->at(1)->str()));

			fprintf(stderr, "%s => %s\n", xcvatr_terms->back(), vep_terms->back());
		} // i_l loop.

		// Following are the VEP column indices.
		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		// 1_11479372_INS_ENST00000294484.7        1:11479372-11479373     N       ENSG00000204624 ENST00000423056 Transcript      splice_region_variant,non_coding_transcript_variant     207-208 -
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		int cds_posn_counter = 30;

		//char cumul_var_name[1000];
		//cumul_var_name[0] = 0;
		//char cumul_var_impacting_transcript_id[1000];
		//cumul_var_impacting_transcript_id[0] = 0;
		//vector<char*>* cumul_var_vep_impacts = new vector<char*>();

		int line_i = 0;

		FILE* f_xcvtr_deletion_annots = open_f(XCVATR_deletion_annotation_fp, "r");
		FILE* f_op = open_f(vep_op_fp, "w");
		while (1)
		{
			char* cur_line = getline(f_xcvtr_deletion_annots);
			if (cur_line == NULL)
			{
				break;
			}

			if (cur_line[0] == '#')
			{
				delete[] cur_line;
				continue;
			}

			//1       65408   65409   A . 1_65408_DEL_ENST00000641515.2       0 + EFFECT_PROMOTOR_OVERLAP ENST00000641515.2
			t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");

			char* cur_var_chrom = cur_line_toks->at(0)->str();
			char* cur_var_posn = cur_line_toks->at(1)->str();
			char* cur_var_impact_str = cur_line_toks->at(6)->str();
			char* cur_var_imp_transcript_id = cur_line_toks->at(7)->str();
			char* cur_var_name = cur_line_toks->at(3)->str();
			t_string_tokens* name_toks = t_string::tokenize_by_chars(cur_var_name, " ");
			char* ref_all = name_toks->at(0)->str();
			char* alt_all = name_toks->at(1)->str();
			t_string_tokens* cur_var_id_toks = name_toks->at(2)->tokenize_by_chars("_");

			// Make sure this is an insert.
			if (!t_string::compare_strings(cur_var_id_toks->at(2)->str(), "DEL"))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_var_id_toks);
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			if (line_i % 100000 == 0)
			{
				fprintf(stderr, "@ %d. annotation: %s                    \r", line_i, cur_var_name);
			}

			line_i++;

			// Following converts from BED to VEP coordinates, it is the position of the first nucleotide that is deleted.
			int cur_del_start_pos = atoi(cur_var_posn) + 1;

			// Make sure what we are matching is what we generated.
			if (!t_string::compare_strings(cur_var_id_toks->at(3)->str(), cur_var_imp_transcript_id))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_var_id_toks);
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			t_string::clean_tokens(cur_var_id_toks);

			//fprintf(stderr, "Alleles: %c -> %c\n", ref_all[0], alt_all[0]);

			cds_posn_counter++;
			if (cds_posn_counter > 500)
			{
				cds_posn_counter = 30;
			}

			// Map the impact.
			int xcvtr_imp_str_i = t_string::get_i_str(xcvatr_terms, cur_var_impact_str);
			if (xcvtr_imp_str_i == xcvatr_terms->size())
			{
				fprintf(stderr, "Could not find the impact string: %s\n", cur_var_impact_str);
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			char* mapped_vep_impact_str = vep_terms->at(xcvtr_imp_str_i);

			char cur_var_id[1000];
			char cur_var_alt_all[1000];
			if (sscanf(cur_var_name, "%*s %s %s", cur_var_alt_all, cur_var_id) != 2)
			{
				fprintf(stderr, "Could not parse variant identifier from %s\n", cur_var_name);
				exit(0);
			}

			//fprintf(stderr, "PROCESSED_LINE: %s (%d tokens)\n", cur_line, cur_line_toks->size());
			// 1_11479369_DEL_ENST00000294484.7        1:11479370      -       ENSG00000204624 ENST00000294484 Transcript      splice_region_variant,5_prime_UTR_variant       216     -       -       -       -       -       IMPACT=LOW;STRAND=1
			char feat_type[] = "transcript";
			fprintf(f_op, "%s\t%s:%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n",
				cur_var_id, cur_var_chrom, cur_del_start_pos, "-",
				cur_var_imp_transcript_id, cur_var_imp_transcript_id, feat_type,
				mapped_vep_impact_str, cds_posn_counter, cds_posn_counter);

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			t_string::clean_tokens(name_toks);
			delete[] cur_line;
		} // xcvtr annotation reading loop.

		close_f(f_xcvtr_deletion_annots, XCVATR_deletion_annotation_fp);
		close_f(f_op, vep_op_fp);
	} // -convert_Deletion_annotations_XCVATR_2_VEP option.
	else if (t_string::compare_strings(argv[1], "-convert_Insertion_annotations_XCVATR_2_VEP"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [XCVATR annotated Insertion path] [XCVATR 2 VEP term mapping list] [Output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* XCVATR_insert_annotation_fp = argv[2];
		char* XCVATR_2_VEP_mapping_list_fp = argv[3];
		char* vep_op_fp = argv[4];

		// Load the mapper.
		fprintf(stderr, "Loading XCVATR->VEP term mapping from %s\n", XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* term_mapper_lines = buffer_file(XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* xcvatr_terms = new vector<char*>();
		vector<char*>* vep_terms = new vector<char*>();
		for (int i_l = 0; i_l < term_mapper_lines->size(); i_l++)
		{
			fprintf(stderr, "Parsing: %s\n", term_mapper_lines->at(i_l));

			t_string_tokens* cur_toks = t_string::tokenize_by_chars(term_mapper_lines->at(i_l), "\t");
			if (cur_toks->size() < 2)
			{
				fprintf(stderr, "Could not correctly parse: %d tokens.\n", cur_toks->size());
				exit(0);
			}

			xcvatr_terms->push_back(t_string::copy_me_str(cur_toks->at(0)->str()));
			vep_terms->push_back(t_string::copy_me_str(cur_toks->at(1)->str()));

			fprintf(stderr, "%s => %s\n", xcvatr_terms->back(), vep_terms->back());
		} // i_l loop.

		// Following are the VEP column indices.
		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		// 1_11479372_INS_ENST00000294484.7        1:11479372-11479373     N       ENSG00000204624 ENST00000423056 Transcript      splice_region_variant,non_coding_transcript_variant     207-208 -
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		int cds_posn_counter = 30;

		//char cumul_var_name[1000];
		//cumul_var_name[0] = 0;
		//char cumul_var_impacting_transcript_id[1000];
		//cumul_var_impacting_transcript_id[0] = 0;
		//vector<char*>* cumul_var_vep_impacts = new vector<char*>();

		int line_i = 0;

		FILE* f_xcvtr_insertion_annots = open_f(XCVATR_insert_annotation_fp, "r");
		FILE* f_op = open_f(vep_op_fp, "w");
		while (1)
		{
			char* cur_line = getline(f_xcvtr_insertion_annots);
			if (cur_line == NULL)
			{
				break;
			}

			if (cur_line[0] == '#')
			{
				delete[] cur_line;
				continue;
			}

			//1       65407   65409   . N 1_65408_INS_ENST00000641515.2       0       +       EFFECT_PROMOTOR_OVERLAP ENST00000641515.2
			t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, "\t");

			char* cur_var_chrom = cur_line_toks->at(0)->str();
			char* cur_var_posn = cur_line_toks->at(1)->str();
			char* cur_var_impact_str = cur_line_toks->at(6)->str();
			char* cur_var_imp_transcript_id = cur_line_toks->at(7)->str();
			char* cur_var_name = cur_line_toks->at(3)->str();
			t_string_tokens* name_toks = t_string::tokenize_by_chars(cur_var_name, " ");
			char* ref_all = name_toks->at(0)->str();
			char* alt_all = name_toks->at(1)->str();
			t_string_tokens* cur_var_id_toks = name_toks->at(2)->tokenize_by_chars("_");

			// Make sure this is an insert.
			if (!t_string::compare_strings(cur_var_id_toks->at(2)->str(), "INS"))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_var_id_toks);
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			if (line_i % 100000 == 0)
			{
				fprintf(stderr, "@ %d. annotation: %s                    \r", line_i, cur_var_name);
			}

			line_i++;

			// Following converts from BED to VEP coordinates, it is the first nucleotide just before insert.
			int cur_ins_start_pos = atoi(cur_var_posn)+1;

			// Make sure what we are matching is what we generated.
			if (!t_string::compare_strings(cur_var_id_toks->at(3)->str(), cur_var_imp_transcript_id))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_var_id_toks);
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			t_string::clean_tokens(cur_var_id_toks);

			//fprintf(stderr, "Alleles: %c -> %c\n", ref_all[0], alt_all[0]);

			cds_posn_counter++;
			if (cds_posn_counter > 500)
			{
				cds_posn_counter = 30;
			}

			// Map the impact.
			int xcvtr_imp_str_i = t_string::get_i_str(xcvatr_terms, cur_var_impact_str);
			if (xcvtr_imp_str_i == xcvatr_terms->size())
			{
				fprintf(stderr, "Could not find the impact string: %s\n", cur_var_impact_str);
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			char* mapped_vep_impact_str = vep_terms->at(xcvtr_imp_str_i);

			char cur_var_id[1000];
			char cur_var_alt_all[1000];
			if (sscanf(cur_var_name, "%*s %s %s", cur_var_alt_all, cur_var_id) != 2)
			{
				fprintf(stderr, "Could not parse variant identifier from %s\n", cur_var_name);
				exit(0);
			}

			//fprintf(stderr, "PROCESSED_LINE: %s (%d tokens)\n", cur_line, cur_line_toks->size());

			char feat_type[] = "transcript";
			fprintf(f_op, "%s\t%s:%d-%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d-%d\n",
				cur_var_id, cur_var_chrom, cur_ins_start_pos, cur_ins_start_pos+1, cur_var_alt_all,
				cur_var_imp_transcript_id, cur_var_imp_transcript_id, feat_type,
				mapped_vep_impact_str, cds_posn_counter, cds_posn_counter, cds_posn_counter+1);

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			t_string::clean_tokens(name_toks);
			delete[] cur_line;
		} // xcvtr annotation reading loop.

		close_f(f_xcvtr_insertion_annots, XCVATR_insert_annotation_fp);
		close_f(f_op, vep_op_fp);
	} // -convert_Indel_annotations_XCVATR_2_VEP option.
	else if (t_string::compare_strings(argv[1], "-convert_SNV_annotations_XCVATR_2_VEP"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s %s [XCVATR annotated SNVs path] [XCVATR 2 VEP term mapping list] [Output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* XCVATR_snv_annotation_fp = argv[2];
		char* XCVATR_2_VEP_mapping_list_fp = argv[3];
		char* vep_op_fp = argv[4];

		// Load the mapper.
		fprintf(stderr, "Loading XCVATR->VEP term mapping from %s\n", XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* term_mapper_lines = buffer_file(XCVATR_2_VEP_mapping_list_fp);
		vector<char*>* xcvatr_terms = new vector<char*>();
		vector<char*>* vep_terms = new vector<char*>();
		for (int i_l = 0; i_l < term_mapper_lines->size(); i_l++)
		{
			fprintf(stderr, "Parsing: %s\n", term_mapper_lines->at(i_l));

			t_string_tokens* cur_toks = t_string::tokenize_by_chars(term_mapper_lines->at(i_l), "\t");
			if (cur_toks->size() < 2)
			{
				fprintf(stderr, "Could not correctly parse: %d tokens.\n", cur_toks->size());
				exit(0);
			}

			xcvatr_terms->push_back(t_string::copy_me_str(cur_toks->at(0)->str()));
			vep_terms->push_back(t_string::copy_me_str(cur_toks->at(1)->str()));

			fprintf(stderr, "%s => %s\n", xcvatr_terms->back(), vep_terms->back());
		} // i_l loop.

		// Following are the VEP column indices.
		// #Uploaded_variation     Location        Allele  Gene    Feature Feature_type    Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  Existing_variation      Extra
		int VAR_ID_col_i = 0;
		int LOCATION_col_i = 1;
		int ALLELE_col_i = 2;
		int GENE_ID_col_i = 3;
		int TRANSCRIPT_ID_col_i = 4;
		int CONSEQUENCE_col_i = 6;
		int CDS_POSN_col_i = 8;

		int cds_posn_counter = 30;

		//char cumul_var_name[1000];
		//cumul_var_name[0] = 0;
		//char cumul_var_impacting_transcript_id[1000];
		//cumul_var_impacting_transcript_id[0] = 0;
		//vector<char*>* cumul_var_vep_impacts = new vector<char*>();

		int line_i = 0;

		FILE* f_xcvtr_snv_annots = open_f(XCVATR_snv_annotation_fp, "r");
		FILE* f_op = open_f(vep_op_fp, "w");
		while (1)
		{
			char* cur_line = getline(f_xcvtr_snv_annots);
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

			char* cur_var_chrom = cur_line_toks->at(0)->str();
			char* cur_var_posn = cur_line_toks->at(2)->str();
			char* cur_var_impact_str = cur_line_toks->at(6)->str();
			char* cur_var_imp_transcript_id = cur_line_toks->at(7)->str();
			char* cur_var_name = cur_line_toks->at(3)->str();
			t_string_tokens* name_toks = t_string::tokenize_by_chars(cur_var_name, " ");
			char* ref_all = name_toks->at(0)->str();
			char* alt_all = name_toks->at(1)->str();
			t_string_tokens* cur_var_id_toks = name_toks->at(2)->tokenize_by_chars("_");

			if (line_i % 100000 == 0)
			{
				fprintf(stderr, "@ %d. annotation: %s                    \r", line_i, cur_var_name);
			}

			line_i++;

			// Make sure what we are matching is what we generated.
			if (!t_string::compare_strings(cur_var_id_toks->at(3)->str(), cur_var_imp_transcript_id))
			{
				delete[] cur_line;
				t_string::clean_tokens(cur_var_id_toks);
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			t_string::clean_tokens(cur_var_id_toks);

			//fprintf(stderr, "Alleles: %c -> %c\n", ref_all[0], alt_all[0]);

			cds_posn_counter++;
			if (cds_posn_counter > 500)
			{
				cds_posn_counter = 30;
			}

			// Map the impact.
			int xcvtr_imp_str_i = t_string::get_i_str(xcvatr_terms, cur_var_impact_str);
			if (xcvtr_imp_str_i == xcvatr_terms->size())
			{
				fprintf(stderr, "Could not find the impact string: %s\n", cur_var_impact_str);
				delete[] cur_line;
				t_string::clean_tokens(cur_line_toks);
				t_string::clean_tokens(name_toks);
				continue;
			}

			char* mapped_vep_impact_str = vep_terms->at(xcvtr_imp_str_i);

			char cur_var_id[1000];
			char cur_var_alt_all[1000];
			if (sscanf(cur_var_name, "%*s %s %s", cur_var_alt_all, cur_var_id) != 2)
			{
				fprintf(stderr, "Could not parse variant identifier from %s\n", cur_var_name);
				exit(0);
			}

			//fprintf(stderr, "PROCESSED_LINE: %s (%d tokens)\n", cur_line, cur_line_toks->size());

			char feat_type[] = "transcript";
			fprintf(f_op, "%s\t%s:%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n",
					cur_var_id, cur_var_chrom, cur_var_posn, cur_var_alt_all,
					cur_var_imp_transcript_id, cur_var_imp_transcript_id, feat_type,
					mapped_vep_impact_str, cds_posn_counter, cds_posn_counter);

	//		// Map the current XCVTR impact to VEP impact.
	//		if (!t_string::compare_strings(cumul_var_name, cur_var_name) ||
	//			!t_string::compare_strings(cumul_var_impacting_transcript_id, cur_var_imp_transcript_id))
	//		{
	//			// If we already have a variant, print it in VEP format.
	//			if (t_string::string_length(cumul_var_name) > 0)
	//			{
	//				fprintf(stderr, ">>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	//				//fprintf(stderr, "Saving the current cumulative variant effect: %s (%s), %d impacts.\n",
	//				//	cumul_var_name,
	//				//	cumul_var_impacting_transcript_id,
	//				//	cumul_var_vep_impacts->size());

	//				// Save the current cumulative impacts for the current variant on the current transcript.
	//				t_string* cumul_vep_impact = new t_string();
	//				for (int imp_i = 0; imp_i < cumul_var_vep_impacts->size(); imp_i++)
	//				{
	//					if (imp_i > 0)
	//					{
	//						cumul_vep_impact->concat_char(',');
	//					}

	//					cumul_vep_impact->concat_string(cumul_var_vep_impacts->at(imp_i));
	//				} // imp_i loop.

	//				//fprintf(stderr, "Generated cumulative VEP impact string: %s\n", cumul_vep_impact->str());

	//				/*
	//int VAR_ID_col_i = 0;
	//int LOCATION_col_i = 1;
	//int ALLELE_col_i = 2;
	//int GENE_ID_col_i = 3;
	//int TRANSCRIPT_ID_col_i = 4;
	//int CONSEQUENCE_col_i = 6;
	//int CDS_POSN_col_i = 8;
	//				*/
	//				char cur_var_id[1000];
	//				char cur_var_alt_all[1000];
	//				if (sscanf(cumul_var_name, "%*s %s %s", cur_var_alt_all, cur_var_id) != 2)
	//				{
	//					fprintf(stderr, "Could not parse variant identifier from %s\n", cur_var_name);
	//					exit(0);
	//				}

	//				char feat_type[] = "transcript";
	//				fprintf(f_op, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n",
	//					cur_var_id, cur_var_posn, cur_var_alt_all,
	//					cumul_var_impacting_transcript_id, cumul_var_impacting_transcript_id, feat_type,
	//					cumul_vep_impact->str(), cds_posn_counter, cds_posn_counter);

	//				fprintf(stderr, "<<<<<<<<<<<<<<<<<<<<<<<<<<\n");

	//				delete cumul_vep_impact;
	//			} // variant existence check.

	//			// Clear the impacts for the previous impact info.
	//			cumul_var_vep_impacts->clear();

	//			// Update the cumulative variant info.
	//			fprintf(stderr, "Updating the current cumulative variant effect: %s (%s)\n", cur_var_name, cur_var_imp_transcript_id);
	//			strcpy(cumul_var_impacting_transcript_id, cur_var_imp_transcript_id);
	//			strcpy(cumul_var_name, cur_var_name);
	//			cumul_var_vep_impacts->push_back(t_string::copy_me_str(mapped_vep_impact_str));
	//		}
	//		else
	//		{
	//			// Add the mapped impact.
	//			cumul_var_vep_impacts->push_back(t_string::copy_me_str(mapped_vep_impact_str));
	//		}

			/*fprintf(stderr, "PROCESSED_LINE: %s (%d tokens)\n", cur_line, cur_line_toks->size());*/

			// Free memory.
			t_string::clean_tokens(cur_line_toks);
			t_string::clean_tokens(name_toks);
			delete[] cur_line;
		} // xcvtr annotation reading loop.

		close_f(f_xcvtr_snv_annots, XCVATR_snv_annotation_fp);
		close_f(f_op, vep_op_fp);
	} // -convert_SNV_annotations_XCVATR_2_VEP option.
	if (t_string::compare_strings(argv[1], "-plain_aggregate_encoded_SNV_genotype_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encoded genotypes directory] [EOI regions BED path] [VOI regions BED path] [Sample id's list file path] [Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encoded_vectorized_snv_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* VOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		plain_aggregate_encoded_SNV_genotype_Full_Matrix(encoded_vectorized_snv_dir,
			EOI_regs_BED_fp,
			VOI_regs_BED_fp,
			sample_ids_list_fp,
			genotype_encoding_type);
	} // -plain_aggregate_encoded_SNV_genotype_matrix option.
	else if (t_string::compare_strings(argv[1], "-plain_aggregate_encoded_Del_genotype_matrix"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Encoded genotypes directory] [EOI regions BED path] [VOI regions BED path] [Sample id's list file path] [Genotype Encoding (0:Existence/1:Allele Count)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* encoded_vectorized_del_dir = argv[2];
		char* EOI_regs_BED_fp = argv[3];
		char* VOI_regs_BED_fp = argv[4];
		char* sample_ids_list_fp = argv[5];
		int genotype_encoding_type = atoi(argv[6]);

		plain_aggregate_encoded_Del_genotype_matrix(encoded_vectorized_del_dir,
													EOI_regs_BED_fp,
													VOI_regs_BED_fp,
													sample_ids_list_fp,
													genotype_encoding_type);
	} // -plain_aggregate_encoded_Del_genotype_matrix option.
	if (t_string::compare_strings(argv[1], "-encode_DEL_genotypes_2_matrix_per_VCF"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Per chrom VCF dir.] [Sample ids list path] [EOI regions BED path] [Genotype Encoding (0:Existence/1:Allele Count)] [Genotypes output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* per_chrom_VCF_dir = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* EOI_regs_BED_fp = argv[4];
		int genotype_encoding_type = atoi(argv[5]);
		char* genotypes_op_dir = argv[6];

		encode_DEL_genotypes_2_matrix_per_VCF(per_chrom_VCF_dir, sample_ids_list_fp,
			EOI_regs_BED_fp,
			genotype_encoding_type, genotypes_op_dir);
	} // -compare_VEP_annotations_with_SVAT_annotations option.
	else if (t_string::compare_strings(argv[1], "-encode_SNV_genotypes_2_matrix_per_VCF"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Per chrom VCF dir.] [Sample ids list path] [EOI regions BED path] [Genotype Encoding (0:Existence/1:Allele Count)] [Genotypes output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* per_chrom_VCF_dir = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* EOI_regs_BED_fp = argv[4];
		int genotype_encoding_type = atoi(argv[5]);
		char* genotypes_op_dir = argv[6];

		encode_SNV_genotypes_2_matrix_per_VCF(per_chrom_VCF_dir, sample_ids_list_fp,
			EOI_regs_BED_fp,
			genotype_encoding_type, genotypes_op_dir);
	} // -compare_VEP_annotations_with_SVAT_annotations option.
	else if (t_string::compare_strings(argv[1], "-compare_VEP_annotations_with_SVAT_annotations"))
	{ 
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [VEP annotated output file path] [VEP impact strings list (full)] \
[VEP impacts 2 focus list file path] [High impacts list file path] [SVAT annotated BED file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vep_annotated_op_fp = argv[2];
		char* VEP_impact_string_ctx_fp = argv[3];
		char* vep_impacts_2_focus_list_fp = argv[4];
		char* hi_impacts_list_fp = argv[5];
		char* svat_annotated_vcf_fp = argv[6];

		compare_VEP_annotations_with_SVAT_annotations(vep_annotated_op_fp,
														VEP_impact_string_ctx_fp,
														vep_impacts_2_focus_list_fp,
														hi_impacts_list_fp,
														svat_annotated_vcf_fp);
	} // -compare_signal_annotated_vs_direct_VEP_annotations option.
	else if (t_string::compare_strings(argv[1], "-compare_ANNOVAR_vs_SVAT_variant_wise"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [ANNOVAR annotated output file path] [ANNOVAR-2-SVAT impact term mapping list (2 columns)] [ANNOVAR transcript id's list path] [SVAT annotated BED file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annovar_annotated_op_fp = argv[2];
		char* annovar_2_svat_matching_impact_terms_list_fp = argv[3];
		char* annovar_transcripts_list_fp = argv[4];
		char* svat_annotated_vcf_fp = argv[5];

		compare_ANNOVAR_vs_SVAT_variant_wise(annovar_annotated_op_fp,
			annovar_2_svat_matching_impact_terms_list_fp,
			annovar_transcripts_list_fp,
			svat_annotated_vcf_fp);
	} // -compare_ANNOVAR_vs_SVAT_variant_wise option.
	else if (t_string::compare_strings(argv[1], "-compare_ANNOVAR_annotations_with_SVAT_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [ANNOVAR annotated output file path] [ANNOVAR-2-SVAT impact term mapping list (2 columns)] [SVAT annotated BED file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annovar_annotated_op_fp = argv[2];
		char* annovar_2_svat_matching_impact_terms_list_fp = argv[3];
		char* svat_annotated_vcf_fp = argv[4];

		compare_ANNOVAR_annotations_with_SVAT_annotations(annovar_annotated_op_fp, 
			annovar_2_svat_matching_impact_terms_list_fp,
			svat_annotated_vcf_fp);
	} // -compare_signal_annotated_vs_direct_VEP_annotations option.

	else if (t_string::compare_strings(argv[1], "-separate_VCF_per_chromosome"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [VCF file path (stdin ok)] [chromosome ids/lengths list file path] [Per chromosome output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vcf_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* per_chrom_op_dir = argv[4];

		separate_VCF_per_chromosome(vcf_fp, chr_ids_lengths_list_fp, per_chrom_op_dir);
	} // -separate_VEP_output_per_chromosome option.
	else if (t_string::compare_strings(argv[1], "-separate_VEP_output_per_chromosome"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [VEP output file path (stdin ok)] [chromosome ids/lengths list file path] [Per chromosome output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* vep_op_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* per_chrom_op_dir = argv[4];

		separate_VEP_output_per_chromosome(vep_op_fp, chr_ids_lengths_list_fp, per_chrom_op_dir);
	} // -separate_VEP_output_per_chromosome option.
	else if (t_string::compare_strings(argv[1], "-generate_simulated_InDels_on_EOI"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Maximum indel length] [Per position variant Probability] [Deletion Probability (<1.0)] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		int l_max_indel = atoi(argv[3]);
		double per_posn_var_prob = atof(argv[4]);
		double del_prob = atof(argv[5]);
		char* bin_seq_dir = argv[6];
		char* op_vcf_fp = argv[7];

		generate_EOI_randomized_Indel_VCF(EOI_bed_fp, l_max_indel, per_posn_var_prob, del_prob, bin_seq_dir, op_vcf_fp);
	} // -generate_simulated_InDels_on_EOI option.
	else if (t_string::compare_strings(argv[1], "-generate_simulated_SNVs_on_EOI"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Per position variant Probability] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		double per_posn_var_prob = atof(argv[3]);
		char* bin_seq_dir = argv[4];
		char* op_vcf_fp = argv[5];

		generate_EOI_randomized_SNV_VCF(EOI_bed_fp, per_posn_var_prob, bin_seq_dir, op_vcf_fp);
	} // -generate_simulated_InDels_on_EOI option.
	else if (t_string::compare_strings(argv[1], "-extract_annotated_SNVs_from_signal"))
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s %s [Annotated SNV signals directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* op_fp = argv[3];

		extract_annotated_SNVs_from_signal(annotation_signal_dir, op_fp);
	} // -extract_annotated_SNVs_from_signal option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_Deletions_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_Deletions_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Deletions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_Insertions_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_Insertions_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Insertions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-translate_annotated_SNVs_from_annotated_signals"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotated Indels signals directory] [Per chromosome Deletions VCF directory] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_fp = argv[4];

		translate_annotated_SNVs_from_annotated_signals(annotation_signal_dir, per_chrom_VCF_dir, op_fp);
	} // -translate_annotated_Insertions_from_annotated_signals option.
	else if (t_string::compare_strings(argv[1], "-multiply_signalized_SNVs_and_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotation signal directory] [Variant signal directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* variant_signal_dir = argv[3];
		char* op_dir = argv[4];

		multiply_SNV_variant_and_annotation_signals(annotation_signal_dir, variant_signal_dir, op_dir);
	} // -multiply_signalized_SNVs_and_annotations option.
	else if (t_string::compare_strings(argv[1], "-multiply_signalized_Indels_and_annotations"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Annotation signal directory] [Variant signal directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* annotation_signal_dir = argv[2];
		char* variant_signal_dir = argv[3];
		char* op_dir = argv[4];

		multiply_Indel_variant_and_annotation_signals(annotation_signal_dir, variant_signal_dir, op_dir);
	} // -multiply_signalized_Indels_and_annotations option.
	else if (t_string::compare_strings(argv[1], "-get_frame_gene_unique_EOIs_BED_per_CDS_regions"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Element IDs list file path] [Sub-element selector: \"CDS\"/\"exon\"] [GFF file path] [Exon extension length (5 bps)] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* element_ids_list_fp = argv[2];
		char* CDS_exon_selector = argv[3];
		char* gff_fp = argv[4];
		int l_exon_ext = atoi(argv[5]);
		char* op_fp = argv[6];

		get_frame_gene_unique_EOIs_BED_per_CDS_regions(element_ids_list_fp, CDS_exon_selector, gff_fp, l_exon_ext, op_fp);
	}
	else if (t_string::compare_strings(argv[1], "-get_EOIs_BED_per_exonic_regions"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Element IDs list file path] [GFF file path] [CDS/Exon selection: \"CDS\"/\"exon\"] [Element type (\"Gene\"/\"Transcript\")] [Exon extension length (5 bps)] [Output file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* element_ids_list_fp = argv[2];
		char* gff_fp = argv[3];
		char* CDS_exon_selector = argv[4];
		char* EOI_element_type_str = argv[5];
		int l_exon_ext = atoi(argv[6]);
		char* op_fp = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "EOIs with genes.\n");
			EOI_element_type = EOI_GENE_ELEMENTS;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "EOIs with transcripts.\n");
			EOI_element_type = EOI_TRANSCRIPT_ELEMENTS;
		}

		get_EOIs_BED_per_exonic_regions_per_GENCODE_GFF(element_ids_list_fp, gff_fp, CDS_exon_selector, EOI_element_type, l_exon_ext, op_fp);
	} // -get_EOIs_BED_per_exonic_regions
	else if (t_string::compare_strings(argv[1], "-generate_EOI_enumerating_SNV_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		char* bin_seq_dir = argv[3];
		char* op_fp = argv[4];

		generate_EOI_enumerating_SNV_VCF(EOI_bed_fp, bin_seq_dir, op_fp);
	} // -generate_Exome_enumerating_VCF option.
	else if (t_string::compare_strings(argv[1], "-generate_EOI_enumerating_Indel_VCF"))
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [EOI regions BED file path] [Maximum indel length] [Genome Sequence directory path] [Output VCF file path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_bed_fp = argv[2];
		int l_max_indel = atoi(argv[3]);
		char* bin_seq_dir = argv[4];
		char* op_fp = argv[5];

		generate_EOI_enumerating_Indel_VCF(EOI_bed_fp, l_max_indel, bin_seq_dir, op_fp);
	} // -generate_Exome_enumerating_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_Insertion_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Per chromosome VCF directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* VCF_fp = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_Insertion_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			VCF_fp,
			op_dir);
	} // -signalize_Indel_genotypes_per_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_Deletion_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Per chromosome VCF directory] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_Deletion_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			per_chrom_VCF_dir,
			op_dir);
	} // -signalize_Indel_genotypes_per_VCF option.
	else if (t_string::compare_strings(argv[1], "-signalize_SNV_genotypes_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Per chromosome VCF directory] [Output directory]\n", argv[0], argv[1]);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* VCF_fp = argv[3];
		char* op_dir = argv[4];

		signalize_VCF_SNV_genotypes_per_EOI_regs(EOI_regs_BED_fp,
			VCF_fp,
			op_dir);
	} // -signalize_VCF_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_non_coding_variants_per_VCF"))
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [VCF file path (stdin not ok)] [Output directory]\n", argv[0], argv[1]);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* per_chrom_VCF_dir = argv[3];
		char* op_dir = argv[4];

		signalize_non_coding_variants_per_VCF(EOI_regs_BED_fp,
			per_chrom_VCF_dir,
			op_dir);
	} // -signalize_VCF_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_SNVs"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP Output file path (stdin not ok)] \
[Genome sequence directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_fp = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "DO NOT USE SUMMARIZATION OVER GENES, THIS OPTION IS NOT TESTED WELL.\n");exit(0);
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_SNVs_per_EOI_regs_element_summarization(EOI_regs_BED_fp, 
																		genome_dir, 
																		EOI_element_type, 
																		vep_op_fp,
																		sorted_impact_value_strings_list_fp,
																		op_dir);
	} // -signalize_VEP_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Inserts"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP Output file path (stdin not ok)] \
[Genome sequence directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_fp = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Inserts_per_EOI_regs_element_summarization(EOI_regs_BED_fp,
			genome_dir, 
			EOI_element_type,
			vep_op_fp,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.



	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Deletes"))
	{
		if (argc != 8)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP per chromosome output directory] \
[Per chromosome sequences directory] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_dir = argv[4];
		char* genome_dir = argv[5];
		char* sorted_impact_value_strings_list_fp = argv[6];
		char* op_dir = argv[7];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization(EOI_regs_BED_fp,
			genome_dir,
			EOI_element_type,
			vep_op_dir,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.
	else if (t_string::compare_strings(argv[1], "-signalize_VEP_annotated_Deletes_w_CDS_vicinity"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "%s %s [Elements of Interest BED file path] [Element type (\"Gene\"/\"Transcript\")] [VEP per chromosome output directory] \
[Per chromosome sequences directory] [CDS BED file path] [Sorted VEP impact value strings list file path] [Output directory]\n", argv[0], argv[1]);
			exit(0);
		}

		char* EOI_regs_BED_fp = argv[2];
		char* EOI_element_type_str = argv[3];
		char* vep_op_dir = argv[4];
		char* genome_dir = argv[5];
		char* cds_bed_fp = argv[6];
		char* sorted_impact_value_strings_list_fp = argv[7];
		char* op_dir = argv[8];

		int EOI_element_type = -1;
		if (t_string::compare_strings_ci(EOI_element_type_str, "gene"))
		{
			fprintf(stderr, "Summarizing impact over genes.\n");
			EOI_element_type = VEP_SIGNALIZE_GENE_SUMMARIZATION;
		}
		else if (t_string::compare_strings_ci(EOI_element_type_str, "transcript"))
		{
			fprintf(stderr, "Summarizing impact over transcripts.\n");
			EOI_element_type = VEP_SIGNALIZE_TRANSCRIPT_SUMMARIZATION;
		}

		signalize_VEP_annotated_Deletes_per_EOI_regs_element_summarization_with_CDS_vicinity(EOI_regs_BED_fp,
			genome_dir, cds_bed_fp,
			EOI_element_type,
			vep_op_dir,
			sorted_impact_value_strings_list_fp,
			op_dir);
	} // -signalize_VEP_annotated_variants option.
	else if (strcmp(argv[1], "-binarize_fasta_per_file") == 0)
	{
		if (argc != 4)
		{
			fprintf(stderr, "%s -binarize_fasta_per_file [FASTA file path] [Output directory]\n", argv[0]);
			exit(0);
		}

		char* fasta_fp = argv[2];
		char* op_dir = argv[3];
		binarize_fasta_file(fasta_fp, op_dir);
	}
	

}
