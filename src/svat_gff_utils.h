#ifndef __GFF_UTILS__
#define __GFF_UTILS__

#include <vector>
using namespace std;

struct t_annot_region;

struct t_gene_annotation_info
{
	int val;
};

struct t_transcript_annotation_info
{
	char* gene_name;
	t_annot_region* gene;
	t_annot_region* fp_UTR;
	t_annot_region* tp_UTR;
	vector<t_annot_region*>* CDSs;
	bool valid_3mer_frame;
	int l_transcript;
};

struct t_intron_annotation_info
{
	char* gene_name;
	t_annot_region* gene;
	char* sequence;
};

struct t_CDS_annotation_info
{
	t_annot_region* gene;

	t_annot_region* transcript;
	char* sequence;

	// Coding frame index: 0,1,2 corresponding to which nucleotide of the frame.
	int coding_frame_i_per_first_nuc;

	int l_transcript_so_far;
};

struct t_UTR_annotation_info
{
	t_annot_region* transcript;
};

struct t_promotor_annotation_info
{
	t_annot_region* transcript;
};

// All the annotation information relating GENCODE GFF is set here.
void load_elements_set_annotation_information_per_GENCODE_GFF(char* gencode_gff_fp, char* genome_bin_seq_dir,
																vector<t_annot_region*>* gene_regs,
																vector<t_annot_region*>* transcript_regs,
																vector<t_annot_region*>* cds_regs,														
																vector<t_annot_region*>* UTR_5p_regs,
																vector<t_annot_region*>* UTR_3p_regs,
																vector<t_annot_region*>* intronic_regs,
																vector<t_annot_region*>* promotor_regs, int l_promotor);

vector<t_annot_region*>* load_GFF(char* gff_fp);
vector<t_annot_region*>* load_GFF_with_line_information(char* gff_fp);

vector<t_annot_region*>* get_gff_regions_per_feature_type(vector<t_annot_region*>* all_regions,
															char* feature_type);

vector<t_annot_region*>* get_gff_regions_per_source_string(vector<t_annot_region*>* all_regions,
															char* source_str);

//vector<t_annot_region*>* get_gff_feature_type(t_annot_region* region);

//char* get_gff_property_per_grp_str(t_annot_region* entry, char* (*get_prop_per_group_str)(char*));

struct t_gff_info
{
	char* feature_type; // Type of the feature: exon, TSS, gene, transcripts, start_codon, ...
	char* group_str; // The information that identifies the entry.
	char* source_str; // This is the experiment id for this entry.
	vector<char*>* prop_ids;
	vector<char*>* prop_vals;
};

// Extract the properties for each entry in the GFF entry list. The callback extracts the property string. This function forms the interface between
// the external functions and the GFF property entries, i.e., group_str entry.
vector<char*>* get_gff_string_properties_per_entry(vector<t_annot_region*>* all_regions,
												char* prop_id);

char* get_prop_val_per_prop_id(t_annot_region* region, char* prop_id);

void parse_group_strs(vector<t_annot_region*>* all_regions, void (*parse_grp_str)(t_annot_region*, vector<char*>*, vector<char*>*));

// Go over all the regions, get the regions with the same property together, choose the largest region that contains the all regions with the selected property.
void form_largest_regions_per_property_type(vector<t_annot_region*>* all_regions,
											vector<t_annot_region*>* merged_selected_regions,
											vector<char*>* sel_region_props, 
											char* prop_id);

void associate_gff_children_2_parents(vector<t_annot_region*>* parent_regions,
	vector<t_annot_region*>* children_regions,
	vector<char*>* (*get_parent_feature_id_per_child_group_str)(t_annot_region*), // Use this to extract the id for the parent using the group string for child.
	char* (*get_parent_feature_id_per_parent_group_str)(t_annot_region*)); // Use this to extract the id for the parent using the group string for parent.

t_annot_region* copy_gff_entry(t_annot_region* gff_entry);

#endif // __GFF_UTILS__

