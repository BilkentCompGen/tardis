#ifndef __VARIANTS
#define __VARIANTS

#include <stdio.h>
#include <stdbool.h>
#include "refgenome.h"

#define add_to_beginning 0
#define  totalNumInd 20
#define  strSize 200 //maximum size of a constant string, such as readName length

enum SVTYPE {
	DEL, INS, INV, TANDUP, MEI, TRANSCHROM
};
enum SVTYPE getEnum(char c);
char * svtypeToChar(enum SVTYPE svt);

typedef struct strvar {
	char *chr_name; /*chromosome name*/
	int outer_start; /* outer start coordinate */
	int inner_start; /* inner start coordinate */
	int outer_end; /* outer end coordinate */
	int inner_end; /* inner end coordinate */
	enum SVTYPE svtype; /* type of SV */
	float avg_edit; /* average edit distance of reads mapped to reference */
	int min_svlen; /* lower bound of SV size */
	int max_svlen; /* upper bound of SV size */
	char *samples; /* list of samples that carry the SV. Might need a separate linked list */
	char *mei_name; /* the name of the mobile element */
	double conf_score; /* confidence score for the called variant */
	bool filtered; /* True if the variation is filtered in the filtering stage */
	bool mei_del;
	int rp[totalNumInd]; /* RP Support for each individual */
	int sr[totalNumInd]; /* SR Support for each individual */
	long depth[totalNumInd]; /* Total Read depth of the variation */
	float cn[totalNumInd];
	float del_likelihood[totalNumInd];
	bool low_qual;
	double homogeneity_score;

	struct strvar *next; /* next pointer for linked list */
	struct strvar **head; /* head pointer points to pointer that points to head element */
}strvar;

extern struct strvar ** vars;
extern int indCount;
extern int sv_count;
extern int sv_lowqual_count;

typedef struct chr_index {
	char *chrom_name;
	int index;
	struct chr_index *next;
}chr_index;

/* functions */
int rd_filtering(struct strvar ** variations, bam_info** in_bams, ref_genome* ref, char* dups_file);
struct strvar ** init_vars(int num_chroms);
struct strvar* new_strvar(char *chrName,int outer_start, int inner_start, int outer_end, int inner_end, enum SVTYPE svtype,
		float avg_edit, int min_svlen, int max_svlen, char *samples, double conf_score, bool filtered, bool,
		char *mei_name, long depth[], float cn[], double del_likelihood[], int rp[], int sr[], double homogeneity_score);
void add_strvar(struct strvar ** variations, struct strvar* sv);
void print_strvar(bam_info** in_bams, parameters* params, struct strvar* sv, FILE* fpOut);
int print_all_vars(struct strvar ** variations, parameters *params, FILE *fpOut);
void free_quick(bam_info** in_bams, parameters *params, ref_genome* ref);
void free_sensitive(bam_info** in_bams, parameters *params, ref_genome* ref);
void free_chr(struct strvar ** variations, char* chr_);
void print_vcf_header( FILE *fpOut, bam_info** in_bams, parameters *params);
int print_chr(struct strvar ** variations, parameters *params, char* chr_, FILE *fpOut);
void print_sv_stats();
#endif
