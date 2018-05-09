#ifndef __VARIANTS
#define __VARIANTS

#include <stdio.h>
#include <stdbool.h>

#include "processrefgen.h"

#define add_to_beginning 0
#define  totalNumInd 20
#define  strSize 200 //maximum size of a constant string, such as readName length


typedef struct strvar {
	char *chr_name; /*chromosome name*/
	int outer_start; /* outer start coordinate */
	int inner_start; /* inner start coordinate */
	int outer_end; /* outer end coordinate */
	int inner_end; /* inner end coordinate */
	char svtype; /* type of SV */
	char *mei_name; /* the name of the mobile element */
	bool filtered; /* True if the variation is filtered in the filtering stage */
	bool mei_del;
	int rp[totalNumInd]; /* RP Support for each individual */
	int sr[totalNumInd]; /* SR Support for each individual */
	double cnv_score[totalNumInd]; // likelihood score
	bool low_qual;
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

struct strvar ** init_vars(int num_chroms);
struct strvar* new_strvar(char *chrName, int outer_start, int inner_start, int outer_end, int inner_end, char svtype,
		bool filtered, bool mei_del, char *mei_name, double cnv_score[], int rp[], int sr[]);
void add_strvar(struct strvar ** variations, struct strvar* sv);
void print_strvar(bam_info** in_bams, parameters* params, struct strvar* sv, FILE* fpOut);
int print_all_vars(struct strvar ** variations, parameters *params, FILE *fpOut);
void free_quick(bam_info** in_bams, parameters *params);
void free_chr(struct strvar ** variations, char* chr_);
void print_vcf_header( FILE *fpOut, bam_info** in_bams, parameters *params);
int print_chr(struct strvar ** variations, parameters *params, char* chr_, FILE *fpOut);
void print_sv_stats();
#endif
