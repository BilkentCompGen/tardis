#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "variants.h"
#include <htslib/faidx.h>

int del_cnt = 0;
int ins_cnt = 0;
int inv_cnt = 0;
int mei_cnt = 0;
int numt_cnt = 0;
int tandup_cnt = 0;
int invdup_cnt = 0;
int interdup_cnt = 0;
int dup_cnt = 0;
int mei_cnt_filtered = 0;

long total_del_length = 0;
long total_ins_length = 0;
long total_inv_length = 0;
long total_mei_length = 0;
long total_numt_length = 0;
long total_tandup_length = 0;
long total_invdup_length = 0;
long total_interdup_length = 0;


int sv_count;
int sv_lowqual_count;
int indCount;

// create a new structure and return address
struct strvar* new_strvar(char *chrName, int outer_start, int inner_start, int outer_end, int inner_end, char svtype,
		bool filtered, bool mei_del, char *mei_name, char *mei_type, double cnv_score[], int rp[], int sr[],
		double homogeneity_score, float weight, int zygosity[], bool imprecise)
{
	int i;
	struct strvar* a_strvar = getMem( sizeof( struct strvar));

	a_strvar->chr_name = chrName;
	a_strvar->outer_start = outer_start;
	a_strvar->inner_start = inner_start;
	a_strvar->outer_end = outer_end;
	a_strvar->inner_end = inner_end;
	a_strvar->svtype = svtype;
	a_strvar->mei_name = mei_name;
	a_strvar->mei_type = mei_type;
	a_strvar->filtered = filtered;
	a_strvar->mei_del = mei_del;
	a_strvar->imprecise = imprecise;
	a_strvar->homogeneity_score = homogeneity_score;
	a_strvar->weight = weight;

	for( i = 0; i < indCount; i++)
	{
		a_strvar->cnv_score[i] = cnv_score[i];
		a_strvar->rp[i] = rp[i];
		a_strvar->sr[i] = sr[i];
		if( zygosity[i] == 0)
			a_strvar->zygosity[i] = "0/0";
		else if( zygosity[i] == 1)
			a_strvar->zygosity[i] = "0/1";
		else if( zygosity[i] == 2)
			a_strvar->zygosity[i] = "1/0";
		else if( zygosity[i] == 3)
			a_strvar->zygosity[i] = "1/1";
	}

	return a_strvar;
}


char* readRefAltSeq( parameters *params, char* chr_name, int start, int end)
{
	int i, min, max, loc_length, chr_index;
	char *ref_seq;
	long bp_cnt = 0;
	faidx_t* ref_fai;

	chr_index = sonic_refind_chromosome_index( params->this_sonic, chr_name);
	min = start;
	max = end;

	ref_fai = fai_load( params->ref_genome);
	ref_seq = faidx_fetch_seq( ref_fai, params->this_sonic->chromosome_names[chr_index], min, max, &loc_length);

	fai_destroy( ref_fai);

	return ref_seq;
}

char* readRefAltSeqMEI( parameters *params, char* chr_name, char *mei_string)
{
	sonic_interval *this_interval;

	char *seq;
	int pos_start, pos_end, chr_index;

	chr_index = sonic_refind_chromosome_index( params->this_sonic, chr_name);

	pos_start = 0;
	pos_end = 499;

	while( pos_end < params->this_sonic->chromosome_lengths[chr_index])
	{
		this_interval = sonic_intersect( params->this_sonic, chr_name, pos_start, pos_end, SONIC_REP);

		if( this_interval != NULL)
		{
			if( strcmp( this_interval->repeat_item->repeat_type, mei_string) == 0)
			{
				seq = readRefAltSeq( params, chr_name, this_interval->start - 1, this_interval->end - 1);
				return seq;
			}
		}
		pos_start += 500;
		pos_end += 500;
	}
	return NULL;
}

char* reverseComplement( char* str)
{
	int i;
	char* str2 = NULL;

	set_str( &str2, str);
	reverse_string( str2);

	for(i = 0; i < strlen(str2); i++)
	{
		char tmp = complement_char(str2[i]);
		str2[i] = tmp;
	}
	return str2;
}


void print_sv_stats()
{
	fprintf(logFile,"\n\nTARDIS found %d SVs total\n", del_cnt + ins_cnt + inv_cnt + tandup_cnt + mei_cnt + interdup_cnt + invdup_cnt);
	fprintf(logFile,"\tDeletion: %d\n", del_cnt);
	fprintf(logFile,"\tInsertion: %d\n", ins_cnt);
	fprintf(logFile,"\tInversion: %d\n", inv_cnt);
	fprintf(logFile,"\tTandem Duplication: %d\n", tandup_cnt);
	fprintf(logFile,"\tInterspersed Duplication: %d\n", interdup_cnt);
	fprintf(logFile,"\tInterspersed (Inverted) Duplication: %d\n", invdup_cnt);
	fprintf(logFile,"\tMEI: %d (%d filtered)\n", mei_cnt, mei_cnt_filtered);
	fprintf(logFile,"\tNUMT: %d\n", numt_cnt);

	fprintf(stderr,"\n\nTARDIS is complete. Found %d SVs total\n", del_cnt + ins_cnt + inv_cnt + tandup_cnt + mei_cnt + interdup_cnt + invdup_cnt);
	fprintf(stderr,"\tDeletion: %d\n", del_cnt);
	fprintf(stderr,"\tInsertion: %d\n", ins_cnt);
	fprintf(stderr,"\tInversion: %d\n", inv_cnt);
	fprintf(stderr,"\tTandem Duplication: %d\n", tandup_cnt);
	fprintf(stderr,"\tInterspersed Duplication: %d\n", interdup_cnt);
	fprintf(stderr,"\tInterspersed (Inverted) Duplication: %d\n", invdup_cnt);
	fprintf(stderr,"\tMEI: %d (%d filtered)\n", mei_cnt, mei_cnt_filtered);
	fprintf(stderr,"\tNUMT: %d\n", numt_cnt);
}

//add the variation in ascending order according to inner_start
void print_strvar( bam_info** in_bams, parameters* params, struct strvar* sv, FILE* fpOut)
{
	int sv_len = 0, i, j, control, ind_id, rp_total = 0, sr_total = 0;
	char* seq = ".", *seq_rev = ".";

	/* Sum the rp support for each individual with GT 0/1 */
	for( i = 0; i < params->num_bams; i++)
	{
		if( in_bams[i]->contribution == true)
		{
			rp_total += sv->rp[i];
			sr_total += sv->sr[i];
		}
	}
	if( rp_total <= 0)
		sv->filtered = true;

	if( sv->svtype == DELETION)
	{
		/* Find ref and alt sequences */
		seq = readRefAltSeq( params, sv->chr_name, sv->inner_start, sv->inner_end);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%c\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, seq, seq[0], ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
		{
			if(sv->mei_del == true)
				fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, seq[0], "<", "DEL:ME:", sv->mei_type, ">",( sv->filtered == false) ? "PASS" : "LowQual");
			else
				fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, seq[0], "<", "DEL", ">",( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else
		{
			if(sv->mei_del == true)
				fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, ".", "<", "DEL:ME:", sv->mei_type, ">",( sv->filtered == false) ? "PASS" : "LowQual");
			else
				fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, ".", "<", "DEL", ">",( sv->filtered == false) ? "PASS" : "LowQual");

		}
		if( seq != NULL)
			free( seq);

		sv_len = abs( ( sv->inner_end - sv->inner_start + 1));
		total_del_length += sv_len;
	}
	else if( sv->svtype == INSERTION)
	{
		/* Find ref and alt sequences */
		seq = readRefAltSeq( params, sv->chr_name, sv->inner_start, sv->inner_end);

		if( params->seq_resolved != 0 && seq != NULL)
		{
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_ins_", ++ins_cnt, seq[0], ".", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_ins_", ++ins_cnt, seq[0], "<", "INS", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_ins_", ++ins_cnt, ".", "<", "INS", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		sv_len = abs( ( sv->inner_end - sv->inner_start + 1));
		total_ins_length += sv_len;
	}
	else if( sv->svtype == INVERSION)
	{
		seq = readRefAltSeq( params, sv->chr_name, sv->outer_start, sv->outer_end);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			/* Find ref and alt sequences */
			seq_rev = reverseComplement( seq);

			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_inv_", ++inv_cnt, seq, seq_rev, 255, ( sv->filtered == false) ? "PASS" : "LowQual");

			if( seq_rev != NULL)
				free( seq_rev);
		}
		else if( seq != NULL)
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_inv_", ++inv_cnt, seq[0], "<", "INV", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_inv_", ++inv_cnt, ".", "<", "INV", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		sv_len = abs( ( sv->inner_end - sv->inner_start + 1));
		total_inv_length += sv_len;
	}
	else if( sv->svtype == MEIFORWARD)
	{
		seq = readRefAltSeqMEI( params, sv->chr_name, sv->mei_name);

		if( seq != NULL)
			sv_len = strlen(seq);

		if( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE)
		{
			if( seq != NULL)
				fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_mei_", ++mei_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "mfilt");
			else
				fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_mei_", ++mei_cnt, ".", ".", 255, ( sv->filtered == false) ? "PASS" : "mfilt");

			if( sv->filtered)
				mei_cnt_filtered++;
		}
		else if( seq != NULL)
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_mei_", ++mei_cnt, seq[0], "<", "INS:ME:", sv->mei_name, ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");
		else
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_mei_", ++mei_cnt, ".", "<", "INS:ME:", sv->mei_name, ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");

		if( seq != NULL)
			free( seq);

		//sv_len = abs( ( sv->outer_end - sv->inner_end + 1));
		total_mei_length += sv_len;
	}
	else if( sv->svtype == MEIREVERSE)
	{
		seq = readRefAltSeqMEI( params, sv->chr_name, sv->mei_name);
		if(seq != NULL)
			sv_len = strlen(seq);

		if( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE)
		{
			if( seq != NULL)
			{
				seq_rev = reverseComplement( seq);
				fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_end) + 1, "vh_mei_", ++mei_cnt, seq[0], seq_rev, 255, ( sv->filtered == false) ? "PASS" : "mfilt");
			}
			else
			{
				seq_rev = NULL;
				fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_end) + 1, "vh_mei_", ++mei_cnt, ".", ".", 255, ( sv->filtered == false) ? "PASS" : "mfilt");
			}
			if( sv->filtered)
				mei_cnt_filtered++;

			if( seq_rev != NULL)
				free( seq_rev);
		}
		else if( seq != NULL)
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_end) + 1, "vh_mei_", ++mei_cnt, seq[0], "<", "INS:ME:", sv->mei_name, ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");
		else
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_end) + 1, "vh_mei_", ++mei_cnt, ".", "<", "INS:ME:", sv->mei_name, ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");

		if( seq != NULL)
			free( seq);

		//sv_len = abs( ( sv->outer_end - sv->inner_end + 1));
		total_mei_length += sv_len;
	}
	else if( sv->svtype == NUMTFORWARD || sv->svtype == NUMTREVERSE)
	{
		if( params->seq_resolved != 0)
		{
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_numt_", ++numt_cnt, ".", ".", 255, ( sv->filtered == false) ? "PASS" : "mfilt");
		}
		fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1, "vh_numt_", ++numt_cnt, ".", "<", "INS:MT", ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");

		sv_len = abs( ( sv->outer_end - sv->inner_end + 1));
		total_numt_length += sv_len;
	}

	else if( sv->svtype == TANDEMDUP)
	{
		seq = readRefAltSeq( params, sv->chr_name, sv->outer_start, sv->outer_end);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut, "%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], "<", "DUP:TANDEM", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, ".", "<", "DUP:TANDEM", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		tandup_cnt++;
		sv_len = abs( ( sv->outer_end - sv->outer_start + 1));
		total_tandup_length += sv_len;
	}
	else if( sv->svtype == INVDUPRIGHT)
	{
		if( sv->inner_start > sv->outer_start)
			seq = readRefAltSeq( params, sv->chr_name, sv->outer_start, sv->inner_start);
		else
			seq = readRefAltSeq( params, sv->chr_name, sv->inner_start, sv->outer_start);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, ".", "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		invdup_cnt++;
		sv_len = abs( ( sv->inner_start - sv->outer_start + 1));
		total_invdup_length += sv_len;
	}
	else if( sv->svtype == INVDUPLEFT )
	{
		if( sv->outer_end > sv->inner_end)
			seq = readRefAltSeq( params, sv->chr_name, sv->inner_end, sv->outer_end);
		else
			seq = readRefAltSeq( params, sv->chr_name, sv->outer_end, sv->inner_end);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, seq[0], "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, ".", "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		invdup_cnt++;
		sv_len = abs( ( sv->inner_end - sv->outer_end + 1));
		total_invdup_length += sv_len;
	}
	else if( sv->svtype == INTERDUPRIGHT)
	{
		if( sv->inner_start > sv->outer_start)
			seq = readRefAltSeq( params, sv->chr_name, sv->outer_start, sv->inner_start);
		else
			seq = readRefAltSeq( params, sv->chr_name, sv->inner_start, sv->outer_start);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq!= NULL)
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, seq[0], "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->outer_start) + 1," vh_dup_", ++dup_cnt, ".", "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		interdup_cnt++;
		sv_len = abs( ( sv->inner_start - sv->outer_start + 1));
		total_interdup_length += sv_len;
	}
	else if( sv->svtype == INTERDUPLEFT )
	{
		if( sv->outer_end > sv->inner_end)
			seq = readRefAltSeq( params, sv->chr_name, sv->inner_end, sv->outer_end);
		else
			seq = readRefAltSeq( params, sv->chr_name, sv->outer_end, sv->inner_end);

		if( ( params->seq_resolved != 0 || sv->imprecise == NOT_IMPRECISE) && seq != NULL)
		{
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, seq[0], seq, 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else if( seq != NULL)
			fprintf( fpOut,"%s\t%i\t%s%d\t%c\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, seq[0], "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		else
			fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_end) + 1," vh_dup_", ++dup_cnt, ".", "<", "DUP:ISP", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");

		if( seq != NULL)
			free( seq);

		interdup_cnt++;
		sv_len = abs( ( sv->inner_end - sv->outer_end + 1));
		total_interdup_length += sv_len;
	}
	else
		fprintf(stderr, "Problemmmm %c\n", sv->svtype);

	if( sv->svtype == MEIFORWARD)
		fprintf( fpOut, "END=%d;SVLEN=%d;MEINFO=%s;RPSUP=%d;SRSUP=%d;", (sv->outer_start) + sv_len + 1, sv_len, sv->mei_name, rp_total, sr_total);
	if( sv->svtype == MEIREVERSE)
		fprintf( fpOut, "END=%d;SVLEN=%d;MEINFO=%s;RPSUP=%d;SRSUP=%d;", (sv->outer_end) + sv_len + 1, sv_len, sv->mei_name, rp_total, sr_total);
	else if( sv->svtype == INVERSION)
		fprintf( fpOut, "END=%d;SVLEN=%d;RPSUP=%d;SRSUP=%d;", (sv->outer_end) + 1, sv_len, rp_total, sr_total);
	else if( sv->svtype == INVDUPLEFT || sv->svtype == INTERDUPLEFT)
		fprintf( fpOut, "END=%d;SVLEN=%d;POS2=%d;RPSUP=%d;SRSUP=%d;", (sv->outer_end) + 1, sv_len, (sv->outer_start), rp_total, sr_total);
	else if( sv->svtype == INVDUPRIGHT || sv->svtype == INTERDUPRIGHT)
		fprintf( fpOut, "END=%d;SVLEN=%d;POS2=%d;RPSUP=%d;SRSUP=%d;", (sv->inner_start) + 1, sv_len, (sv->outer_end), rp_total, sr_total);
	else if( sv->svtype == TANDEMDUP)
		fprintf( fpOut, "END=%d;SVLEN=%d;RPSUP=%d;SRSUP=%d;", (sv->outer_end) + 1, sv_len, rp_total, sr_total);
	else if( sv->svtype == DELETION)
		fprintf( fpOut, "END=%d;SVLEN=%d;RPSUP=%d;SRSUP=%d;", (sv->inner_end) + 1, -sv_len, rp_total, sr_total);
	else
		fprintf( fpOut, "END=%d;SVLEN=%d;RPSUP=%d;SRSUP=%d;", (sv->inner_end) + 1, sv_len, rp_total, sr_total);

	/* If SV is imprecise */
	if( sv->imprecise == IMPRECISE)
		fprintf( fpOut, "CIEND=0,%d;CIPOS=-%d,0;IMPRECISE;", ( sv->outer_end - sv->inner_end), ( sv->inner_start - sv->outer_start));
	else if( sv->imprecise == START_IMPRECISE)
		fprintf( fpOut, "CIEND=0,0;CIPOS=-%d,0;IMPRECISE;", ( sv->inner_start - sv->outer_start));
	else if( sv->imprecise == END_IMPRECISE)
		fprintf( fpOut, "CIEND=0,%d;CIPOS=0,0;IMPRECISE;", ( sv->outer_end - sv->inner_end));
	else
		fprintf( fpOut, "CIEND=0,0;CIPOS=0,0;PRECISE;");

	/* print TARDIS version info */
	fprintf( fpOut, "SVALG=TARDIS_v%s;", TARDIS_VERSION);

	/* print SV type */
	if( ( sv->svtype == DELETION) && ( sv->mei_del == true))
		fprintf( fpOut, "SVTYPE=DEL:ME:%s\t", sv->mei_type);
	else if( sv->svtype == DELETION)
		fprintf( fpOut, "SVTYPE=DEL\t");
	else if( sv->svtype == MEIFORWARD || sv->svtype == MEIREVERSE)
		fprintf( fpOut, "SVTYPE=INS:ME:%s\t", sv->mei_type);
	else if( sv->svtype == NUMTFORWARD || sv->svtype == NUMTREVERSE)
		fprintf( fpOut, "SVTYPE=INS:MT\t");
	else if( sv->svtype == TANDEMDUP)
		fprintf( fpOut, "SVTYPE=DUP:TANDEM\t");
	else if( sv->svtype == INVDUPLEFT || sv->svtype == INVDUPRIGHT)
		fprintf( fpOut, "SVTYPE=DUP:ISP;ISINV\t");
	else if( sv->svtype == INTERDUPLEFT || sv->svtype == INTERDUPRIGHT)
		fprintf( fpOut, "SVTYPE=DUP:ISP\t");
	else if( sv->svtype == INSERTION)
		fprintf( fpOut, "SVTYPE=INS\t");
	else if( sv->svtype == INVERSION)
		fprintf( fpOut, "SVTYPE=INV\t");

	/* Format field */
	if (params->ten_x || params->output_hs)
		fprintf( fpOut, "GT:CNVL:RP:SR:HS:WE");
	else
		fprintf( fpOut, "GT:CNVL:RP:SR");

	control = 0;
	fprintf( fpOut, "\t");
	for( j = 0; j < params->num_bams; j++)
	{
		if( in_bams[j]->contribution == false)
		{
			if (params->ten_x || params->output_hs)
				fprintf( fpOut, "0/0:%2.6lf:%d:%d:%8.6f:%8.10f\t", sv->cnv_score[j], sv->rp[j], sv->sr[j], sv->homogeneity_score, sv->weight);
			else 
				fprintf( fpOut, "0/0:%2.6lf:%d:%d\t", sv->cnv_score[j], sv->rp[j], sv->sr[j]);
		}
		else
		{
			if (params->ten_x || params->output_hs)
				fprintf( fpOut, "%s:%2.6lf:%d:%d:%8.6f:%8.10f\t", sv->zygosity[j], sv->cnv_score[j], sv->rp[j], sv->sr[j], sv->homogeneity_score, sv->weight);
			else
				fprintf( fpOut, "%s:%2.6f:%d:%d\t", sv->zygosity[j], sv->cnv_score[j], sv->rp[j], sv->sr[j]);
		}
	}
	fprintf( fpOut, "\n");

	if( sv->filtered == true)
		sv_lowqual_count++;
	else
		sv_count++;
}


void print_vcf_header( FILE *fpOut, bam_info** in_bams, parameters *params)
{
	int i;
	/*char header_info[]="##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n"
			"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n"
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
			"##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample ID\">\n"
			"##INFO=<ID=SVALG,Number=1,Type=String,Description=\"SV discovery algorithm\">\n"
			"##INFO=<ID=RPSUP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n"
			"##INFO=<ID=SRSUP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n"
			"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n"
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
			"##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n"
			"##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME\">\n"
			"##INFO=<ID=METRANS,Number=4,Type=String,Description=\"Mobile element transduction info of the form CHR\">\n"
			"##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=\"Indicates a novel structural variation\">\n"
			"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
			"##INFO=<ID=SVMETHOD,Number=.,Type=String,Description=\"Type of approach used to detect SV: RP (read pair), RD (read depth), SR (split read), or AS (assembly)\">\n"
			"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	 */
	char header_info[]="##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n"
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
			"##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample ID\">\n"
			"##INFO=<ID=SVALG,Number=1,Type=String,Description=\"SV discovery algorithm\">\n"
			"##INFO=<ID=RPSUP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n"
			"##INFO=<ID=SRSUP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n"
			"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n"
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
			"##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n"
			"##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME\">\n"
			"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
			"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
			"##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Insertion loci of the duplicated segment\">\n"
			"##INFO=<ID=ISINV,Number=0,Type=Flag,Description=\"Duplicated segment is in inverted orientation\">\n";


	char header_filter[]="##FILTER=<ID=LowQual,Description=\"Genotype call confidence below LOD 1.3\">\n"
			"##FILTER=<ID=mfilt,Description=\"Mobile element filter\">\n";

	char header_format[] = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=CNVL,Number=1,Type=Float,Description=\"CNV Likelihood\">\n"
			"##FORMAT=<ID=HS,Number=1,Type=Float,Description=\"10x Barcode Homogeneity Score\">\n"
			"##FORMAT=<ID=WE,Number=1,Type=Float,Description=\"SV Weight\">\n"
			"##FORMAT=<ID=RP,Number=1,Type=Integer,Description=\"Read Pair Support\">\n"
			"##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Split Read Support\">\n";

	char header_alt[]="##ALT=<ID=DEL,Description=\"Deletion\">\n"
			"##ALT=<ID=DEL:ME,Description=\"Deletion of mobile element\">\n"
			"##ALT=<ID=INS,Description=\"Novel sequence insertion\">\n"
			"##ALT=<ID=INS:ME,Description=\"Mobile element insertion\">\n"
			"##ALT=<ID=MT,Description=\"Nuclear mitochondrial insertion\">\n"
			"##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">\n"
			"##ALT=<ID=DUP:ISP,Description=\"Interspersed duplication\">\n"
			"##ALT=<ID=INV,Description=\"Inversion\">\n";
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	fprintf(fpOut,"##fileformat=VCFv4.3\n");
	fprintf(fpOut, "##fileDate=%d%s%d%s%d\n", timeinfo->tm_year+1900, (timeinfo->tm_mon+1<10 ? "0" : ""), timeinfo->tm_mon+1, (timeinfo->tm_mday<10 ? "0" : ""), timeinfo->tm_mday);
	//fprintf(fpOut,"##reference=%s\n", params->ref_genome);
	fprintf(fpOut,"##reference=%s\n", params->sonic_file);
	fprintf(fpOut, "%s%s%s%s", header_info,header_filter,header_format,header_alt);

	/* Contig field */
	for( i = params->first_chr; i <= params->last_chr; i++)
		fprintf(fpOut, "##contig=<ID=%s,length=%d>\n", params->this_sonic->chromosome_names[i], params->this_sonic->chromosome_lengths[i]);

	fprintf(fpOut, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s","#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");

	for( i = 0; i < params->num_bams; i++)
	{
		fprintf( fpOut, "\t%s",in_bams[i]->sample_name);
	}
	fprintf( fpOut, "\n");
}
