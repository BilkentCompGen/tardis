#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "variants.h"
#include "vh_intervalhandler.h"

int del_cnt = 0;
int ins_cnt = 0;
int inv_cnt = 0;
int mei_cnt = 0;
int tandup_cnt = 0;
int mei_cnt_filtered = 0;

long total_del_length = 0;
long total_ins_length = 0;
long total_inv_length = 0;
long total_mei_length = 0;
long total_tandup_length = 0;

int sv_count;
int sv_lowqual_count;
int indCount;

// create a new structure and return address
struct strvar* new_strvar(char *chrName, int outer_start, int inner_start, int outer_end, int inner_end, char svtype,
		float avg_edit, int min_svlen, int max_svlen, char *samples, double conf_score, bool filtered, bool mei_del,
		char *mei_name, long depth[], float cn[], double del_likelihood[], int rp[], int sr[], double homogeneity_score,
		float weight)
{
	int i;
	struct strvar* a_strvar = getMem( sizeof( struct strvar));

	a_strvar->chr_name = chrName;
	a_strvar->outer_start = outer_start;
	a_strvar->inner_start = inner_start;
	a_strvar->outer_end = outer_end;
	a_strvar->inner_end = inner_end;
	a_strvar->svtype = svtype;
	a_strvar->avg_edit = avg_edit;
	a_strvar->min_svlen = min_svlen;
	a_strvar->max_svlen = max_svlen;
	a_strvar->samples = samples;
	a_strvar->mei_name = mei_name;
	a_strvar->conf_score = conf_score;
	a_strvar->next = NULL;
	a_strvar->head = NULL;
	a_strvar->filtered = filtered;
	a_strvar->mei_del = mei_del;
	a_strvar->homogeneity_score = homogeneity_score;
	a_strvar->weight = weight;

	for( i = 0; i < indCount; i++)
	{
		a_strvar->depth[i] = depth[i];
		a_strvar->cn[i] = cn[i];
		a_strvar->del_likelihood[i] = del_likelihood[i];
		a_strvar->rp[i] = rp[i];
		a_strvar->sr[i] = sr[i];
	}

	return a_strvar;
}


void print_sv_stats()
{
	fprintf(logFile,"\n\nTARDIS found %d SVs total\n", del_cnt + ins_cnt + inv_cnt + tandup_cnt + mei_cnt);
	fprintf(logFile,"\tDeletion: %d\t sv length = %li\n", del_cnt, total_del_length);
	fprintf(logFile,"\tInsertion: %d\t sv length = %li\n", ins_cnt, total_ins_length);
	fprintf(logFile,"\tInversion: %d\t sv length = %li\n", inv_cnt, total_inv_length);
	fprintf(logFile,"\tTandem Duplication: %d\t sv length = %li\n", tandup_cnt, total_tandup_length);
	fprintf(logFile,"\tMEI: %d (%d filtered)\t sv length = %li\n", mei_cnt, mei_cnt_filtered, total_mei_length);
}

//add the variation in ascending order according to inner_start
void print_strvar( bam_info** in_bams, parameters* params, struct strvar* sv, FILE* fpOut)
{
	int sv_len, i, j, control, ind_id, rp_total = 0, sr_total = 0;

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

	sv_len = abs( ( sv->inner_end - sv->inner_start + 1));
	if( sv->svtype == DELETION)
	{
		if(sv->mei_del == true)
		{
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, ".", "<", "DEL:ME:", sv->mei_name, ">",( sv->filtered == false) ? "PASS" : "LowQual");
		}
		else
			fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t255\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_del_", ++del_cnt, ".", "<", "DEL", ">",( sv->filtered == false) ? "PASS" : "LowQual");
		total_del_length += sv_len;
	}

	else if( sv->svtype == INSERTION)
	{
		fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_ins_", ++ins_cnt, ".", "<", "INS", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		total_ins_length += sv_len;
	}
	else if( sv->svtype == INVERSION)
	{
		fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_inv_", ++inv_cnt, ".", "<", "INV", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		total_inv_length += sv_len;
	}
	else if( sv->svtype == MEIFORWARD || sv->svtype == MEIREVERSE)
	{
		fprintf( fpOut, "%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1, "vh_mei_", ++mei_cnt, ".", "<", "MEI", ">", 255, ( sv->filtered == false) ? "PASS" : "mfilt");
		if( sv->filtered)
		{
			mei_cnt_filtered++;
			total_mei_length += sv_len;
		}
	}
	else if( sv->svtype == TANDEMDUP)
	{
		fprintf( fpOut,"%s\t%i\t%s%d\t%s\t%s%s%s\t%d\t%s\t", sv->chr_name, ( sv->inner_start) + 1," vh_dup_", ++tandup_cnt, ".", "<", "DUP:TANDEM", ">", 255, ( sv->filtered == false) ? "PASS" : "LowQual");
		total_tandup_length += sv_len;
	}

	if( sv->svtype == MEIFORWARD || sv->svtype == MEIREVERSE)
		fprintf( fpOut, "END=%d;SVLEN=%d;TYPE=%s;RPSUP=%d;SRSUP=%d;", (sv->inner_end) + 1, sv_len, sv->mei_name, rp_total, sr_total);
	else
		fprintf( fpOut, "END=%d;SVLEN=%d;RPSUP=%d;SRSUP=%d;", (sv->inner_end) + 1, sv_len, rp_total, sr_total);

	/* If SV is imprecise */
	if( sv->inner_start != sv->inner_end)
		fprintf( fpOut, "CIEND=0,%d;CIPOS=-%d,0;IMPRECISE;", ( sv->outer_end - sv->inner_end), ( sv->inner_start - sv->outer_start));

	if( sv->svtype == DELETION)
		fprintf( fpOut, "SVTYPE=DEL\t");
	else if( sv->svtype == MEIFORWARD || sv->svtype == MEIREVERSE)
		fprintf( fpOut, "SVTYPE=MEI\t");
	else if( sv->svtype == TANDEMDUP)
		fprintf( fpOut, "SVTYPE=DUP\t");
	else if( sv->svtype == INSERTION)
		fprintf( fpOut, "SVTYPE=INS\t");
	else if( sv->svtype == INVERSION)
		fprintf( fpOut, "SVTYPE=INV\t");

	/* Format field */
	if (params->ten_x || params->output_hs)
		fprintf( fpOut, "GT:DL:RD:CN:RP:SR:HS:WE");
	else
		fprintf( fpOut, "GT:DL:RD:CN:RP:SR");

	control = 0;
	fprintf( fpOut, "\t");
	for( j = 0; j < params->num_bams; j++)
	{
		if( in_bams[j]->contribution == false)
		{
			if (params->ten_x || params->output_hs)
				fprintf( fpOut, "0/0:%.1f:%li:%.1f:%d:%d:%8.6f:%8.10f\t", sv->del_likelihood[j], sv->depth[j], sv->cn[j], sv->rp[j], sv->sr[j], sv->homogeneity_score, sv->weight);
			else 
				fprintf( fpOut, "0/0:%.1f:%li:%.1f:%d:%d\t", sv->del_likelihood[j], sv->depth[j], sv->cn[j], sv->rp[j], sv->sr[j]);
		}
		else
		{
			if (params->ten_x || params->output_hs)
				fprintf( fpOut, "0/1:%.1f:%li:%.1f:%d:%d:%8.6f:%8.10f\t", sv->del_likelihood[j], sv->depth[j], sv->cn[j], sv->rp[j], sv->sr[j], sv->homogeneity_score, sv->weight);
			else
				fprintf( fpOut, "0/1:%.1f:%li:%.1f:%d:%d\t", sv->del_likelihood[j], sv->depth[j], sv->cn[j], sv->rp[j], sv->sr[j]);
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
	char header_info[]="##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">\n"
			"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n"
			"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"
			"##INFO=<ID=DBRIPID,Number=1,Type=String,Description=\"ID of this element in DBRIP\">\n"
			"##INFO=<ID=DBVARID,Number=1,Type=String,Description=\"ID of this element in DBVAR\">\n"
			"##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample ID\">\n"
			"##INFO=<ID=SVALG,Number=1,Type=String,Description=\"SV discovery algorithm\">\n"
			"##INFO=<ID=RPSUP,Number=1,Type=Integer,Description=\"Number of supporting read pairs\">\n"
			"##INFO=<ID=SRSUP,Number=1,Type=Integer,Description=\"Number of supporting split reads\">\n"
			"##INFO=<ID=DGVID,Number=1,Type=String,Description=\"ID of this element in Database of Genomic Variation\">\n"
			"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n"
			"##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">\n"
			"##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">\n"
			"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
			"##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME\">\n"
			"##INFO=<ID=METRANS,Number=4,Type=String,Description=\"Mobile element transduction info of the form CHR\">\n"
			"##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=\"Indicates a novel structural variation\">\n"
			"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
			"##INFO=<ID=SVMETHOD,Number=.,Type=String,Description=\"Type of approach used to detect SV: RP (read pair), RD (read depth), SR (split read), or AS (assembly)\">\n"
			"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
			"##INFO=<ID=VALIDATED,Number=0,Type=Flag,Description=\"Validated by PCR, Assembly or Microarray\">\n"
			"##INFO=<ID=VALMETHOD,Number=.,Type=String,Description=\"Type of validation: CGH, PCR, SAV (superarray), CAP (capture-array), or ASM (assembly)\">\n";

	char header_filter[]="##FILTER=<ID=LowQual,Description=\"Genotype call confidence below LOD 1.3\">\n"
			"##FILTER=<ID=dpr5,Description=\"Read Depth probability below 5%\">\n";

	char header_format[] = "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n"
			"##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">\n"
			"##FORMAT=<ID=DL,Number=1,Type=String,Description=\"Deletion Likelihood\">\n"
			"##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Per-sample genotype filter, PASS for called genotypes or list of excluding filters\">\n"
			"##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Genotype Likelihoods\">\n"
			"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=HS,Number=1,Type=Float,Description=\"10x Barcode Homogeneity Score\">\n"
			"##FORMAT=<ID=RD,Number=1,Type=String,Description=\"Read Depth\">\n"
			"##FORMAT=<ID=RP,Number=1,Type=String,Description=\"Read Pair Support\">\n"
			"##FORMAT=<ID=SR,Number=1,Type=String,Description=\"Split Read Support\">\n";

	char header_alt[]="##ALT=<ID=DEL,Description=\"Deletion\">\n";
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	fprintf(fpOut,"##fileformat=VCFv4.1\n");
	fprintf(fpOut, "##fileDate=%d%s%d%s%d\n", timeinfo->tm_year+1900, (timeinfo->tm_mon+1<10 ? "0" : ""), timeinfo->tm_mon+1, (timeinfo->tm_mday<10 ? "0" : ""), timeinfo->tm_mday);
	/* TODO. Fix this with the SONIC info field */
	fprintf(fpOut,"##reference=1000GenomesPilot-NCBI37\n");
	fprintf(fpOut, "%s%s%s%s", header_info,header_filter,header_format,header_alt);
	fprintf(fpOut, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s","#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT");

	for( i = 0; i < params->num_bams; i++)
	{
		fprintf( fpOut, "\t%s",in_bams[i]->sample_name);
	}
	fprintf( fpOut, "\n");
}


void free_sensitive(bam_info** in_bams, parameters *params, ref_genome* ref)
{
	int lib_index, i, j;

	if ( params->run_rd == 1)
	{
		/* Free the RD array */
		for( i = 0; i < params->num_bams; i++)
		{
			for( lib_index = 0; lib_index < in_bams[i]->num_libraries; lib_index++)
			{
				for( j = 0; j < params->this_sonic->number_of_chromosomes; j++)
					free( in_bams[i]->read_depth[j]);
				free( in_bams[i]->read_depth);
			}
		}
	}

	for( i = 0; i < params->num_bams; i++)
	{
		free( in_bams[i]->sample_name);
		for( j = 0; j < in_bams[i]->num_libraries; j++)
		{
			free( in_bams[i]->libraries[j]->divet);
			free( in_bams[i]->libraries[j]->fastq1);
			free( in_bams[i]->libraries[j]->fastq2);
			free( in_bams[i]->libraries[j]->libname);
			free( in_bams[i]->libraries[j]->listSoftClip);
		}
		free( in_bams[i]);
	}
	free( in_bams);

	free( params->bam_files);
	free( params->bam_list_path);
	free( params->dups);
	free( params->gaps);
	free( params->mei);
	free( params->outprefix);
	free( params->ref_genome);
	free( params->reps);
	free( params->sonic_file);
	free( params);
}

void free_quick(bam_info** in_bams, parameters *params, ref_genome* ref)
{
	int lib_index, i, j;

	/* Free refgenome struct */
	free( ref->ref_name);
	free( ref->in_bam);
	free( ref);

	/* Free bams and related libraries */
	for( i = 0; i < params->num_bams; i++)
	{
		free( in_bams[i]->sample_name);
		for( j = 0; j < in_bams[i]->num_libraries; j++)
		{
			free( in_bams[i]->libraries[j]->divet);
			free( in_bams[i]->libraries[j]->fastq1);
			free( in_bams[i]->libraries[j]->fastq2);
			free( in_bams[i]->libraries[j]->libname);
			free( in_bams[i]->libraries[j]->listSoftClip);
		}
		free( in_bams[i]);
	}
	free( in_bams);

	/* Free params struct */
	free( params->bam_files);
	free( params->bam_list_path);
	free( params->dups);
	free( params->gaps);
	free( params->mei);
	free( params->outprefix);
	free( params->reps);
	free( params->sonic_file);
	free( params);
}


