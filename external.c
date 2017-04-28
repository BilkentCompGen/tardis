#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "external.h"

int remap_mrfast( parameters *params, bam_info ** in_bams, configuration *cfg)
{
	int i, j;
	char cmdline[4096];
	int return_value;

	fprintf( stderr, "\nRemapping with mrFAST. Allons-y!\n\n");
	for( i = 0; i < params->num_bams; i++)
	{
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
			fprintf( stderr, "\nRemapping:\n\t\tSample: %s\n\t\tLibrary: %s\n", in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			if ( params->threads == 1)
			{
				sprintf( cmdline, "%s --search %s --pe --min %d --max %d --sample %s --lib %s --rg %s --seq1 %s --seq2 %s -o %s-%s.sam -u %s-%s.unmapped.fastq --seqcomp --outcomp",
						cfg->path_mrfast, params->ref_genome, in_bams[i]->libraries[j]->conc_min, in_bams[i]->libraries[j]->conc_max,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname, in_bams[i]->libraries[j]->libname,
						in_bams[i]->libraries[j]->fastq1, in_bams[i]->libraries[j]->fastq2,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			}
			else
			{
				sprintf( cmdline, "%s --search %s --pe --min %d --max %d --sample %s --lib %s --rg %s --seq1 %s --seq2 %s -o %s-%s.sam -u %s-%s.unmapped.fastq --seqcomp --outcomp --threads %d",
						cfg->path_mrfast, params->ref_genome, in_bams[i]->libraries[j]->conc_min, in_bams[i]->libraries[j]->conc_max,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname, in_bams[i]->libraries[j]->libname,
						in_bams[i]->libraries[j]->fastq1, in_bams[i]->libraries[j]->fastq2,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname,
						in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname, params->threads);
			}

			if ( TARDIS_DEBUG == 1)
			{
				fprintf( stderr, "[MRFAST COMMAND LINE] %s\n", cmdline);
			}

			return_value = system( cmdline);

			if (WIFSIGNALED(return_value) && (WTERMSIG(return_value) == SIGINT || WTERMSIG(return_value) == SIGQUIT))
			{
				fprintf( stderr, "mrFAST remapping failed.\n");
				return RETURN_ERROR;
			}

			/* to be implemented.
				params->threads will be used for multithreading option of mrFAST
			 */
		}
	}  

	return RETURN_SUCCESS;
}
