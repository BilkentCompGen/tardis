#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "external.h"
#include "processfq.h"

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
		}
	}  

	return RETURN_SUCCESS;
}


void plot_histogram  ( parameters *params, configuration *cfg, char *sample_name, char *libname, int sample_size, int *fragment_lengths, float average, float stdev)
{
  FILE *hist;
  char histogram_file_name[MAX_LENGTH];
  char plot_file_name[MAX_LENGTH];
  int i, current;
  int cnt;
  char cmdline[4096];
  int return_value;
  
  if (cfg->path_gnuplot == NULL){
    fprintf(stderr, "[TARDIS EXTERNAL] Cannot plot because gnuplot is missing.\n");
    return;
  }

  sprintf(histogram_file_name, "%s%s-%s.hist", params->outprefix, sample_name, libname);
  hist = fopen(histogram_file_name, "w");

  current = -1819289; /* magic number */
  cnt = 1;

  for (i = 0; i < sample_size; i++){
    if (fragment_lengths[i] == current)
      cnt++;
    else{
      if (current != -1819289)
	fprintf( hist, "%d\t%d\n", current, cnt);
      cnt = 1;
      current = fragment_lengths[i];
    }
  }
  fclose(hist);

  sprintf(plot_file_name, "%s%s-%s.plot", params->outprefix, sample_name, libname);
  hist = fopen(plot_file_name, "w");
  
  fprintf(hist, "set xlabel \"Fragment size\"\n");
  fprintf(hist, "set ylabel \"Number of pairs\"\n");
  fprintf(hist, "set title \'Sample: %s - Library: %s - Average: %.2f - Stdev: %.2f\' font \"Arial,10\"\n", sample_name, libname, average, stdev);

  fprintf(hist, "set xrange [0:%d]\n", (int)(average+3*stdev));
  fprintf(hist, "set terminal pdf enhanced color\n");
  fprintf(hist, "set output \"%s-%s-%s.pdf\"\n", params->outprefix, sample_name, libname);
  
  fprintf(hist, "plot \"%s\" with boxes ti \"fragment\"\n", histogram_file_name);

  fclose(hist);

  sprintf(cmdline, "%s %s", cfg->path_gnuplot, plot_file_name);

      
  if ( TARDIS_DEBUG == 1)
    {
      fprintf( stderr, "[GNUPLOT COMMAND LINE] %s\n", cmdline);
    }
  
  return_value = system( cmdline);
  
  if (WIFSIGNALED(return_value) && (WTERMSIG(return_value) == SIGINT || WTERMSIG(return_value) == SIGQUIT))
    {
      fprintf( stderr, "Histogram plotting failed.\n");
      return;
    }

  remove(plot_file_name);
  remove(histogram_file_name);

}
