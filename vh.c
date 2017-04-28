#include <stdio.h>
#include "vh.h"

int run_vh( ref_genome* ref, parameters *params, bam_info ** in_bams){
	int i, j;
	double preProsPrune= 0.001;
	int overMapLimit=500;
	char divetfile[MAX_SEQ];
	char outputfile[MAX_SEQ];
	char outputread[MAX_SEQ];
	char svfile[MAX_SEQ] = {};

	/* Print all structural variations in .vcf format */
	sprintf( svfile, "%s.vcf", params->outprefix);
	FILE* fpVcf = safe_fopen( svfile, "w");

	sprintf( outputfile,"%s.clusters", params->outprefix);
	sprintf( outputread,"%s.name", params->outprefix);

	for( i = 0; i < params->num_bams; i++)
	{
		for ( j = 0; j < in_bams[i]->num_libraries; j++)
		{
			sprintf( divetfile, "%s-%s.sam_DIVET.vh",  in_bams[i]->sample_name, in_bams[i]->libraries[j]->libname);
			set_str( &( in_bams[i]->libraries[j]->divet), divetfile);
		}
	}
	print_vcf_header(fpVcf, in_bams, params);

	vh_logInfo( "Calculating maximal clusters.");
	if ( TARDIS_DEBUG && !params->skip_vhcluster) // this parameter is only intended for debugging purposes. End users shouldn't use this
		vh_clustering( in_bams, ref, params, preProsPrune, outputfile, outputread, overMapLimit);

	vh_logInfo( "Applying SET-COVER approximation to find putative structural variation.");
	vh_setcover( in_bams, params, ref, outputread, outputfile, fpVcf);

	fprintf( stderr, "\nTARDIS is complete. Found %d SVs. Results are in the %s file.", sv_count, svfile);

	if ( !TARDIS_DEBUG)
	{
		remove( outputfile);
		remove( outputread);
	}
	return RETURN_SUCCESS;
}
