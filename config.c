#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"

void load_config( configuration* cfg)
{
	FILE* config;
	char* next_line = NULL;
	char config_filename[MAX_LENGTH];
	size_t len = 0;
	int bytes_read;
	int i;

	/* Initialize the executable paths to empty */
	cfg->path_samtools = NULL;
	cfg->path_bcftools = NULL;
	cfg->path_mrfast = NULL;
	cfg->path_gnuplot = NULL;
	cfg->path_megablast = NULL;

	/* Combine the home directory path with the name of the configuration file */
	sprintf( config_filename, "%s/%s", getenv( "HOME"), CONFIG_FILE);

	/* Open the configuration file for reading. Do not use safe_fopen here, we will create if it doesn't exist */
	config = fopen( config_filename, "r");
	if( config == NULL)
	{
		/* Create new config file */
		create_config( cfg, config_filename);
	}
	else
	{
		/* Get the paths from the pre-existing configuration file */
		i = 0;
		while( ( bytes_read = getline( &next_line, &len, config)) != -1)
		{
			/* Ignore comments, which begin with the '#' character */
			if( next_line[0] != '#')
			{
				next_line[strlen(next_line)-1] = 0; // get rid of \n

				if( strstr( next_line, "SAMTOOLS"))
				{
					set_str( &( cfg->path_samtools), next_line + strlen("SAMTOOLS = "));
				}
				else if( strstr( next_line, "BCFTOOLS"))
				{
					set_str( &( cfg->path_bcftools), next_line + strlen("BCFTOOLS = "));
				}
				else if( strstr( next_line, "MRFAST"))
				{
					set_str( &( cfg->path_mrfast), next_line + strlen("MRFAST = "));
				}
				else if( strstr( next_line, "GNUPLOT"))
				{
					set_str( &( cfg->path_gnuplot), next_line + strlen("GNUPLOT = "));
				}				
				else if( strstr( next_line, "MEGABLAST"))
				{
					set_str( &( cfg->path_megablast), next_line + strlen("MEGABLAST = "));
				}
				else
				{
					fprintf( stderr, "Configuration file has wrong format or unspecified external tools used.\n");
				}
				i = i + 1;
			} 
		}

		/* Free memory allocated internally by getline */
		free( next_line);

		/* If the paths are still NULL, then they are either not installed, or not in the PATH. */

		/*
		if( cfg->path_samtools == NULL)
		{
			fprintf( stderr, "Warning: samtools path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "samtools path: %s\n", cfg->path_samtools);
		}
		 */

		/*
		if( cfg->path_bcftools == NULL)
		{
			fprintf( stderr, "Warning: bcftools path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "bcftools path: %s\n", cfg->path_bcftools);
		}
		 */

		if( cfg->path_mrfast == NULL)
		{
			fprintf( stderr, "Warning: mrfast path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "mrfast path: %s\n", cfg->path_mrfast);
		}

		if( cfg->path_gnuplot == NULL)
		{
			fprintf( stderr, "Warning: gnuplot path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "gnuplot path: %s\n", cfg->path_gnuplot);
		}

		/*
		if( cfg->path_megablast == NULL)
		{
			fprintf( stderr, "Warning: megablast path is not in the configuration file.\n");
		}
		else
		{
			fprintf( stderr, "megablast path: %s\n", cfg->path_megablast);
		}
		 */
	}
}

void create_config( configuration* cfg, char* config_filename)
{
	FILE* config;
	FILE* pipe;
	char executable_path[MAX_LENGTH];

	/* popen(...) executes shell commands */
	/* "which" finds the path to the specified executable if it exists */
	/* "2>/dev/null" redirects the output into a virtual black hole */
	pipe = popen( "which samtools 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "samtools not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			set_str( &( cfg->path_samtools), executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which bcftools 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "bcftools not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			set_str( &( cfg->path_bcftools), executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which mrfast 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "mrfast not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			set_str( &( cfg->path_mrfast), executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which gnuplot 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "gnuplot not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			set_str( &( cfg->path_gnuplot), executable_path);
		}
		pclose( pipe);
	}

	pipe = popen( "which megablast 2>/dev/null", "r");
	if( pipe == NULL)
	{
		fprintf( stderr, "Error opening pipe\n");
	}
	else
	{
		if( fgets( executable_path, MAX_LENGTH, pipe) == NULL)
		{
			fprintf( stderr, "megablast not found in PATH. Install it or manually configure the %s file.\n", config_filename);
		}
		else
		{
			set_str( &( cfg->path_megablast), executable_path);
		}
		pclose( pipe);
	}

	config = safe_fopen( config_filename, "w");
	if( cfg->path_samtools != NULL)
	{
		fprintf( config, "SAMTOOLS = %s", cfg->path_samtools);
	}

	if( cfg->path_bcftools != NULL)	
	{
		fprintf( config, "BCFTOOLS = %s", cfg->path_bcftools);
	}

	if( cfg->path_mrfast != NULL)
	{
		fprintf( config, "MRFAST = %s", cfg->path_mrfast);
	}

	if( cfg->path_gnuplot != NULL)
	{
		fprintf( config, "GNUPLOT = %s", cfg->path_gnuplot);
	}

	if( cfg->path_megablast != NULL)
	{
		fprintf( config, "MEGABLAST = %s", cfg->path_megablast);
	}
	fclose( config);

	fprintf( stderr,"\n");
	fprintf( stderr, "*************************************************************************\n");
	fprintf( stderr, "*           Config file $HOME/.tardis_config is created.                *\n"); 
	fprintf( stderr, "*            Check the configuration file for any errors,               *\n");
	fprintf( stderr, "* then run TARDIS again to load the configuration and normal operation. *\n");
	fprintf( stderr, "*************************************************************************\n");
}
