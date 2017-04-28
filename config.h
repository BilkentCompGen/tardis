#ifndef __CONFIG
#define __CONFIG

/* Name of the configuration file can be changed here */
#define CONFIG_FILE ".tardis_config"


#include "common.h"

/* External tool executable paths */
typedef struct _configuration
{
	char* path_samtools;
	char* path_bcftools;
	char* path_mrfast;
	char* path_gnuplot;
	char* path_megablast;
} configuration;

/* Function Prototypes */
void load_config( configuration* cfg);
void create_config( configuration* cfg, char* config_filename);

#endif
