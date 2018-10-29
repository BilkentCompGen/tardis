#ifndef __COMMANDLINE
#define __COMMANDLINE

#include "common.h"

int parse_command_line( int, char**, parameters*);
int parse_bam_list( parameters** params);
int parse_bam_array( parameters* params, char *optarg);
void print_help( void);

#endif
