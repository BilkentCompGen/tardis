#ifndef __COMMANDLINE
#define __COMMANDLINE

#include "common.h"

int parse_command_line( int, char**, parameters*);
int parse_bam_list( parameters** params);
void print_help( void);

#endif
