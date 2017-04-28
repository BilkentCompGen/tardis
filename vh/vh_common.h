#ifndef VH_COMMON_H__
#define VH_COMMON_H__

#ifndef NULL
#define NULL 0
#endif

//This define lines selects the main function that should be run
//#define MAIN_HASH     //For testing the hash file
//#define MAIN_MAIN     //The main Main function of the program

//#define MAIN_DELETION_CLUSTER

#include <stdlib.h>
#include "vh_logger.h"
#include "../common.h"

#define EXIT_CODE_SUCCESS 0	//Exit code for exiting program successfully
#define EXIT_CODE_ARG_ERROR 1	//Exit code 1 is for error in main's argument list
#define EXIT_CODE_DIVET_ERROR 2	//2 is for any error regarding the input divet files
#define EXIT_CODE_MEMORY_ERROR 3	//3 is for memory allocation problems
#define EXIT_CODE_GAP_ERROR	4

char g_error_message[500];
extern int g_maxListBrkPointIntr;	//defined in vh_createMaxClusterDeletion.cpp

extern struct LibraryInfo *g_libInfo;	//Defined in vh_divethandler.cpp
extern FILE *fileOutput;	//Defined in vh_main.cpp

#define maxInversionLen 5000000
#define maxDeletionLen 1000000
#define maxDuplicationLen 1000000

void quitProgram (int exitCode);


#endif
