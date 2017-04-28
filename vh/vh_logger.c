#include "vh_logger.h"

#include <stdio.h>

char tempMsg[500];
FILE *g_logOutputFile=NULL;
static void init_streams(void) __attribute__((constructor)); //constructor for initializing stdout

static void init_streams(void)
{
	g_logOutputFile = stdout;
}
int g_currentLogLevel = LOG_LEVEL_ALL;
char g_loggerMsgBuffer[400];

void vh_log (char *message, int logLevel)
{
	if (!(logLevel & g_currentLogLevel))
		return;

	fprintf (g_logOutputFile, "%s", message);
}

void vh_initLogger (FILE * logOutputFile, int logLevel)
{
	g_logOutputFile = logOutputFile;
	g_currentLogLevel = logLevel;
}

void vh_logError (char *message)
{
	if (!(LOG_LEVEL_LOG_ERROR & g_currentLogLevel))
		return;

	sprintf (tempMsg, "****ERROR: \t%s\n", message);
	vh_log (tempMsg, LOG_LEVEL_LOG_ERROR);
}

void vh_logOutput (char *message)
{
	vh_log (message, LOG_LEVEL_LOG_OUTPUT);
}

void vh_logInfo (char *message)
{
	if (!(LOG_LEVEL_LOG_INFO & g_currentLogLevel))
		return;

	sprintf (tempMsg, "INFO: \t%s\n", message);
	vh_log (tempMsg, LOG_LEVEL_LOG_INFO);
}

void vh_logWarning (char *message)
{
	if (!(LOG_LEVEL_LOG_WARNING & g_currentLogLevel))
		return;

	sprintf (tempMsg, "WARN: \t%s\n", message);
	vh_log (tempMsg, LOG_LEVEL_LOG_WARNING);
}

void vh_logDebug (char *message)
{
	if (!(LOG_LEVEL_DEBUG_INFO & g_currentLogLevel))
		return;

	sprintf (tempMsg, "DEBUG \t%s\n", message);
	vh_log (tempMsg, LOG_LEVEL_DEBUG_INFO);
}

void vh_logTime ()
{
	vh_log ("Log time to be implemented", LOG_LEVEL_LOG_WARNING);
}
