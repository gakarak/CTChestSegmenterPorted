
#include "stdafx.h"

#define LOGIFLE "chestsegmenter.log"

///////////////////////////////////////////////////////////

void LOG(const char * mess, bool onScreen, bool toFile){
	if (onScreen){
		printf("%s\n", mess);
	}

	if (toFile){
		time_t t = time(NULL);
		char s_time[30];
		strcpy(s_time, ctime(&t));
		s_time[strlen(ctime(&t)) - 1] = 0;

		FILE * f = fopen(LOGIFLE, "at");
		fprintf(f, "%s  --  %s\n", s_time, mess);
		fclose(f);
	}
}

void LOG_Begin(){
	time_t t = time(NULL);
	char s_time[30];
	strcpy(s_time, ctime(&t));
	s_time[strlen(ctime(&t)) - 1] = 0;

	FILE * f = fopen(LOGIFLE, "at");
	fprintf(f, "\n\n");
	fprintf(f, "*****************************************************\n");
	fprintf(f, "* SEEKING AGENTS                                    *\n");
	fprintf(f, "*****************************************************\n");
	fprintf(f, "%s  --  START\n", s_time);
	fclose(f);
}