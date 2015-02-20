/********************************************************
*
* Some non-specific useful functions
*
**********************************************************/

#include "stdafx.h"

char str_buffer[1024];

double PI(){
	return 3.1415926535;
}

char * buf(){
	return str_buffer;
}

float sqr(float x){
	return x*x;
}

float dist2(float x1, float y1, float x2, float y2){
	return sqrt(sqr(x1 - x2) + sqr(y1 - y2));
}

float dist3(float x1, float y1, float z1, float x2, float y2, float z2){
	return sqrt(sqr(x1 - x2) + sqr(y1 - y2) + sqr(z1 - z2));
}

int signum(int x){
	return (x > 0) ? 1 : (x < 0) ? -1 : 0;
}

void order2intAsc(int* p1, int* p2){
	if ((*p1) > (*p2)) {
		int tmp = (*p2);
		(*p2) = (*p1);
		(*p1) = tmp;
	}
}

void ClearFile(const char * filename){
	FILE * f = fopen(filename, "wt");
	fclose(f);
}

void CropFilenameExtention(char* s){
	for (int i = strlen(s) - 1; i > 0; i--){
		if (s[i] == '.'){
			s[i] = 0;
			return;
		}
	}
}


Mat dlmread(const char * filename, char delim){

	CvMLData ml;

	if (access(filename, 0) == -1){
		throw "Specified file not found";
	}

	ml.set_delimiter(delim);
	ml.read_csv(filename);

	Mat m(ml.get_values(), true);

	return m; 
}

Mat colread(const char * filename){

	if (access(filename, 0) == -1){
		throw "Specified file not found";
	}

	char buf[100];

	FILE * f = fopen(filename, "rt");

	int numel = 0;
	while (!feof(f)){
		buf[0] = 0;
		fgets(buf, 100, f);
		if (strlen(buf) > 1)
			numel++;
	}

	Mat m(numel, 1, CV_32F, Scalar(0));

	fseek(f, 0, SEEK_SET);
	numel = 0;
	while (!feof(f)){
		buf[0] = 0;
		fgets(buf, 100, f);
		if (strlen(buf) > 1){
			float fl;
			sscanf(buf, "%f\n", &fl);
			m.at<float>(numel, 0) = fl;
			numel++;
		}	
	}

	fclose(f);

	return m;
}

void dlmwrite(const char * filename, Mat &m, char delim){
	FILE * f = fopen(filename, "wt");

	if (f == NULL){
		throw "Failed to append file for writing";
	}

	for (int r = 0; r < m.rows; r++){
		Mat toSave;
		m.row(r).convertTo(toSave, CV_32F);

		MatIterator_<float> it, end; 
		for (it = toSave.begin<float>(), end = toSave.end<float>(); it != end; ++it){
			if (it != end - 1)
				fprintf(f, "%f%c", *it, delim);
			else
				fprintf(f, "%f\n", *it);
		}
	}

	fclose(f);
}