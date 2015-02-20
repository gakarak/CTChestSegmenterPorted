/********************************************************
*
* Functions to handle Analyze75 file format
*
*********************************************************/

#include "../stdafx.h"

using namespace cv;


void Analyze75ShowVector(Vector<Mat> &v, short * dims, double alpha_, double beta_);


bool NoAnalyzeExt(const char * filename){
	bool noExt = true;
	const char * pch;

	pch = strstr(filename, ".img");
	noExt = noExt && (pch == NULL || pch != (filename + strlen(filename) - 4));

	pch = strstr(filename, ".hdr");
	noExt = noExt && (pch == NULL || pch != (filename + strlen(filename) - 4));

	return noExt;
}

void MakeExtFilename(const char * filename, char * extFilename, const char * ext){
	if (NoAnalyzeExt(filename))
        sprintf(extFilename, "%s.%s\n", filename, ext);
	else{
		strncpy(extFilename, filename, strlen(filename) - 4);
        extFilename[strlen(filename) - 4] = '\n';
		strcat(extFilename, ".");
		strcat(extFilename, ext);
	}
}

void MakeHDRFilename(const char * filename, char * hdrFilename){
	MakeExtFilename(filename, hdrFilename, "hdr");
}

void MakeIMGFilename(const char * filename, char * imgFilename){
	MakeExtFilename(filename, imgFilename, "img");
}

HDRInfo Analyze75Info(const char * filename)
{
	char hdrFilename[256];
	MakeExtFilename(filename, hdrFilename, "hdr");

	FILE* pFile = fopen(hdrFilename, "rb");

	if (pFile == NULL){
		throw "Specified file not found";
	}

	HDRInfo hdrInfo;

	fread(&hdrInfo, 348, 1, pFile);
	fclose(pFile);

	return hdrInfo;
}

void Analyze75Read(const char * filename, AnalyzeImage * im){
	char str[1024], buff[1024];

	sprintf(buff, "Reading analyze '%s'... ", filename);
	LOG(buff);

	MakeHDRFilename(filename, str);
	HDRInfo info = Analyze75Info(str);
	short d2 = info.dims.dim[2];
	info.dims.dim[2] = info.dims.dim[1];
	info.dims.dim[1] = d2;
	im->hdrInfo = info;

	if (info.dims.datatype != DT_SIGNED_SHORT)
		throw "Analyze data type is not DT_SIGNED_SHORT";

	short * dims = info.dims.dim;
	im->slices = Vector<Mat>(dims[3]);

	MakeIMGFilename(filename, str);
	FILE * f = fopen(str, "rb");
	if (f == NULL){
		sprintf(buf(), "File not found '%s'", str);
		throw buf();
	}

	short sbuf[262144];
	for (int k = 0; k < dims[3]; k++){
		fread(sbuf, sizeof(short), dims[1] * dims[2], f);

        Mat m(dims[1], dims[2], CV_16S, sbuf, sizeof(short) * dims[2]);
		m.copyTo(im->slices[k]);
		
		/*im->slices[k] = Mat(dims[1], dims[2], CV_16S);
		for (int ii = 0; ii < dims[1] * dims[2]; ii++){
			int i = ii / dims[1];
			int j = ii % dims[2];

			im->slices[k].at<short>(i, j) = sbuf[ii];
		}*/
	}

	fclose(f);

	LOG("OK");
}

void Analyze75Show(AnalyzeImage * im, double alpha_, double beta_){
	short * dims = im->hdrInfo.dims.dim;
	Analyze75ShowVector(im->slices, dims, alpha_, beta_);
}

void Analyze75ShowVector(Vector<Mat> &v, short * dims, double alpha_, double beta_){
	Mat m;

	int i = dims[3] / 2, key = -1;
	while (key != 27){
//		v[i].convertTo(m, CV_32F, alpha_, beta_);
        double vMin, vMax;
        cv::minMaxLoc(v[i],&vMin,&vMax,NULL,NULL);
        cv::normalize(v[i],m, 0,255, CV_MINMAX, CV_8U);
		sprintf(buf(), "%ix%ix%i", dims[1], dims[2], dims[3]);
        std::stringstream ss;
        ss << "slice #" << i << "/" << dims[3] << ", min/max = (" << vMin << "/" << vMax << ")";
        cv::putText(m, ss.str(), cv::Point(20,20), CV_FONT_HERSHEY_PLAIN, 1.2, cv::Scalar::all(255), 3);
        cv::putText(m, ss.str(), cv::Point(20,20), CV_FONT_HERSHEY_PLAIN, 1.2, cv::Scalar::all(0  ), 1);
        imshow(buf(), m);

		key = waitKey(20);

        if (key == 'q') {
			i = (i - 1 + dims[3]) % dims[3];
//            std::cout << "key = " << key << ", slice #" << i << std::endl;
        }
        if (key == 'w') {
			i = (i + 1) % dims[3];
//            std::cout << "key = " << key << ", slice #" << i  << std::endl;
        }
		if (key == 27){
			destroyWindow(buf());
			return;
		}
	}
}
