#include "src/stdafx.h"

void CopyFileA(const char* from, const char* to, bool isBreakIfExists=true) {
    bool isFoutExist = (access(to,F_OK)==0);
    if(isFoutExist && isBreakIfExists) {
        std::cout << "[DEBUG] file [" << to << "] exist, skip copy..." << std::endl;
    } else {
        std::ifstream fFrom(from, std::ios::binary);
        std::ofstream fTo  (to,   std::ios::binary|std::ios::trunc);
        fTo << fFrom.rdbuf();
    }
}

//////////////////////////////////////////////////
int main(int argc, char* argv[]) {

	LOG_Begin();
	srand(time(NULL));
	LOG("v1.0");

	/*argc = 3;
	argv[1] = "id005.hdr";
	argv[2] = "-s";*/

	char * image_fn, * shape_fn;
	bool writefile = true;
    char defaultSurface[] = "sphere3.obj";

	if (argc < 2) {
		printf("<< CT CHEST SEGMENTER >>\n");
		printf("A tool for coarse segmentation of pathological lungs on CT [Liauchuk V., 2014]\n\n");
		printf("Usage: CTChestSegmenter.exe <analyze_image.hdr> [surface.obj] [-s]\n");
		printf("       if -s is specified then the segmentation result is not saved but\n");
		printf("       displayed only\n");
		return 0;
	}
	else{
		image_fn = argv[1];
	}

	if (argc < 3)
        shape_fn = defaultSurface; //"sphere3.obj";
	else {
        writefile = strcmp(argv[2], "-s") != 0;
        shape_fn = writefile ? argv[2] : &defaultSurface[0];
	}

	if (writefile)
        writefile = (argc < 4) ? true : (strcmp(argv[3], "-s") != 0);
	
	char buf[10000];
	try {
		AlgorithmData * alg = getAlgorithmData();
		Analyze75Read(image_fn, &(alg->im));

		ReadObj(shape_fn, alg->sphere.vertices, alg->sphere.faces, alg->sphere.edges, 
			&alg->sphere.nv, &alg->sphere.nf, &alg->sphere.ne);
		sprintf(buf, "Surface read : %i vertices, %i faces, %i edges\n", alg->sphere.nv, alg->sphere.nf, alg->sphere.ne);
		LOG(buf);

		LOG("Preprocessing Image");
		PreprocessImage();

		alg->center = Point3f(0.6, 0.5, 0.45);

		LOG("Setting initial iteration...");
		Mat x0 = InitialIteration();
		sprintf(buf, "Cost (initial) = %f", CostFunction(x0));
		LOG(buf);

		LOG("Gradient descent minimization...");
		Mat x = GradientDescend(x0, 10000);
		sprintf(buf, "Cost (final) = %f", CostFunction(x));
		LOG(buf);

		Surface surf = MakeSegmentingShape(x, alg->im.hdrInfo.dims.dim);
		LOG("Making raster image");
		ShapeToSegmentedAnalyze(alg->im, surf.vertices, surf.faces, surf.nv, surf.nf);
		
		if (writefile){
			LOG("Saving the result");
			char fn[1000];
			strcpy(fn, image_fn);
			CropFilenameExtention(fn);
			sprintf(fn, "%s_chsegm.img", fn);
			FILE * f = fopen(fn, "wb");
			for (int k = 0; k < alg->im.slices.size(); k++){
				fwrite(alg->im.slices[k].data, sizeof(short), alg->im.slices[k].rows * alg->im.slices[k].cols, f);
			}
			fclose(f);

			MakeHDRFilename(image_fn, buf);
			strcpy(fn, image_fn);
			CropFilenameExtention(fn);
			sprintf(fn, "%s_chsegm.hdr", fn);
            CopyFileA(buf, fn, false);
		}
		else
			Analyze75Show(&alg->im, 1. / 2000.);
			
	}
	catch (char * err){
		char errbuf[2048];
		sprintf(errbuf, "*** ERROR : %s\n", err);
		LOG(errbuf);
		fcloseall();
		return 1;
	}
	catch (...){
		LOG("*** UNRECOGNIZED ERROR\n");
		fcloseall();
		return 2;
	}

	LOG("FINISHED");

	return 0;
}
