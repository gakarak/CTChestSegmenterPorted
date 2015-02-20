#include "src/stdafx.h"

#include <nifti1_io.h>

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
//        Analyze75Read(image_fn, &(alg->im));
        nifti_image * nim = NULL;
        nim = nifti_image_read(image_fn, 1);
        if( !nim ) {
           std::stringstream ss;
           ss << "** Failed to read NIfTI image from file [" << image_fn << "]";
           throw ss.str();
        }

        if(nim->datatype!=4) {
            nifti_image_free(nim);
            std::stringstream ss;
            ss << "Invalid data format for 3D-image [" << image_fn << "]";
            throw ss.str();
        }
////////////////////////////////
        short int *ptrData = (short int*)nim->data;
        int nX, nY, nZ;
        nX = nim->nx;
        nY = nim->ny;
        nZ = nim->nz;
        short int siScale     = static_cast< short int >( nim->scl_slope );
        short int siIntercept = static_cast< short int >( nim->scl_inter );
        size_t numVoxels = nX*nY*nZ;
        if( abs( siScale ) > 0 ) {
            for( size_t i = 0; i < numVoxels; i++ ) {
                ptrData[i] *= siScale;
            }
        }
        if( abs( siIntercept ) > 0 ) {
            for( size_t i = 0; i < numVoxels; i++ ) {
                ptrData[ i ] += siIntercept;
            }
        }
////////////////////////////////
        HDRInfo hdrInfo;
        int numDimNifti = ARRAY_LEN(nim->dim);
        int numDimAnalz = ARRAY_LEN(hdrInfo.dims.dim);
        int numMin = cv::min(numDimAnalz,numDimNifti);
        for(size_t ii=0; ii<numMin-1; ii++) {
            hdrInfo.dims.dim[ii]=nim->dim[ii];
        }
        hdrInfo.dims.datatype = (short)nim->datatype;
        alg->im.hdrInfo = hdrInfo;
        int numZ = nim->dim[3];
        int sliceShift = nim->dim[1]*nim->dim[2];
        alg->im.slices = Vector<cv::Mat>(nim->dim[3]);
        for(int ii=0; ii<numZ; ii++) {
            Mat m(nim->dim[1], nim->dim[2], CV_16S, &ptrData[ii*sliceShift]); //, sizeof(short) * nim->dim[2]);
            m.copyTo(alg->im.slices[ii]);
            /*
            double vMin, vMax;
            cv::minMaxLoc(m, &vMin, &vMax, NULL, NULL);
            printf("[%d/%d] min/max = (%0.2f/%0.2f)\n", ii, numZ, vMin, vMax);
            */
        }

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
		    std::stringstream ssfout;
            ssfout << image_fn;
            ssfout << "_segm.nii.gz";
            for (int k = 0; k < alg->im.slices.size(); k++) {
                short* tmpPtr = (short*)alg->im.slices[k].data;
                std::copy(  tmpPtr,
                            tmpPtr+sliceShift,
                            &ptrData[k*sliceShift]);
            }
            if( nifti_set_filenames(nim, ssfout.str().c_str(), 1, 1) ) {
                nifti_image_free(nim);
                std::stringstream ss;
                ss << "** Failed to save 3D-image  [" << ssfout.str() << "]";
                throw ss.str();
            }
            nifti_image_write( nim );
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
