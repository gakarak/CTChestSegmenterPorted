#include "stdafx.h"

AnalyzeImage FilterAnalyze(struct AnalyzeImage & im, Mat & flt, bool viaFourier = false);

////////////////////////////////////////////////////////////////////////////////////////////////////

AlgorithmData alg;

AlgorithmData * getAlgorithmData(){
	return &alg;
}

AnalyzeImage FilterAnalyzeConv(struct AnalyzeImage & im, Mat & flt) 
{
	AnalyzeImage filtered;
	filtered.slices.reserve(im.slices.size());

	filtered.hdrInfo = im.hdrInfo;
	int dk = flt.size[2] / 2;
	int dj = flt.size[1] / 2;
	int di = flt.size[0] / 2;
	for (int k = 0; k < im.slices.size(); k++){
        printf("%d / %lu\t\t\r", k + 1, im.slices.size());
		Mat m(im.slices[0].rows, im.slices[0].cols, CV_32F, Scalar::all(0));
		for (int j = 0; j < im.slices[0].cols; j++){
			for (int i = 0; i < im.slices[0].rows; i++){
				for (int kk = 0; kk < flt.size[2]; kk++){
					for (int jj = 0; jj < flt.size[1]; jj++){
						for (int ii = 0; ii < flt.size[0]; ii++){
							if (k + kk - dk >= 0 && k + kk - dk < im.slices.size() &&
								j + jj - di >= 0 && j + jj - dj < im.slices[0].cols &&
								i + ii - dj >= 0 && i + ii - di < im.slices[0].rows)
								m.at<float>(i, j) += flt.at<float>(ii, jj, kk) * im.slices[k + kk - dk].at<float>(i + ii - di, j + jj - dj);
						}
					}
				}
			}
		}
		filtered.slices.push_back(m);
	}

	return filtered;
}

AnalyzeImage FilterAnalyzeFourier(struct AnalyzeImage & im, Mat & flt) 
{
	Mat planes[] = {flt, flt * 0};
	Mat cflt;	// complex-type filter
	merge(planes, 2, cflt);

	Mat wi = AnalyzeFFT(im);
	int N3[] = {wi.size[0], wi.size[1], wi.size[2]};
	Mat wf = Mat3FFT(cflt, false, N3);

	// element-wise mutiplification
	for (int i = 0; i < wf.size[0]; i++)
	for (int j = 0; j < wf.size[1]; j++)
	for (int k = 0; k < wf.size[2]; k++){
        std::complex<float> c1 = wf.at<std::complex<float> >(i, j, k);
        std::complex<float> c2 = wi.at<std::complex<float> >(i, j, k);
        wi.at<std::complex<float> >(i, j, k) = c1 * c2;
	}

	// inverse transform
	Mat im2 = Mat3FFT(wi, true, N3);

	int m[3] = {flt.size[0] / 2, flt.size[1] / 2, flt.size[2] / 2};
	AnalyzeImage filtered;
	filtered.hdrInfo = im.hdrInfo;
	filtered.slices.reserve(im.slices.size());
	for (int k = 0; k < im.slices.size(); k++){
		Mat sl0(im2.size[0], im2.size[1], CV_32FC2, Scalar::all(0)), sl;
		if (k + m[2] < im2.size[2]){
			Range rng[3] = {Range::all(), Range::all(), Range(k + m[2], k + m[2] + 1)};
		
			Mat plns[2];
			im2(rng).copyTo(sl);
			sl.copySize(sl0);
			split(sl, plns);
			magnitude(plns[0], plns[1], plns[0]);

			plns[1] = plns[1] * 0;
			plns[0](Range(2 * m[0], plns[0].size[0]), Range(2 * m[1], plns[0].size[1])).
				copyTo(plns[1](Range(m[0], plns[0].size[0] - m[0]), Range(m[1], plns[0].size[1] - m[1])));

			filtered.slices.push_back(plns[1]);
		}
		else {
			filtered.slices.push_back(filtered.slices[0] * 0);
		}
	}
	return filtered;
}

AnalyzeImage FilterAnalyze(struct AnalyzeImage & im, Mat & flt, bool viaFourier){
	if (viaFourier){
		return FilterAnalyzeFourier(im, flt);
	}
	else {
		return FilterAnalyzeConv(im, flt);
	}
}

Mat MakeFilter3D(float rmm, float * pixdim) 
{
	float rh = (ceil(rmm / pixdim[1]));
	float rv = (ceil(rmm / pixdim[3]));

	// std = rmm / 2
	int mdims[] = {2 * int(rh) + 1, 2 * int(rh) + 1, 2 * int(rv) + 1};
	Mat blurFilt(3, mdims, CV_32F);
	for (int k = 0; k < blurFilt.size[2]; k++){
		for (int j = 0; j < blurFilt.size[1]; j++){
			for (int i = 0; i < blurFilt.size[0]; i++){
				int idx[] = {i, j, k};
				float x2 = sqr((i - rh) / rh) + sqr((j - rh) / rh) + sqr((k - rv) / rv);
				blurFilt.at<float>(idx) = x2 < 1 ? exp(- x2 * 4) : 0;
			}
		}
	}
	blurFilt /= sum(blurFilt)[0];
	return blurFilt;
}

AnalyzeImage BinarizeResizeImage(AnalyzeImage &im) 
{
	AnalyzeImage bin;
	bin.hdrInfo = im.hdrInfo;

	float thr = 1350.f;
	// binarizing, resizing
	for (int k = 1; k < im.slices.size(); k += 2){
		Mat m1, m2;
		im.slices[k - 1].convertTo(m1, CV_32F);
		im.slices[k].convertTo(m2, CV_32F);
		threshold(m1, m1, thr, 1., THRESH_BINARY);
		threshold(m2, m2, thr, 1., THRESH_BINARY);
		Mat m = m1 + m2;
		resize(m, m, Size(0, 0), 0.5, 0.5);
		threshold(m, m, 0.5, 1.0, THRESH_BINARY);
		bin.slices.push_back(m);
	}
	bin.hdrInfo.dims.dim[1] = bin.hdrInfo.dims.dim[1] / 2;
	bin.hdrInfo.dims.dim[2] = bin.hdrInfo.dims.dim[2] / 2;
	bin.hdrInfo.dims.dim[3] = bin.slices.size();
	bin.hdrInfo.dims.pixdim[1] = bin.hdrInfo.dims.pixdim[1] * 2;
	bin.hdrInfo.dims.pixdim[2] = bin.hdrInfo.dims.pixdim[2] * 2;
	bin.hdrInfo.dims.pixdim[3] = bin.hdrInfo.dims.pixdim[3] * 2;

	return bin;
}

AnalyzeImage TruncateImage(AnalyzeImage &im, float val) 
{
	AnalyzeImage bin;
	bin.hdrInfo = im.hdrInfo;

	for (int k = 0; k < im.slices.size(); k++){
		Mat m;
		im.slices[k].convertTo(m, CV_32F);
		threshold(m, m, val, val, THRESH_TRUNC);
		bin.slices.push_back(m);
	}
	return bin;
}

void CalcGradients() {
	AlgorithmData * alg = getAlgorithmData();

	float data0[3][3][3] = {{{-1, -2, -1}, {-2, -4, -2}, {-1, -2, -1}}, 
	{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
	float data1[3][3][3] = {{{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}, 
	{{-2, -4, -2}, {0, 0, 0}, {2, 4, 2}}, {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}}};
	float data2[3][3][3] = {{{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}, 
	{{-2, 0, 2}, {-4, 0, 4}, {-2, 0, 2}}, {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}}};
	int dims[] = {3, 3, 3};
	Mat sobelFilt0(3, dims, CV_32F, data0);
	Mat sobelFilt1(3, dims, CV_32F, data1);
	Mat sobelFilt2(3, dims, CV_32F, data2);
	sobelFilt0 = sobelFilt0 * (1.f / 32.);
	sobelFilt1 = sobelFilt1 * (1.f / 32.);
	sobelFilt2 = sobelFilt2 * (1.f / alg->z2xy / 32.);

	alg->gr[0] = FilterAnalyze(alg->blurred, sobelFilt0);
	alg->gr[1] = FilterAnalyze(alg->blurred, sobelFilt1);
	alg->gr[2] = FilterAnalyze(alg->blurred, sobelFilt2);
}

void PreprocessImage(){
	AlgorithmData * alg = getAlgorithmData();

	AnalyzeImage bin = BinarizeResizeImage(alg->im);

	float rmm = 15;
	float * pixdim = bin.hdrInfo.dims.pixdim;
	Mat blurFilt = MakeFilter3D(rmm, pixdim);

	alg->blurred = bin;
	alg->blurred = FilterAnalyze(bin, blurFilt, true);
	alg->blurred = TruncateImage(alg->blurred, 0.1f);
	alg->blurred = FilterAnalyze(alg->blurred, blurFilt, true);
	alg->z2xy = pixdim[3] / pixdim[1];

	CalcGradients();	
}

bool PointWithin(struct Point3D p, struct Point3D p0, struct Point3D p1, struct Point3D p2, struct Point3D p3){
	Point3D a = Point3D_(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
	Point3D b = Point3D_(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);
	Point3D c = Point3D_(p3.x - p0.x, p3.y - p0.y, p3.z - p0.z);
	Point3D d = Point3D_(p.x - p0.x, p.y - p0.y, p.z - p0.z);

	float abc = MixProd(a, b, c);

	if (abc * MixProd(d, b, c) < 0)
		return false;
	if (abc * MixProd(a, d, c) < 0) 
		return false;
	if (abc * MixProd(a, b, d) < 0) 
		return false;

	Point3D pp = Point3D_(p1.x - p.x, p1.y - p.y, p1.z - p.z);
	Point3D ab = Point3D_(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
	Point3D ac = Point3D_(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);

	if (MixProd(a, ab, ac) * MixProd(pp, ab, ac) < 0)
		return false;

	return true;
}

Point3i VerticeToCoordinates(Point3D v, float x, short * dims) 
{
	AlgorithmData * alg = getAlgorithmData();
	Point3f rel;
	rel.x = v.x * x * dims[1] / 2 + dims[1] * alg->center.x;
	rel.y = v.y * x * dims[1] / 2 + dims[2] * alg->center.y;
	rel.z = v.z * x * dims[1] / 2 / alg->z2xy + dims[3] * alg->center.z;
	Point3i crd(int(floor(0.5 + rel.x)), int(floor(0.5 + rel.y)), int(floor(0.5 + rel.z)));
	return crd;
}

float CostFunction(Mat x, Mat * drv, float alpha_tension){
	float defaultAlphaTension = 6.f;
	alpha_tension = alpha_tension < -0.5f ? defaultAlphaTension : alpha_tension;

	AlgorithmData * alg = getAlgorithmData();
	short * dims = alg->blurred.hdrInfo.dims.dim;

	float bones = 0.f;
	Mat bonesDrv(x.rows, x.cols, CV_32F, Scalar::all(0));
	for (int i = 0; i < x.cols; i++){
		Point3D v = alg->sphere.vertices[i];
		Point3i crd = VerticeToCoordinates(v, x.at<float>(0, i), dims);

		if (crd.x >= 0 && crd.x < dims[1] && crd.y >= 0 && 
			crd.y < dims[2] && crd.z >= 0 && crd.z < dims[3]){
			bones += alg->blurred.slices[crd.z].at<float>(crd.x, crd.y);
		}

		if (drv){
			if (crd.x >= 0 && crd.x < dims[1] && crd.y >= 0 && 
				crd.y < dims[2] && crd.z >= 0 && crd.z < dims[3]){
					float tmp1 = dims[1] * (alg->gr[0].slices[crd.z].at<float>(crd.x, crd.y) * v.x + 
						alg->gr[1].slices[crd.z].at<float>(crd.x, crd.y) * v.y + 
						alg->gr[2].slices[crd.z].at<float>(crd.x, crd.y) * v.z);
					//printf("bdrv(%i) = %.5f\n", i, tmp1);
					bonesDrv.at<float>(0, i) = tmp1;
			}		
		}
	}
	bones = bones / x.cols;
	bonesDrv = bonesDrv / x.cols;

	if (alpha_tension == 0.f){
		if (drv){
			Mat toCopy = - bonesDrv;
			toCopy.copyTo(*drv);
		}
		return -bones;
	}

	float tension = 0.f;
	Mat tensionDrv(x.rows, x.cols, CV_32F, Scalar::all(0));
	Point3D * vrt = alg->sphere.vertices;
	for (int j = 0; j < alg->sphere.ne; j++){
		Edge e = alg->sphere.edges[j];
		
		Point3D v1 = vrt[e.verticeIndex1];
		Point3D v2 = vrt[e.verticeIndex2];
		float x1 = x.at<float>(0, e.verticeIndex1);
		float x2 = x.at<float>(0, e.verticeIndex2);
		Point3D p1 = Point3D_(v1.x * x1, v1.y * x1, v1.z * x1);
		Point3D p2 = Point3D_(v2.x * x2, v2.y * x2, v2.z * x2);

		float len2 = sqr(p1.x - p2.x) + sqr(p1.y - p2.y) + sqr(p1.z - p2.z);
		tension += len2;

		if (drv){
			tensionDrv.at<float>(0, e.verticeIndex1) += - 2 * (vrt[e.verticeIndex1].x * (p2.x - p1.x) + 
				vrt[e.verticeIndex1].y * (p2.y - p1.y) + vrt[e.verticeIndex1].z * (p2.z - p1.z));
			tensionDrv.at<float>(0, e.verticeIndex2) += - 2 * (vrt[e.verticeIndex2].x * (p1.x - p2.x) + 
				vrt[e.verticeIndex2].y * (p1.y - p2.y) + vrt[e.verticeIndex2].z * (p1.z - p2.z));
		}
	}
	tension = tension / alg->sphere.ne;
	tensionDrv = tensionDrv / alg->sphere.ne;

	float cost = - bones + alpha_tension * tension;

	if (drv){
		/*FileStorage fs("1.xml", FileStorage::WRITE);
		fs << "bonesDrv" << bonesDrv << "tensionDrv" << tensionDrv;*/

		Mat resDrv = - bonesDrv + tensionDrv * alpha_tension;
		resDrv.copyTo(*drv);
	}

	return cost;
}

Mat InitialIteration(){
	AlgorithmData * alg = getAlgorithmData();
	Mat x0(1, alg->sphere.nv, CV_32F, Scalar::all(0.5)), x0best;
	x0.copyTo(x0best);
	
	int i;
	/* without tension */
	float xx[] = {0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
	for (i = 0; i < x0.cols; i++){
		x0best.copyTo(x0);
		int ixx = 0;
		float bestCost = 1e6;
		for (int j = ARRAY_LEN(xx) - 1; j >= 0; j--){
			x0.at<float>(0, i) = xx[j];
			float cf = CostFunction(x0, NULL, 0.f);
			if (cf < bestCost){
				ixx = j;
				bestCost = cf;
			}
		}
		x0best.at<float>(0, i) = xx[ixx];
	}
	x0best.copyTo(x0);
	
	/* with tension */
	Vector<float> xx2;
	float xxv = 0.25;
	while (xxv < 1.25)
		xx2.push_back(xxv += 0.05);

	for (i = 0; i < x0.cols; i++){
		x0best.copyTo(x0);
		int ixx = 0;
		float bestCost = 1e6;
		for (int j = xx2.size() - 1; j >= 0; j--){
			x0.at<float>(0, i) = xx2[j];
			float cf = CostFunction(x0, NULL);
			if (cf < bestCost){
				ixx = j;
				bestCost = cf;
			}
		}
		x0best.at<float>(0, i) = xx2[ixx];
	}
	x0best.copyTo(x0);

	return x0;
}

Mat GradientDescend(Mat x0, int nIterations){
	Mat x, drv(x0.rows, x0.cols, CV_32F, Scalar::all(0));
	x0.copyTo(x);
	float step = 0.05f;

	char buf[1000];
	for (int i = 0; i < nIterations; i++){
		float c = CostFunction(x, &drv);
		if (i < 5 || i % (nIterations / 10) == 0){
			sprintf(buf, "Iter = %i,    Cost = %.5f", i, c);
			LOG(buf, false);
		}
		x = x - drv * step;
	}

	return x;
}

Surface MakeSegmentingShape(Mat x, short * dims){
	Surface surface;
	AlgorithmData * alg = getAlgorithmData();

	Point3D bias = Point3D_(-3, -1, 0);
	surface.ne = 0;
	surface.nv = alg->sphere.nv + 1;
	for (int iv = 0; iv < alg->sphere.nv; iv++){
		Point3i crd = VerticeToCoordinates(alg->sphere.vertices[iv], x.at<float>(iv), dims);
		surface.vertices[iv] = Point3D_(crd.x + bias.x, crd.y + bias.y, crd.z + bias.z);
	}
	Point3i p0 = VerticeToCoordinates(Point3D_(0, 0, 0), 0, dims);
	surface.vertices[surface.nv - 1] = Point3D_(p0.x + bias.x, p0.y + bias.y, p0.z + bias.z);

	surface.nf = alg->sphere.nf;
	for (int ifc = 0; ifc < surface.nf; ifc++){
		surface.faces[ifc] = alg->sphere.faces[ifc];
	}

	return surface;
}

void ShapeToSegmentedAnalyze(AnalyzeImage & im, Point3D* vertices, Face3* faces, int verticesCount, int facesCount){

	printf("(%i) : ", int(im.slices.size() / 2));
	Point3D p0 = vertices[verticesCount - 1]; // reference point
	for (int k = 0; k < im.slices.size() - 1; k = k + 2){
		printf((k + 2) % 20 == 0 ? "X" : "|");
		for (int j = 0; j < im.slices[0].cols - 3; j = j + 4){
			for (int i = 0; i < im.slices[0].cols - 3; i = i + 4){
				Point3D p = Point3D_(i, j, k);

				bool within = false;
				bool finish = false;
				int ifc = 0;
				while (!finish){
					Point3D p1 = vertices[faces[ifc].verticeIndex1];
					Point3D p2 = vertices[faces[ifc].verticeIndex2];
					Point3D p3 = vertices[faces[ifc].verticeIndex3];
					Point3D pmin = Point3D_(min(min(p0.x, p1.x), min(p2.x, p3.x)), min(min(p0.y, p1.y), min(p2.y, p3.y)), min(min(p0.z, p1.z), min(p2.z, p3.z)));
					Point3D pmax = Point3D_(max(max(p0.x, p1.x), max(p2.x, p3.x)), max(max(p0.y, p1.y), max(p2.y, p3.y)), max(max(p0.z, p1.z), max(p2.z, p3.z)));
					if (p.x >= pmin.x && p.x <= pmax.x && p.y >= pmin.y && p.y <= pmax.y && p.z >= pmin.z && p.z <= pmax.z){
						within = within || PointWithin(p, p0, p1, p2, p3);
					}

					ifc++;
					finish = within || ifc == facesCount;
				} 
				
				if (!within){
					for (int dk = 0; dk < 2; dk++)
					for (int dj = 0; dj < 4; dj++)
					for (int di = 0; di < 4; di++){
						im.slices[k + dk].at<short>(i + di, j + dj) = -2000;
					}
				}
			}
		}
	}

	if (im.slices.size() % 2 == 1){
		im.slices[im.slices.size() - 1] = Scalar::all(-2000);
	}
	/* Mark buoys
	for (int iv = 0; iv < verticesCount; iv++){
		Point3D p = vertices[iv];
		if (p.x > 0 && p.y > 0 && p.z > 0 && p.x < im.slices[0].rows - 4 && p.y < im.slices[0].cols - 4 && p.z < im.slices.size() - 2){
			for (int dk = 0; dk < 2; dk++)
			for (int dj = 0; dj < 4; dj++)
			for (int di = 0; di < 4; di++){
				im.slices[int(p.z) + dk].at<short>(int(p.x) + di, int(p.y) + dj) = 2000;
			}
		}
	}*/
	printf("\n");
}
