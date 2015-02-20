
struct Surface {
	struct Point3D vertices[3000];
	struct Face3 faces[3000];
	struct Edge edges[3000];
	int nv, nf, ne;
};

struct AlgorithmData {
	struct AnalyzeImage im;			// initial image
	struct AnalyzeImage blurred;	// binarized and blurred image
	struct AnalyzeImage gr[3];		// three components 
	struct Surface sphere;			// the initial spherical surface for segmentation
	cv::Point3f center;				// reference point
	float z2xy;						// pixel sizes ratio = dim_Z / dim_X
};

AlgorithmData * getAlgorithmData();

void ShapeToSegmentedAnalyze(struct AnalyzeImage & im, struct Point3D* vertices, struct Face3* facets, int verticesCount, int facesCount);

void PreprocessImage();

float CostFunction(cv::Mat x, cv::Mat * drv = NULL, float alpha_tension = -1);

cv::Mat InitialIteration();

struct Surface MakeSegmentingShape(cv::Mat x, short * dims);

Mat GradientDescend(Mat x0, int nIterations = 2000);