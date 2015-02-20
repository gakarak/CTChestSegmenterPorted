
double PI();

char * buf();

float sqr(float x);

float dist2(float x1, float y1, float x2, float y2);

float dist3(float x1, float y1, float z1, float x2, float y2, float z2);

int signum(int x);

void order2intAsc(int* p1, int* p2);

void ClearFile(const char * filename);

void CropFilenameExtention(char* s);

cv::Mat dlmread(const char * filename, char delim = '\t');

void dlmwrite(const char * filename, cv::Mat &m, char delim = '\t');

cv::Mat colread(const char * filename);