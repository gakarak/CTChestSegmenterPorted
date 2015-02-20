
void FourierTest();

cv::Mat AnalyzeFFT(struct AnalyzeImage &x, bool inverse = false);

cv::Mat Mat3FFT(cv::Mat &x, bool inverse = false, int * N3 = 0);