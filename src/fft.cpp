#include "stdafx.h"

void ditfft2(float x[4][4], std::complex<float> y[4][4], int startx[], int starty[], int Ns[], int step[]){
	printf("Ns = {%i, %i}\tstartx = {%i, %i}\tstarty = {%i, %i}\tstep = {%i, %i}\n", 
		Ns[0], Ns[1], startx[0], startx[1], starty[0], starty[1], step[0], step[1]);
	if (Ns[0] == 1 && Ns[1] == 1){
		y[starty[0]][starty[1]] = x[startx[0]][startx[1]];
	}
	else if (Ns[0] > Ns[1]){	// Ns[0] is devided by 2
		int startx1[2] = {startx[0], startx[1]};
		int starty1[2] = {starty[0], starty[1]};
		int step1[2] = {step[0] * 2, step[1]};
		int Ns1[2] = {Ns[0] / 2, Ns[1]};
		ditfft2(x, y, startx1, starty1, Ns1, step1);

		int startx2[2] = {startx[0] + step[0], startx[1]};
		int starty2[2] = {starty[0] + Ns[0] / 2, starty[1]};
		int step2[2] = {step[0] * 2, step[1]};
		int Ns2[2] = {Ns[0] / 2, Ns[1]};
		ditfft2(x, y, startx2, starty2, Ns2, step2);

		for (int i = 0; i < Ns[1]; i++){
			for (int k = 0; k < Ns[0] / 2; k++){
				std::complex<float> t = y[starty[0] + k][starty[1] + i];
				std::complex<float> e(0., 1);
				e = std::exp(e * std::complex<float>((-2) * PI() * k / float(Ns[0]), 0));
				y[starty[0] + k][starty[1] + i] = t + e * y[starty[0] + k + Ns[0] / 2][starty[1] + i];
				y[starty[0] + k + Ns[0] / 2][starty[1] + i] = t - e * y[starty[0] + k + Ns[0] / 2][starty[1] + i];
			}
		}
	}
	else{	// Ns[1] is devided by 2
		int startx1[2] = {startx[0], startx[1]};
		int starty1[2] = {starty[0], starty[1]};
		int step1[2] = {step[0], step[1] * 2};
		int Ns1[2] = {Ns[0], Ns[1] / 2};
		ditfft2(x, y, startx1, starty1, Ns1, step1);

		int startx2[2] = {startx[0], startx[1] + step[1]};
		int starty2[2] = {starty[0], starty[1] + Ns[1] / 2};
		int step2[2] = {step[0], step[1] * 2};
		int Ns2[2] = {Ns[0], Ns[1] / 2};
		ditfft2(x, y, startx2, starty2, Ns2, step2);

		for (int i = 0; i < Ns[0]; i++){
			for (int k = 0; k < Ns[1] / 2; k++){
				std::complex<float> t = y[starty[0] + i][starty[1] + k];
				std::complex<float> e(0., 1);
				e = std::exp(e * std::complex<float>((-2) * PI() * k / float(Ns[1]), 0));
				y[starty[0] + i][starty[1] + k] = t + e * y[starty[0] + i][starty[1] + k + Ns[1] / 2];
				y[starty[0] + i][starty[1] + k + Ns[1] / 2] = t - e * y[starty[0] + i][starty[1] + k + Ns[1] / 2];
			}
		}
	}
}

void FourierTest(){
	float x[4][4] = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 8, 7, 6}, {5, 4, 3, 2}};
	std::complex<float> y[4][4];

	int zrs[2] = {0, 0};
	int ons[2] = {1, 1};
	int frs[2] = {4, 4};
	ditfft2(x, y, zrs, zrs, frs, ons);

	Mat m(4, 4, CV_32FC2, y);
	FileStorage fs("fft.xml", FileStorage::WRITE);
	fs << "m" << m;
}

void ditfft3(AnalyzeImage &x, Mat &y, int startx[3], int starty[3], int N[3], int step[3], float sign){
	if (N[0] == 64 && N[1] == 64 && N[2] == 64)
		printf("|");
	if (N[0] == 1 && N[1] == 1 && N[2] == 1){
		std::complex<float> c(0, 0);
		if (startx[0] < x.slices[0].rows && startx[1] < x.slices[0].cols && startx[2] < x.slices.size())
			c = std::complex<float>(x.slices[startx[2]].at<float>(startx[0], startx[1]), 0);
        y.at<std::complex<float> >(starty[0], starty[1], starty[2]) = c;
	}
	else {
		if (N[0] >= N[1] && N[0] >= N[2]){	// Ns[0] is devided by 2
			int startx1[3] = {startx[0], startx[1], startx[2]};
			int starty1[3] = {starty[0], starty[1], starty[2]};
			int step1[3] = {step[0] * 2, step[1], step[2]};
			int N1[3] = {N[0] / 2, N[1], N[2]};
			ditfft3(x, y, startx1, starty1, N1, step1, sign);

			int startx2[3] = {startx[0] + step[0], startx[1], startx[2]};
			int starty2[3] = {starty[0] + N[0] / 2, starty[1], starty[2]};
			int step2[3] = {step[0] * 2, step[1], step[2]};
			int N2[3] = {N[0] / 2, N[1], N[2]};
			ditfft3(x, y, startx2, starty2, N2, step2, sign);

			for (int i = 0; i < N[0] / 2; i++)
				for (int j = 0; j < N[1]; j++)
					for (int k = 0; k < N[2]; k++){
                        std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
						std::complex<float> e(0., 1);
						e = std::exp(e * std::complex<float>(sign * (-2) * PI() * i / float(N[0]), 0));
                        y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                            t + e * y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k);
                        y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k) =
                            t - e * y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k);
					}
		}
		else if (N[1] >= N[0] && N[1] >= N[2]){	// Ns[1] is devided by 2
			int startx1[3] = {startx[0], startx[1], startx[2]};
			int starty1[3] = {starty[0], starty[1], starty[2]};
			int step1[3] = {step[0], step[1] * 2, step[2]};
			int N1[3] = {N[0], N[1] / 2, N[2]};
			ditfft3(x, y, startx1, starty1, N1, step1, sign);

			int startx2[3] = {startx[0], startx[1] + step[1], startx[2]};
			int starty2[3] = {starty[0], starty[1] + N[1] / 2, starty[2]};
			int step2[3] = {step[0], step[1] * 2, step[2]};
			int N2[3] = {N[0], N[1] / 2, N[2]};
			ditfft3(x, y, startx2, starty2, N2, step2, sign);

			for (int i = 0; i < N[0]; i++)
				for (int j = 0; j < N[1] / 2; j++)
					for (int k = 0; k < N[2]; k++){
                        std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
						std::complex<float> e(0., 1);
						e = std::exp(e * std::complex<float>(sign * (-2) * PI() * j / float(N[1]), 0));
                        y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                            t + e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k);
                        y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k) =
                            t - e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k);
					}
		}
		else {	// Ns[2] is devided by 2
			int startx1[3] = {startx[0], startx[1], startx[2]};
			int starty1[3] = {starty[0], starty[1], starty[2]};
			int step1[3] = {step[0], step[1], step[2] * 2};
			int N1[3] = {N[0], N[1], N[2] / 2};
			ditfft3(x, y, startx1, starty1, N1, step1, sign);

			int startx2[3] = {startx[0], startx[1], startx[2] + step[2]};
			int starty2[3] = {starty[0], starty[1], starty[2] + N[2] / 2};
			int step2[3] = {step[0], step[1], step[2] * 2};
			int N2[3] = {N[0], N[1], N[2] / 2};
			ditfft3(x, y, startx2, starty2, N2, step2, sign);

			for (int i = 0; i < N[0]; i++)
				for (int j = 0; j < N[1]; j++)
					for (int k = 0; k < N[2] / 2; k++){
                        std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
						std::complex<float> e(0., 1);
						e = std::exp(e * std::complex<float>(sign * (-2) * PI() * k / float(N[2]), 0));
                        y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                            t + e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2);
                        y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2) =
                            t - e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2);
					}
		}
	}
}

void ditfft3(Mat &x, Mat &y, int startx[3], int starty[3], int N[3], int step[3], float sign){
	if (N[0] == 64 && N[1] == 64 && N[2] == 64)
		printf("|");
	if (N[0] == 1 && N[1] == 1 && N[2] == 1){
		std::complex<float> c(0, 0);
		if (startx[0] < x.size[0] && startx[1] < x.size[1] && startx[2] < x.size[2])
            c = x.at<std::complex<float> >(startx[0], startx[1], startx[2]);
        y.at<std::complex<float> >(starty[0], starty[1], starty[2]) = c;
	}
	else if (N[0] >= N[1] && N[0] >= N[2]){	// Ns[0] is devided by 2
		int startx1[3] = {startx[0], startx[1], startx[2]};
		int starty1[3] = {starty[0], starty[1], starty[2]};
		int step1[3] = {step[0] * 2, step[1], step[2]};
		int N1[3] = {N[0] / 2, N[1], N[2]};
		ditfft3(x, y, startx1, starty1, N1, step1, sign);

		int startx2[3] = {startx[0] + step[0], startx[1], startx[2]};
		int starty2[3] = {starty[0] + N[0] / 2, starty[1], starty[2]};
		int step2[3] = {step[0] * 2, step[1], step[2]};
		int N2[3] = {N[0] / 2, N[1], N[2]};
		ditfft3(x, y, startx2, starty2, N2, step2, sign);

		for (int i = 0; i < N[0] / 2; i++)
			for (int j = 0; j < N[1]; j++)
				for (int k = 0; k < N[2]; k++){
                    std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
					std::complex<float> e(0., 1);
					e = std::exp(e * std::complex<float>(sign * (-2) * PI() * i / float(N[0]), 0));
                    y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                        t + e * y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k);
                    y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k) =
                        t - e * y.at<std::complex<float> >(starty[0] + i + N[0] / 2, starty[1] + j, starty[2] + k);
				}
	}
	else if (N[1] >= N[0] && N[1] >= N[2]){	// Ns[1] is devided by 2
		int startx1[3] = {startx[0], startx[1], startx[2]};
		int starty1[3] = {starty[0], starty[1], starty[2]};
		int step1[3] = {step[0], step[1] * 2, step[2]};
		int N1[3] = {N[0], N[1] / 2, N[2]};
		ditfft3(x, y, startx1, starty1, N1, step1, sign);

		int startx2[3] = {startx[0], startx[1] + step[1], startx[2]};
		int starty2[3] = {starty[0], starty[1] + N[1] / 2, starty[2]};
		int step2[3] = {step[0], step[1] * 2, step[2]};
		int N2[3] = {N[0], N[1] / 2, N[2]};
		ditfft3(x, y, startx2, starty2, N2, step2, sign);

		for (int i = 0; i < N[0]; i++)
			for (int j = 0; j < N[1] / 2; j++)
				for (int k = 0; k < N[2]; k++){
                    std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
					std::complex<float> e(0., 1);
					e = std::exp(e * std::complex<float>(sign * (-2) * PI() * j / float(N[1]), 0));
                    y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                        t + e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k);
                    y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k) =
                        t - e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j + N[1] / 2, starty[2] + k);
				}
	}
	else {	// Ns[2] is devided by 2
		int startx1[3] = {startx[0], startx[1], startx[2]};
		int starty1[3] = {starty[0], starty[1], starty[2]};
		int step1[3] = {step[0], step[1], step[2] * 2};
		int N1[3] = {N[0], N[1], N[2] / 2};
		ditfft3(x, y, startx1, starty1, N1, step1, sign);

		int startx2[3] = {startx[0], startx[1], startx[2] + step[2]};
		int starty2[3] = {starty[0], starty[1], starty[2] + N[2] / 2};
		int step2[3] = {step[0], step[1], step[2] * 2};
		int N2[3] = {N[0], N[1], N[2] / 2};
		ditfft3(x, y, startx2, starty2, N2, step2, sign);

		for (int i = 0; i < N[0]; i++)
			for (int j = 0; j < N[1]; j++)
				for (int k = 0; k < N[2] / 2; k++){
                    std::complex<float> t = y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k);
					std::complex<float> e(0., 1);
					e = std::exp(e * std::complex<float>(sign * (-2) * PI() * k / float(N[2]), 0));
                    y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k) =
                        t + e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2);
                    y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2) =
                        t - e * y.at<std::complex<float> >(starty[0] + i, starty[1] + j, starty[2] + k + N[2] / 2);
				}
	}
}

int nextpow2(int num){
	int n = 1;
	while (n < num)
		n *= 2;
	return n;
}

Mat AnalyzeFFT(AnalyzeImage &x, bool inverse){
	char buf[1000];
	if (~inverse)
        sprintf(buf, "FFT on %dx%dx%lu Analyze image", x.slices[0].rows, x.slices[0].cols, x.slices.size());
	else
        sprintf(buf, "iFFT on %dx%dx%lu Analyze image", x.slices[0].rows, x.slices[0].cols, x.slices.size());
	LOG(buf, false, true);

	int dims[3] = {nextpow2(x.slices[0].rows), nextpow2(x.slices[0].cols), nextpow2(x.slices.size())};
	Mat y(3, dims, CV_32FC2, Scalar::all(0));

	int zrs[3] = {0, 0, 0};
	int ons[3] = {1, 1, 1};
	printf("(%i) : ", (int)sqr(dims[0] / 64));
	ditfft3(x, y, zrs, zrs, dims, ons, inverse ? -1. : 1);
	printf("\n");

	if (inverse)
		y = y / (dims[0] * dims[1] * dims[2]);

	return y;
}

Mat Mat3FFT(Mat &x, bool inverse, int * N3){
	char buf[1000];
	if (~inverse)
		sprintf(buf, "FFT on %ix%ix%i complex Mat", x.size[0], x.size[1], x.size[2]);
	else
		sprintf(buf, "iFFT on %ix%ix%i complex Mat", x.size[0], x.size[1], x.size[2]);
	LOG(buf, false, true);

	int dims[3];
	if (N3 == 0){
		dims[0] = nextpow2(x.size[0]);
		dims[1] = nextpow2(x.size[1]);
		dims[2] = nextpow2(x.size[2]);
	}
	else {
		dims[0] = N3[0];
		dims[1] = N3[1];
		dims[2] = N3[2];
	}

	Mat y(3, dims, CV_32FC2, Scalar::all(0));

	int zrs[3] = {0, 0, 0};
	int ons[3] = {1, 1, 1};
	printf("(%i) : ", (int)sqr(dims[0] / 64));
	ditfft3(x, y, zrs, zrs, dims, ons, inverse ? -1. : 1);
	printf("\n");

	if (inverse)
		y = y / (dims[0] * dims[1] * dims[2]);

	return y;
}
