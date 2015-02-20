#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <io.h>
//#include <Windows.h>

#include "opencv2/opencv.hpp"

#include "common.h"
#include "semicommon.h"
#include "srfio.h"
#include "analyzecv/analyzecv.h"
#include "chsegm.h"
#include "fft.h"

#include <unistd.h>

#define ARRAY_LEN(x) (sizeof(x)/sizeof(x[0]))
