#pragma once

#if defined(MTHELPERS_RECURSES)
#error Recursive header files inclusion detected in MThelpers.h
#else // defined(MTHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MTHELPERS_RECURSES

#if !defined MTHELPERS_h
/** Prevents repeated inclusion of headers. */
#define MTHELPERS_RECURSES

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <vector>
#include <algorithm>
#include "ImaGene/base/Proxy.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/dgeometry2d/C4CSegmentPencil.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/timetools/Clock.h"
#include "ImaGene/dgeometry2d/C4CTangentialCover.h"
#include "ImaGene/helper/C4CTangentialCoverGeometry.h"
#include "ImaGene/dgeometry2d/FreemanChain.h"
#include "ImaGene/dgeometry2d/FreemanChainTransform.h"
#include "ImaGene/dgeometry2d/C4CIteratorOnFreemanChain.h"
#include "ImaGene/digitalnD/C4CIteratorOnFreemanChainSurface.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/C4CIteratorOnSurface.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/helper/ScaleProfile.h"
#include "ImaGene/helper/MultiscaleProfile.h"
#include "ImaGene/helper/DrawingXFIG.h"
#include "ImaGene/helper/CurveVariationsHelper.h"
#include "ImaGene/helper/ContourHelper.h"

#include "ImaGene/mathutils/SimpleLinearRegression.h"
#include "ImaGene/dgeometry2d/BlurredSegmentTgtCover.h"

using namespace std;
using namespace ImaGene;

vector<double>
meaningfulThicknessEstimator(const vector<Vector2D>& polygon,
                             bool& isOpenPolygon,
                             double samplingSizeMax,
                             double samplingStep);

void
meaningfulThicknessEstimator(string infile,
                             double samplingSizeMax,
                             double samplingStep,
                             string outfile);

#endif // !defined MTHELPERS_h

#undef MTHELPERS_RECURSES
#endif // else defined(MTHELPERS_RECURSES)
