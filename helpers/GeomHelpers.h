#pragma once

#if defined(GEOMHELPERS_RECURSES)
#error Recursive header files inclusion detected in Geomhelpers.h
#else // defined(GEOMHELPERS_RECURSES)
/** Prevents recursive inclusion of headers. */
#define GEOMHELPERS_RECURSES

#if !defined GEOMHELPERS_h
/** Prevents repeated inclusion of headers. */
#define GEOMHELPERS_h

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include <numeric>
#include <algorithm>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <exception>
#include <stdio.h>
#include <cstdlib>

using namespace std;
using namespace DGtal;

typedef Viewer3D<> MyViewer;
typedef AlphaThickSegmentComputer< Z2i::Point > AlphaThickSegmentComputer2D;
typedef AlphaThickSegmentComputer<Z2i::RealPoint> AlphaThickSegmentComputer2DD;

//https://stackoverflow.com/questions/1119627/how-to-test-if-a-point-is-inside-of-a-convex-polygon-in-2d-integer-coordinates
bool isInConvexPolygon(Z2i::Point testPoint, std::vector<Z2i::Point> polygon);

void writePointList(const vector<Z3i::Point>& v, const char* filename);
void writeRealPointList(const vector<Z3i::RealPoint>& v, const char* filename);
void writeVectorList(const vector<double> &v, const char *filename);

vector<Z3i::Point> readPointList(const char* filename);
vector<Z3i::RealPoint> readRealPointList(const char* filename);
vector<double> readVectorList(const char *filename);

Color getColor(int i);

bool projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz, vector<std::tuple<bool, bool, bool> >& projection);
bool projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz, bool& validOxy, bool& validOxz, bool& validOyz);
void projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz);

vector<double> getMeaningfulThickness(const string ImaGeneDir, const vector<Z2i::Point> aContour, double maxMT=10, double stepMT=1);
vector<double> sortMeaningfulThickness(const vector<double> vecMT);
vector<double> labelMeaningfulThickness(const vector<Z2i::Point>& aContour, const vector<double>& vecMT);

void getParameters (const AlphaThickSegmentComputer2D& mySeg_Pl1,const AlphaThickSegmentComputer2D& mySeg_Pl2, const bool validOxy, const bool validOxz, const bool validOyz, Z3i::RealPoint& direction, Z3i::RealPoint& intercept, Z3i::RealPoint& thickness);

void drawAsBoundingBox(MyViewer& display, const vector<Z3i::Point>& aContour, size_t idBegin, size_t idEnd, const AlphaThickSegmentComputer2D& mySeg_Pl1, const AlphaThickSegmentComputer2D& mySeg_Pl2, const bool validOxy, const bool validOxz, const bool validOyz);

void getTangentCoverDifferentNoise(const vector<double>& vecThickness, const vector<Z2i::Point>& aContour_Pl1, const vector<Z2i::Point>& aContour_Pl2,
        vector<vector<AlphaThickSegmentComputer2D> >& vecTC_Pl1, vector<vector<AlphaThickSegmentComputer2D>>&vecTC_Pl2);

void getSegmentsAdaptiveTangentCover(const vector<Z3i::Point>& aContour, const vector<double>& vecThickness, const vector<double>& vecMT3D, const vector<vector<AlphaThickSegmentComputer2D> >& vecTC_Pl1, const vector<vector<AlphaThickSegmentComputer2D>>&vecTC_Pl2,
        vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl1, vector<vector<AlphaThickSegmentComputer2D>>& vecATC_Pl2, const bool validOxy, const bool validOxz, const bool validOyz);

void buildAdaptiveTangentCover(const vector<Z3i::Point>& aContour, const vector<double>& vecThickness, const vector<double>& vecMT3D, const vector<Z2i::Point>& aContour_Pl1, const vector<Z2i::Point>& aContour_Pl2,
        vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl1, vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl2, const bool validOxy, const bool validOxz, const bool validOyz);

void drawAdaptiveTangentCover(MyViewer& display, const vector<Z3i::Point>& aContour, const vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl1, const vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl2, const bool validOxy, const bool validOxz, const bool validOyz);

void drawBoundingBox(MyViewer& display, Z3i::RealPoint pf, Z3i::RealPoint pl, Z3i::RealPoint pf2, Z3i::RealPoint pl2, Z3i::RealPoint pf3, Z3i::RealPoint pl3, Z3i::RealPoint pf4, Z3i::RealPoint pl4);

#endif // !defined GEOMHELPERS_h

#undef GEOMHELPERS_RECURSES
#endif // else defined(GEOMHELPERS_RECURSES)

