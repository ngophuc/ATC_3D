////////////////////////////////////////////////////////////////////////////
// Display the Noise level obtained after multi scales analysis
//////////////////////////////////////////////////////////////////////////////
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

static Arguments args;
static const int RESOLUTION=1200;





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

vector<double>
meaningfulThicknessEstimator(const vector<Vector2D>& polygon, bool& isOpenPolygon, double samplingSizeMax, double samplingStep)
{
  BlurredSegmentTgtCover tgc;
  bool loop = ContourHelper::containsLoop(polygon);
  if(loop) {
    cerr << "contains loop aborting..." << endl;
    //exit(0);
  }else{
    cerr << "no loop ok" << endl;
  }
  isOpenPolygon = ContourHelper::isOpenPolygon(polygon,3.0);
  cerr << "contour considered as :" << (isOpenPolygon ? "Open" :"Closed") << endl <<
  "polygon size " << polygon.size() << endl;
  tgc.init(polygon, !isOpenPolygon);
  
  vector<double> vectScales;
  for(double i=1.0; i<samplingSizeMax; i=i+samplingStep){
    vectScales.push_back((double)i );
  }
  
  vector<double> noiseLevels = tgc.getNoiseLevels(vectScales, 1, -0.0);
  vector<Vector2D> contour = tgc.getPointsContour();
  
  vector<double> vecMT;
  for(int i=0; i< contour.size(); i++){
    Vector2D ptA = contour.at(i);
    //if(isOpenPolygon && (i+1==contour.size()) ) break;
    int noiseAindex = (noiseLevels.at(i)==0)? (vectScales.size()-1): noiseLevels.at(i) ;
    double noiseLevelA= (vectScales.at(noiseAindex-1)/2.0);
    vecMT.push_back(noiseLevelA);
  }
  if(isOpenPolygon)
    vecMT.front() = vecMT.back();
     
  return vecMT;
}

int
main( int argc, char** argv )
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addIOArgs( args, true, false );
  
  args.addOption("-srcPolygon", "-srcPolygon <contour.sdp> <posX> <posY>", " ", "0", "1" );
  args.addOption("-setSampling", "-setSampling <maxScale> <samplingStep>: set the maximal scale and sampling step used for contour analysis (default: maxScale=15.0, samplingStep=1.0)",
                 "15.0", "1.0");
  args.addOption("-exportNoiseLevel", "-exportNoiseLevel <filename> export noise level", "noise.dat" );
  if ( ( argc <= 0 )  || ! args.readArguments( argc, argv ) || !args.check("-srcPolygon")  ) {
    cerr << args.usage( "displayNoiseBS", "displayNoiseBS -srcPolygon.","" ) << endl;
    return 1;
  }
  
  vector<Vector2D> polygon;
  if(args.check("-srcPolygon")){
    string fileName = args.getOption("-srcPolygon")->getValue(0);
    uint posX = args.getOption("-srcPolygon")->getIntValue(1);
    uint posY = args.getOption("-srcPolygon")->getIntValue(2);
    string closedOption =args.getOption("-srcPolygon")->getValue(3);
    ifstream ifs ( fileName.c_str() , ifstream::in );
    polygon =  ContourHelper::getPolygonFromStream(ifs, posX, posY);
  }
  
  bool exportNoise = args.check("-exportNoiseLevel");
  string nameExport = args.getOption("-exportNoiseLevel")->getValue(0);
  ofstream outNoiseLevel;
  if(exportNoise){
    outNoiseLevel.open( nameExport.c_str() , ofstream::out );
  }
  
  double samplingSizeMax = args.getOption("-setSampling")->getFloatValue(0);
  double samplingStep = args.getOption("-setSampling")->getFloatValue(1);
  
  cerr << "Computing Meaningful Thickness" << endl;
  
  vector<double> vectScales;
  
  for(double i=1.0; i<samplingSizeMax; i=i+samplingStep){
    vectScales.push_back((double)i );
  }
  
  bool isOpenPolygon;
  vector<double> noiseLevels= meaningfulThicknessEstimator(polygon, isOpenPolygon, samplingSizeMax, samplingStep);
  assert(polygon.size()==noiseLevels.size());
    
  if(exportNoise){
    outNoiseLevel << "# Noise levels exported from displayNoiseBS (source ImaGene): " << std::endl;
    outNoiseLevel << "# Format X Y noiseLevel " << std::endl;
    for(int i=0; i< polygon.size(); i++){
      Vector2D ptA = polygon.at(i);
      //if(isOpenPolygon && (i+1==polygon.size())) break;
      outNoiseLevel << ptA.x() << " " << ptA.y() << " "<<  noiseLevels.at(i) << std::endl;
    }
    outNoiseLevel.close();
  }

  return 0;
}
