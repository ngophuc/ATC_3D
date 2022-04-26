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
meaningfulThicknessEstimator(vector<Vector2D> polygon, double samplingSizeMax, double samplingStep)
{
  BlurredSegmentTgtCover tgc;
  bool loop = ContourHelper::containsLoop(polygon);
  if(loop) {
    cerr << "contains loop aborting..." << endl;
    exit(0);
  }else{
    cerr << "no loop ok" << endl;
  }
  bool isOpenPolygon = ContourHelper::isOpenPolygon(polygon,3.0);
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
     if(isOpenPolygon && (i+1==contour.size()) )
       break;
     int noiseAindex = (noiseLevels.at(i)==0)? (vectScales.size()-1): noiseLevels.at(i) ;
     double noiseLevelA= (vectScales.at(noiseAindex-1)/2.0);
    vecMT.push_back(noiseLevelA);
  }
  
  return vecMT;
}

int
main( int argc, char** argv )
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  StandardArguments::addIOArgs( args, true, false );

  args.addOption("-srcPolygon", "-srcPolygon <contour.sdp> <posX> <posY> <CLOSED|OPEN> (default CLOSED)", " ", "0", "1", "CLOSED" );
  args.addBooleanOption("-estimClosedContour", "-estimClosedContour estimate is the contour is closed or not" );
  args.addOption("-setSampling", "-setSampling <maxScale> <samplingStep>: set the maximal scale and sampling step used for contour analysis (default: maxScale=15.0, samplingStep=1.0)",
     "15.0", "1.0");
  args.addOption("-exportNoiseLevel", "-exportNoiseLevel <filename> export noise level", "noise.dat" );
    if ( ( argc <= 0 )  || ! args.readArguments( argc, argv ) || !args.check("-srcPolygon")  ) {
      cerr << args.usage( "displayNoiseBS",
        "displayNoiseBS  -srcFC or -srcPolygon."
        ,"" ) << endl;
      return 1;
    }
  
  BlurredSegmentTgtCover tgc;

  bool isOpenPolygon;
  if(args.check("-srcPolygon")){
    string fileName = args.getOption("-srcPolygon")->getValue(0);
    uint posX = args.getOption("-srcPolygon")->getIntValue(1);
    uint posY = args.getOption("-srcPolygon")->getIntValue(2);
    string closedOption =args.getOption("-srcPolygon")->getValue(3);
    ifstream ifs ( fileName.c_str() , ifstream::in );
    vector<Vector2D> polygon =  ContourHelper::getPolygonFromStream(ifs, posX, posY);
    bool loop = ContourHelper::containsLoop(polygon);
    if(loop) {
      cerr << "contains loop aborting..." << endl;
      return 0;
    }else{
      cerr << "no loop ok" << endl;
    }
    if(args.check("-estimClosedContour")){
      isOpenPolygon = ContourHelper::isOpenPolygon(polygon,3.0);
    }else{
      isOpenPolygon=!(closedOption=="CLOSED");
    }
    cerr << "contour considered as :" << (isOpenPolygon ? "Open" :"Closed") << endl <<
      "polygon size " << polygon.size() << endl;
    tgc.init(polygon, !isOpenPolygon);
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

  vector<double> noiseLevels= tgc.getNoiseLevels(vectScales, 1, -0.0);
  vector<Vector2D> contour = tgc.getPointsContour();
  
  if(exportNoise){
    outNoiseLevel << "# Noise levels exported from displayNoiseBS (source ImaGene): " << std::endl;
    outNoiseLevel << "# Format X Y noiseLevel " << std::endl;
  }
     
  
  for(int i=0; i< contour.size(); i++){
     Vector2D ptA = contour.at(i);
     if(isOpenPolygon && (i+1==contour.size()) )
       break;
     
     //Vector2D ptB = contour.at((i+1)%contour.size());

     int noiseAindex = (noiseLevels.at(i)==0)? (vectScales.size()-1): noiseLevels.at(i) ;
     //int noiseBindex = (noiseLevels.at((i+1)%contour.size())==0)?(vectScales.size()-1): noiseLevels.at((i+1)%contour.size())  ;
     
     double noiseLevelA= (vectScales.at(noiseAindex-1)/2.0);
     //double noiseLevelB= (vectScales.at(noiseBindex-1)/2.0);

     if(exportNoise){
       outNoiseLevel << ptA.x() << " " << ptA.y() << " "<<  noiseLevelA << std::endl;
     }
  }
  if(exportNoise){
    outNoiseLevel.close();
  }
  
}















