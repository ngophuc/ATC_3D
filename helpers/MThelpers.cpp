////////////////////////////////////////////////////////////////////////////
// Display the Noise level obtained after multi scales analysis
//////////////////////////////////////////////////////////////////////////////
#include "MThelpers.h"

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

void
meaningfulThicknessEstimator(string infile, double samplingSizeMax, double samplingStep, string outfile)
{
  vector<Vector2D> polygon;
  ifstream ifs ( infile.c_str() , ifstream::in );
  uint posX = 0;
  uint posY = 1;
  polygon =  ContourHelper::getPolygonFromStream(ifs, posX, posY);

  ofstream outNoiseLevel;
  outNoiseLevel.open( outfile.c_str() , ofstream::out );
  
  cerr << "Computing Meaningful Thickness" << endl;
  
  vector<double> vectScales;
  for(double i=1.0; i<samplingSizeMax; i=i+samplingStep)
    vectScales.push_back((double)i );
  
  bool isOpenPolygon;
  vector<double> noiseLevels= meaningfulThicknessEstimator(polygon, isOpenPolygon, samplingSizeMax, samplingStep);
  assert(polygon.size()==noiseLevels.size());
    
  outNoiseLevel << "# Noise levels exported from displayNoiseBS (source ImaGene): " << std::endl;
  outNoiseLevel << "# Format X Y noiseLevel " << std::endl;
  for(int i=0; i< polygon.size(); i++){
    Vector2D ptA = polygon.at(i);
    //if(isOpenPolygon && (i+1==polygon.size())) break;
    outNoiseLevel << ptA.x() << " " << ptA.y() << " "<<  noiseLevels.at(i) << std::endl;
  }
  outNoiseLevel.close();
}
