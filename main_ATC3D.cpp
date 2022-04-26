#include <DGtal/io/readers/PointListReader.h>
#include <random>

#ifdef WITH_VISU3D_QGLVIEWER
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#endif

#include "helpers/GeomHelpers.h"
#include "helpers/maximalblurredsegment3D.h"
#include "helpers/GeomHelpers.h"

#include "CLI11.hpp"

/* Build TC 3D with a given thickness */
vector<MaximalBlurredSegment3D> buildTC(std::vector<Z3i::Point>& aContour,
                                        const std::vector<Z2i::Point>& aContourXY,
                                        const std::vector<Z2i::Point>& aContourXZ,
                                        const std::vector<Z2i::Point>& aContourYZ,
                                        double thickness) {
  vector<MaximalBlurredSegment3D> TC;
  Iterator e = aContour.end();
  int it=0;
  int idS, idE;
  for (Iterator i = aContour.begin() ; i != e && idE+1!= aContour.size(); ++i) {
    AlphaThickSegmentComputer3D seg(thickness, 3); //2 or 3 valid planes
    seg.init(i); 
    idS = it;
    idE = it;
    bool isOk = true, okXY, okYZ, okXZ;
    while ( (seg.end() != e) && isOk && seg.isExtendableFront() ) {
      okXY = isInConvexPolygon(aContourXY.at(idE), seg.alphaThickSeg2dXY().getConvexHull());
      okYZ = isInConvexPolygon(aContourYZ.at(idE), seg.alphaThickSeg2dYZ().getConvexHull());
      okXZ = isInConvexPolygon(aContourXZ.at(idE), seg.alphaThickSeg2dXZ().getConvexHull());
      isOk = okXY && okYZ && okXZ;
      if(isOk) {
        seg.extendFront();
        idE++;
      }
    }
    //std::cout<<"End extendFront : ("<<idS<<"-"<<idE<<")"<<std::endl;
    if(it==0){
      MaximalBlurredSegment3D bs(idS, idE, thickness, seg);
      TC.push_back(bs);
    }
    else {
      if(TC.back().getIdEnd()<idE) {//remove seg by inclusion
        MaximalBlurredSegment3D bs(idS, idE, thickness, seg);
        TC.push_back(bs);
      }
    }
    it++;
  }
  return TC;
}
/* Build TC 3D with a given thickness */

/* Associate thickness to 3D points from its projections */
double max_functor(double x, double y, double z) {
  return std::max(x, std::max(y, z));
}
double min_functor(double x, double y, double z) {
  return std::min(x, std::min(y, z));
}
double mean_functor(double x, double y, double z) {
  return (x+y+z)/3.0;
}
double mid_functor(double x, double y, double z) {
  return std::max(std::min(x,y), std::min(std::max(x,y),z));
}
std::vector<double> getThicknessAssociation(std::vector<double> thicknessXY,
                                            std::vector<double> thicknessXZ,
                                            std::vector<double> thicknessYZ,
                                            std::function<double(double, double, double)> functor) {
  std::vector<double> vecThick;
  double thickness=0;
  for(size_t it=0; it<thicknessXY.size(); it++) {
    thickness = functor(thicknessXY.at(it), thicknessXZ.at(it), thicknessYZ.at(it));
    vecThick.push_back(thickness);
  }
  return vecThick;
}

/* Compute ATC 3D */
double maxThicknessSegment(const MaximalBlurredSegment3D& seg, const vector<double>& thckVect)
{
  int idStart=seg.getIdBegin();
  int idEnd=seg.getIdEnd();
  double maxThickness=-1;
  for(int i=idStart; i<=idEnd; i++) {
    double thicknessPoint=thckVect.at(i%thckVect.size());
    if(thicknessPoint>maxThickness)
      maxThickness=thicknessPoint;
  }
  return maxThickness;
}

std::vector<MaximalBlurredSegment3D> buildATC(std::vector<Z3i::Point>& aContour,
                                              const std::vector<double>& vecMT,
                                              bool isClosed = true) {
  //1. Find vector of thickness element
  vector<double> meaningThicknessElement;
  meaningThicknessElement.push_back(vecMT.front());
  for(vector<double>::const_iterator it=vecMT.begin()+1; it!=vecMT.end(); it++) {
    double m=(*it);
    if(std::find(meaningThicknessElement.begin(), meaningThicknessElement.end(),m)==meaningThicknessElement.end())
      meaningThicknessElement.push_back(m);
  }
  std::sort(meaningThicknessElement.begin(),meaningThicknessElement.end());
  for(vector<double>::const_iterator it=meaningThicknessElement.begin(); it!=meaningThicknessElement.end(); it++)
    cout<<"meaningThicknessElement : "<<*it<<endl;
  
  //2. Compute different thickness tangent covers(blurred segments)
  std::vector<std::vector<MaximalBlurredSegment3D> > meaningThicknessTangentCover(meaningThicknessElement.size());
  int index=0;
  for(vector<double>::const_iterator it=meaningThicknessElement.begin(); it!=meaningThicknessElement.end(); it++) {
    double thickness=(*it);//*sqrt(2);
    vector<Z2i::Point> aContourXY, aContourXZ, aContourYZ;
    bool validOxy=false,validOxz=false,validOyz=false;
    projectionPoints3D(aContour,aContourXY,aContourXZ,aContourYZ,validOxy,validOxz,validOyz);
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=buildTC(aContour,aContourXY,aContourXZ,aContourYZ,thickness);
    cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl;
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++)
      meaningThicknessTangentCover[index].push_back(*it_bis);
    index++;
  }
  
  //3. Update thickness of points with tangent covers
  vector<double> vecMTmodified;
  for(vector<double>::const_iterator it=vecMT.begin(); it!=vecMT.end(); it++)
    vecMTmodified.push_back(*it);
  for(int it=(int)meaningThicknessTangentCover.size()-1; it>=0; it--) {
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=meaningThicknessTangentCover.at(it);
    double thickness=meaningThicknessElement.at(it);
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++) {
      int idStart=(*it_bis).getIdBegin();
      int idEnd=(*it_bis).getIdEnd();
      if(idStart!=-1 && idEnd!=-1) {
        double maxThickness=maxThicknessSegment(*it_bis, vecMT);
        for(int i=idStart; i<=idEnd; i++) { //FIXME : idStart+1 => justify !!!
          if(maxThickness==thickness)
            vecMTmodified.at(i%aContour.size())=maxThickness;
        }
      }
      else
        cout<<"negatif"<<endl;
    }
  }
  
  //4. Travel over the tangent covers and select the segments w.r.t the associated thickness of points
  vector<vector<MaximalBlurredSegment3D> > adaptiveMeaningThicknessTangentCover;
  int idCover=0;
  for(vector<vector<MaximalBlurredSegment3D> >::const_iterator it=meaningThicknessTangentCover.begin(); it!=meaningThicknessTangentCover.end(); it++) {
    adaptiveMeaningThicknessTangentCover.push_back(vector<MaximalBlurredSegment3D>());
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=*it;
    vector<MaximalBlurredSegment3D> AdaptiveFuzzySegmentSet;
    int idSeg=0;
    double thickness=meaningThicknessElement.at(idCover);
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++) {
      int idStart=(*it_bis).getIdBegin();
      int idEnd=(*it_bis).getIdEnd();
      if(idStart!=-1 && idEnd!=-1) {
        bool isGoodMTmodif=false, isGoodMT=true;//true
        for(int i=idStart; i<=idEnd; i++) {
          double thicknessMT=vecMT.at(i%aContour.size()); //all elt have same meaningful thickness value(dont contain other meaningful thickness)
          if(thicknessMT==thickness)
            isGoodMT=true;
          double thicknessMTmodif=vecMTmodified.at(i%aContour.size());
          if(thicknessMTmodif==thickness) //there exist at least one elt in modif having meaningful thickness value
            isGoodMTmodif=true;
        }
        if(isGoodMTmodif && isGoodMT)
          AdaptiveFuzzySegmentSet.push_back(*it_bis);
      }
      idSeg++;
    }
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=AdaptiveFuzzySegmentSet.begin();it_bis!=AdaptiveFuzzySegmentSet.end();it_bis++)
      adaptiveMeaningThicknessTangentCover[idCover].push_back(*it_bis);
    idCover++;
  }
  
  // 4 cont. : selection of Seg
  for(vector<vector<MaximalBlurredSegment3D> >::reverse_iterator it1=adaptiveMeaningThicknessTangentCover.rbegin(); it1!=adaptiveMeaningThicknessTangentCover.rend(); ++it1) {
    vector<MaximalBlurredSegment3D>& segmentSet1=*it1;
    for(vector<vector<MaximalBlurredSegment3D> >::reverse_iterator it2=it1+1; it2!=adaptiveMeaningThicknessTangentCover.rend(); ++it2) {
      vector<MaximalBlurredSegment3D>& segmentSet2=*it2;
      for(vector<MaximalBlurredSegment3D>::iterator itSeg1=segmentSet1.begin();itSeg1!=segmentSet1.end();itSeg1++) {
        //get idStart and idEnd of current seg
        int idCurrentStart=(*itSeg1).getIdBegin();
        int idCurrentEnd=(*itSeg1).getIdEnd();
        
        for(vector<MaximalBlurredSegment3D>::iterator itSeg2=segmentSet2.begin();itSeg2!=segmentSet2.end();itSeg2++) {
          int idStart=(*itSeg2).getIdBegin();
          int idEnd=(*itSeg2).getIdEnd();
          if(idCurrentStart<=idStart && idCurrentEnd>=idEnd) {
            segmentSet2.erase(itSeg2);
            itSeg2--;
          }
        }
      }
    }
  }
  
  //5. Reorder the multi-thickness tangent cover
  vector<MaximalBlurredSegment3D> adaptiveTangentCover;
  int seg=0,nbSeg=0;
  vector<int> idThicknessCover; //stock idSeg of the last seg at idThicknessCover
  for(size_t it=0; it<meaningThicknessElement.size();it++)
    idThicknessCover.push_back(0);
  for(size_t it=0; it<adaptiveMeaningThicknessTangentCover.size(); it++)
    nbSeg+=(adaptiveMeaningThicknessTangentCover.at(it)).size();
  
  while(seg<nbSeg) {
    int idMinStart,idMinEnd;
    if(isClosed) {
      idMinStart=2*aContour.size();
      idMinEnd=2*aContour.size();//x2 : due to the closed curve can have ending index over contour.size()
    }
    else {
      idMinStart=aContour.size();
      idMinEnd=aContour.size();
    }
    int idMin=-1, idSeg=-1;
    //scan adaptiveMeaningThicknessTangentCover
    for(size_t it=0; it<adaptiveMeaningThicknessTangentCover.size(); it++) {//thickness level=it
      //current seg of thickness level idThicknessCover.at(i)
      int idCurrentSeg=idThicknessCover.at(it);
      if(idCurrentSeg<(adaptiveMeaningThicknessTangentCover.at(it)).size()) {
        //get idStart and idEnd of seg
        int idCurrentStart=(adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg).getIdBegin();
        int idCurrentEnd=(adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg).getIdEnd();
        if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size()) {
          //find min idCurrentStart
          if(idMinStart==idCurrentStart && idMinEnd<idCurrentEnd) {
            if(idThicknessCover.at(it)<(int)(adaptiveMeaningThicknessTangentCover.at(it)).size()-1) {
              idThicknessCover.at(idMin)=idThicknessCover.at(idMin)+1;
              seg++;
            }
            idSeg=idCurrentSeg;
            idMin=it;
            idMinStart=idCurrentStart;
            idMinEnd=idCurrentEnd;
          }
          else if(idMinStart>idCurrentStart && idMinEnd>=idCurrentEnd) {
            idSeg=idCurrentSeg;
            idMin=it;
            idMinStart=idCurrentStart;
            idMinEnd=idCurrentEnd;
          }
        }
      }
    }
    adaptiveTangentCover.push_back((adaptiveMeaningThicknessTangentCover.at(idMin)).at(idSeg));
    idThicknessCover.at(idMin)=idThicknessCover.at(idMin)+1;
    seg++;
  }
  return adaptiveTangentCover;
}
/* Compute ATC 3D */

/* Compute ATS 3D */
double statThicknessSegment(const MaximalBlurredSegment3D& seg, const vector<double>& thckVect)
{
  int idStart=seg.getIdBegin();
  int idEnd=seg.getIdEnd();
  vector<double> segThick;
  for(int i=idStart; i<=idEnd; i++) {
    double thicknessPoint=thckVect.at(i%thckVect.size());
    segThick.push_back(thicknessPoint);
  }
  
  double maxThick=-1;
  int countMax=-1, countThick=-1;
  for(vector<double>::const_iterator it=segThick.begin(); it!=segThick.end(); it++) {
    countThick=std::count(segThick.begin(), segThick.end(), *it);
    if(countThick>countMax) {
      maxThick=*it;
      countMax=countThick;
    }
  }
  return maxThick;
}

vector<MaximalBlurredSegment3D> buildModifiedATC(vector<Z3i::Point>& aContour,
                                                 const vector<double>& vecMT,
                                                 bool isClosed=true) {
  //1. Find vector of thickness element
  vector<double> meaningThicknessElement;
  meaningThicknessElement.push_back(vecMT.front());
  for(vector<double>::const_iterator it=vecMT.begin()+1; it!=vecMT.end(); it++) {
    double m=(*it);
    if(std::find(meaningThicknessElement.begin(), meaningThicknessElement.end(),m)==meaningThicknessElement.end())
      meaningThicknessElement.push_back(m);
  }
  std::sort(meaningThicknessElement.begin(),meaningThicknessElement.end());
  for(vector<double>::const_iterator it=meaningThicknessElement.begin(); it!=meaningThicknessElement.end(); it++)
    cout<<"meaningThicknessElement : "<<*it<<endl;
  
  //2. Compute different thickness tangent covers(blurred segments)
  vector<vector<MaximalBlurredSegment3D> > meaningThicknessTangentCover(meaningThicknessElement.size());
  int index=0;
  for(vector<double>::const_iterator it=meaningThicknessElement.begin(); it!=meaningThicknessElement.end(); it++) {
    double thickness=(*it)*sqrt(2);//(*it);//
    vector<Z2i::Point> aContourXY, aContourXZ, aContourYZ;
    bool validOxy=false,validOxz=false,validOyz=false;
    projectionPoints3D(aContour,aContourXY,aContourXZ,aContourYZ,validOxy,validOxz,validOyz);
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=buildTC(aContour,aContourXY,aContourXZ,aContourYZ,thickness);
    cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl;
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++)
      meaningThicknessTangentCover[index].push_back(*it_bis);
    index++;
  }
  
  //3. Update thickness of points with tangent covers
  vector<double> vecMTmodified;
  for(vector<double>::const_iterator it=vecMT.begin(); it!=vecMT.end(); it++)
    vecMTmodified.push_back(*it);
  for(int it=(int)meaningThicknessTangentCover.size()-1; it>=0; it--) {
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=meaningThicknessTangentCover.at(it);//*it;
    double thickness=meaningThicknessElement.at(it);
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++) {
      int idStart=(*it_bis).getIdBegin();
      int idEnd=(*it_bis).getIdEnd();
      if(idStart!=-1 && idEnd!=-1) {
        double maxThickness=statThicknessSegment(*it_bis, vecMT);
        for(int i=idStart; i<=idEnd; i++) {
          if(maxThickness==thickness)
            vecMTmodified.at(i%aContour.size())=maxThickness;
        }
      }
      else
        cout<<"negatif"<<endl;
    }
  }
  
  //4. Travel over the tangent covers and select the segments w.r.t the associated thickness of points
  vector<vector<MaximalBlurredSegment3D> > adaptiveMeaningThicknessTangentCover;
  int idCover=0;
  for(vector<vector<MaximalBlurredSegment3D> >::const_iterator it=meaningThicknessTangentCover.begin(); it!=meaningThicknessTangentCover.end(); it++) {
    adaptiveMeaningThicknessTangentCover.push_back(vector<MaximalBlurredSegment3D>());
    vector<MaximalBlurredSegment3D> fuzzySegmentSet=*it;
    vector<MaximalBlurredSegment3D> AdaptiveFuzzySegmentSet;
    int idSeg=0;
    double thickness=meaningThicknessElement.at(idCover);
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=fuzzySegmentSet.begin();it_bis!=fuzzySegmentSet.end();it_bis++) {
      int idStart=(*it_bis).getIdBegin();
      int idEnd=(*it_bis).getIdEnd();
      if(idStart!=-1 && idEnd!=-1) {
        bool isGoodMTmodif=false, isGoodMT=false;//true
        for(int i=idStart; i<=idEnd; i++) {
          double thicknessMT=vecMT.at(i%aContour.size()); //all elt have same meaningful thickness value(dont contain other meaningful thickness)
          if(thicknessMT==thickness)
            isGoodMT=true;
          double thicknessMTmodif=vecMTmodified.at(i%aContour.size());
          if(thicknessMTmodif==thickness) //there exist at least one elt in modif having meaningful thickness value
            isGoodMTmodif=true;
        }
        if(isGoodMTmodif && isGoodMT)
          AdaptiveFuzzySegmentSet.push_back(*it_bis);
      }
      idSeg++;
    }
    for(vector<MaximalBlurredSegment3D>::const_iterator it_bis=AdaptiveFuzzySegmentSet.begin();it_bis!=AdaptiveFuzzySegmentSet.end();it_bis++)
      adaptiveMeaningThicknessTangentCover[idCover].push_back(*it_bis);
    
    idCover++;
  }
  
  // 4 cont. : selection of Seg
  for(vector<vector<MaximalBlurredSegment3D> >::reverse_iterator it1=adaptiveMeaningThicknessTangentCover.rbegin(); it1!=adaptiveMeaningThicknessTangentCover.rend(); ++it1) {
    vector<MaximalBlurredSegment3D>& segmentSet1=*it1;
    for(vector<vector<MaximalBlurredSegment3D> >::reverse_iterator it2=it1+1; it2!=adaptiveMeaningThicknessTangentCover.rend(); ++it2) {
      vector<MaximalBlurredSegment3D>& segmentSet2=*it2;
      for(vector<MaximalBlurredSegment3D>::iterator itSeg1=segmentSet1.begin();itSeg1!=segmentSet1.end();itSeg1++) {
        //get idStart and idEnd of current seg
        int idCurrentStart=(*itSeg1).getIdBegin();
        int idCurrentEnd=(*itSeg1).getIdEnd();
        
        for(vector<MaximalBlurredSegment3D>::iterator itSeg2=segmentSet2.begin();itSeg2!=segmentSet2.end();itSeg2++) {
          int idStart=(*itSeg2).getIdBegin();
          int idEnd=(*itSeg2).getIdEnd();
          if(idCurrentStart<=idStart && idCurrentEnd>=idEnd) {
            segmentSet2.erase(itSeg2);
            itSeg2--;
          }
        }
      }
    }
  }
  
  //5. Reorder the multi-thickness tangent cover
  vector<MaximalBlurredSegment3D> adaptiveTangentCover;
  int seg=0,nbSeg=0;
  vector<int> idThicknessCover; //stock idSeg of the last seg at idThicknessCover
  for(size_t it=0; it<meaningThicknessElement.size();it++)
    idThicknessCover.push_back(0);
  for(size_t it=0; it<adaptiveMeaningThicknessTangentCover.size(); it++)
    nbSeg+=(adaptiveMeaningThicknessTangentCover.at(it)).size();
  
  while(seg<nbSeg) {
    int idMinStart,idMinEnd;
    if(isClosed) {
      idMinStart=2*aContour.size();
      idMinEnd=2*aContour.size();//x2 : due to the closed curve can have ending index over contour.size()
    }
    else {
      idMinStart=aContour.size();
      idMinEnd=aContour.size();
    }
    int idMin=-1, idSeg=-1;
    //scan adaptiveMeaningThicknessTangentCover
    for(size_t it=0; it<adaptiveMeaningThicknessTangentCover.size(); it++) {//thickness level=it
      //current seg of thickness level idThicknessCover.at(i)
      int idCurrentSeg=idThicknessCover.at(it);
      if(idCurrentSeg<(adaptiveMeaningThicknessTangentCover.at(it)).size()) {
        //get idStart and idEnd of seg
        int idCurrentStart=(adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg).getIdBegin();
        int idCurrentEnd=(adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg).getIdEnd();
        if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size()) {
          //find min idCurrentStart
          if(idMinStart==idCurrentStart && idMinEnd<idCurrentEnd) {
            if(idThicknessCover.at(it)<(int)(adaptiveMeaningThicknessTangentCover.at(it)).size()-1) {
              idThicknessCover.at(idMin)=idThicknessCover.at(idMin)+1;
              seg++;
            }
            idSeg=idCurrentSeg;
            idMin=it;
            idMinStart=idCurrentStart;
            idMinEnd=idCurrentEnd;
          }
          else if(idMinStart>idCurrentStart && idMinEnd>=idCurrentEnd) {
            idSeg=idCurrentSeg;
            idMin=it;
            idMinStart=idCurrentStart;
            idMinEnd=idCurrentEnd;
          }
        }
      }
    }
    adaptiveTangentCover.push_back((adaptiveMeaningThicknessTangentCover.at(idMin)).at(idSeg));
    idThicknessCover.at(idMin)=idThicknessCover.at(idMin)+1;
    seg++;
  }
  return adaptiveTangentCover;
}
/* Compute ATS 3D */

int main(int argc, char** argv)
{
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  std::string inputFile = "../data/sinus_noise.dat";
  string mtDir = "../MeaningfulThickness/";
  
  app.description("Tangential cover for 3D irregular noisy digital curves.\n Example:\n \t ATC3D --input <FileName> --mt <MeaningfulThicknessDir> \n");
  app.add_option("-i,--input,1",inputFile,"Input file.")->required()->check(CLI::ExistingFile);
  app.add_option("-m,--mt",mtDir,"MeaningfulThickness directory for noise detection (default ../MeaningfulThickness/).");
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------
  
  QApplication application(argc,argv);
  vector<Z3i::Point> aContour = PointListReader<Z3i::Point>::getPointsFromFile(inputFile.c_str());
  
  //Project into 2D planes
  vector<Z2i::Point> aContourXY, aContourXZ, aContourYZ;
  bool validOxy=false,validOxz=false,validOyz=false;
  bool bijective=projectionPoints3D(aContour,aContourXY,aContourXZ,aContourYZ,validOxy,validOxz,validOyz);
  
  vector<double> vecMTXYtmp, vecMTXZtmp, vecMTYZtmp, vecMT;
  vecMTXYtmp=getMeaningfulThickness(mtDir,aContourXY,10,1);
  vecMTXZtmp=getMeaningfulThickness(mtDir,aContourXZ,10,1);
  vecMTYZtmp=getMeaningfulThickness(mtDir,aContourYZ,10,1);
  vecMT=getThicknessAssociation(vecMTXYtmp, vecMTXZtmp, vecMTYZtmp, &mid_functor);
  
  //Step 1: get all different noise-levels then sort them
  std::vector<double> vecNoiseLevel;
  vecNoiseLevel.push_back(vecMT.front());
  for(size_t it=1; it<vecMT.size(); it++) {
    if (!(std::find(vecNoiseLevel.begin(), vecNoiseLevel.end(), vecMT.at(it)) != vecNoiseLevel.end()))
      vecNoiseLevel.push_back(vecMT.at(it));
  }
  sort(vecNoiseLevel.begin(), vecNoiseLevel.end());
  std::cout<<"Noise level detected: ";
  for(size_t it=0; it<vecNoiseLevel.size(); it++)
    std::cout<<vecNoiseLevel.at(it)<<" ";
  std::cout<<std::endl;
  
  //Step 2: built ATC
  vector<MaximalBlurredSegment3D> ATC = buildATC(aContour, vecMT);
  vector<int> color;
  vector<AlphaThickSegmentComputer3D> segATC;
  for(size_t it=0; it<ATC.size(); it++) {
    segATC.push_back(ATC.at(it).getSegment());
    color.push_back(int(ATC.at(it).getThickness()/sqrt(2.0)));
  }
  
  /* Visualisation part */
  MyViewer viewer;
  viewer.camera()->setType(Camera::ORTHOGRAPHIC);
  viewer.show();
  
  //show 3D points
  for (size_t i =0; i< aContour.size(); i++) {
    Z3i::Point p (aContour[i][0],aContour[i][1], aContour[i][2]);
    viewer << SetMode3D( p.className(), "Grid" );
    viewer << p;
  }
  viewer << CustomColors3D(Color(250,0,0,100),Color(250,0,0,100));
  for (size_t i =0; i< aContour.size()-1; i++) {
    Z3i::Point p1 (aContour[i][0],aContour[i][1], aContour[i][2]);
    Z3i::Point p2 (aContour[i+1][0],aContour[i+1][1], aContour[i+1][2]);
    viewer.addLine(p1, p2, 3);
  }
  //show projections
  double x,y,z,d=0.5;
  for (size_t i =0; i< aContour.size(); i++) {
    viewer << CustomColors3D( Color::Blue,  Color::Blue );
    x=aContourXY[i][0]-d;y=aContourXY[i][1]-d;
    viewer.addQuad(Z3i::RealPoint(x,y,0),Z3i::RealPoint(x+1,y,0),Z3i::RealPoint(x+1,y+1,0),Z3i::RealPoint(x,y+1,0));
    viewer << CustomColors3D( Color::Red,  Color::Red );
    x=aContourXZ[i][0]-d;z=aContourXZ[i][1]-d;
    viewer.addQuad(Z3i::RealPoint(x,0,z),Z3i::RealPoint(x+1,0,z),Z3i::RealPoint(x+1,0,z+1),Z3i::RealPoint(x,0,z+1));
    viewer << CustomColors3D( Color::Green,  Color::Green );
    y=aContourYZ[i][0]-d;z=aContourYZ[i][1]-d;
    viewer.addQuad(Z3i::RealPoint(0,y,z),Z3i::RealPoint(0,y+1,z),Z3i::RealPoint(0,y+1,z+1),Z3i::RealPoint(0,y,z+1));
  }
  viewer << CustomColors3D( Color::Black,  Color::Black );
  for (size_t i =0; i< aContour.size()-1; i++) {
    x=aContourXY[i][0];y=aContourXY[i][1];
    Z3i::Point p1(x,y,0);
    x=aContourXY[i+1][0];y=aContourXY[i+1][1];
    Z3i::Point p2(x,y,0);
    viewer.addLine(p1, p2);
    
    x=aContourXZ[i][0];z=aContourXZ[i][1];
    p1 = Z3i::Point(x,0,z);
    x=aContourXZ[i+1][0];z=aContourXZ[i+1][1];
    p2 = Z3i::Point(x,0,z);
    viewer.addLine(p1, p2);
    
    y=aContourYZ[i][0];z=aContourYZ[i][1];
    p1 = Z3i::Point(0,y,z);
    y=aContourYZ[i+1][0];z=aContourYZ[i+1][1];
    p2 = Z3i::Point(0,y,z);
    viewer.addLine(p1, p2);
  }
  
  // Draw ATC 3D
  viewer << SetMode3D(segATC.front().className(), "BoundingBox");
  for(size_t it=0; it<segATC.size(); it++) {
    switch (color.at(it)) {
      case 1:
        viewer << CustomColors3D(DGtal::Color::Purple,DGtal::Color::Purple);
        break;
      case 2:
        viewer << CustomColors3D(DGtal::Color::Cyan,DGtal::Color::Cyan);
        break;
      case 3:
        viewer << CustomColors3D(DGtal::Color::Yellow,DGtal::Color::Yellow);
        break;
      case 4:
        viewer << CustomColors3D(DGtal::Color::Black,DGtal::Color::Black);
        break;
      case 5:
        viewer << CustomColors3D(DGtal::Color::Green,DGtal::Color::Green);
        break;
      case 6:
        viewer << CustomColors3D(DGtal::Color::Blue,DGtal::Color::Blue);
        break;
      default:
        viewer << CustomColors3D(DGtal::Color::Magenta,DGtal::Color::Magenta);
        break;
    }
    viewer << segATC.at(it);
  }
  viewer << MyViewer::updateDisplay;
  return application.exec();
  return 0;
}

