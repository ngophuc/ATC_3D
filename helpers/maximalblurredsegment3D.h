#ifndef MAXIMALBLURREDSEGMENT3D_H
#define MAXIMALBLURREDSEGMENT3D_H

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/AlphaThickSeg3DComputer.h"

using namespace std;
using namespace DGtal;

typedef std::vector<Z3i::Point>::iterator Iterator;
typedef std::vector<Z3i::Point>::const_iterator ConstIterator;
typedef AlphaThickSeg3DComputer<Iterator,int> AlphaThickSegmentComputer3D;


class MaximalBlurredSegment3D
{
private:
  int idBegin;
  int idEnd;
  double thickness;
  AlphaThickSegmentComputer3D segment;
  
public:
  MaximalBlurredSegment3D();
  MaximalBlurredSegment3D(int idB, int idE, double t, const AlphaThickSegmentComputer3D& seg);
  
  const int getIdBegin() const {return idBegin;}
  const int getIdEnd() const {return idEnd;}
  const double getThickness() const {return thickness;}
  const AlphaThickSegmentComputer3D getSegment() const {return segment;}
  
  void setIdBegin(int id);
  void setIdEnd(int id);
  void setThickness(int t);
  void setSegment(const AlphaThickSegmentComputer3D& seg);
};
  
#endif // MAXIMALBLURREDSEGMENT3D_H
