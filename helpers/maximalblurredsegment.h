#ifndef MAXIMALBLURREDSEGMENT_H
#define MAXIMALBLURREDSEGMENT_H

#include "GeomHelpers.h"

class MaximalBlurredSegment
{
private:
  int idBegin;
  int idEnd;
  double thickness;
  AlphaThickSegmentComputer2D segment;
  
public:
  MaximalBlurredSegment();
  MaximalBlurredSegment(int idB, int idE, double t, const AlphaThickSegmentComputer2D& seg);
  
  const int getIdBegin() const {return idBegin;}
  const int getIdEnd() const {return idEnd;}
  const double getThickness() const {return thickness;}
  const AlphaThickSegmentComputer2D getSegment() const {return segment;}
  
  void setIdBegin(int id);
  void setIdEnd(int id);
  void setSegment(const AlphaThickSegmentComputer2D& seg);
};

#endif // MAXIMALBLURREDSEGMENT_H
