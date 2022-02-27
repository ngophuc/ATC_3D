#include "maximalblurredsegment3D.h"

MaximalBlurredSegment3D::MaximalBlurredSegment3D()
{
  idBegin=-1;
  idEnd=-1;
  thickness = 0;
}

MaximalBlurredSegment3D::MaximalBlurredSegment3D(int idB, int idE, double t, const AlphaThickSegmentComputer3D& seg)
{
  idBegin=idB;
  idEnd=idE;
  thickness=t;
  segment=seg;
}


void MaximalBlurredSegment3D::setIdBegin(int id)
{
  idBegin=id;
}

void MaximalBlurredSegment3D::setIdEnd(int id)
{
  idEnd=id;
}

void MaximalBlurredSegment3D::setThickness(int t)
{
  thickness = t;
}

void MaximalBlurredSegment3D::setSegment(const AlphaThickSegmentComputer3D& seg)
{
  segment=seg;
}
