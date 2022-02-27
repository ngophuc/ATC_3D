#include "maximalblurredsegment.h"

MaximalBlurredSegment::MaximalBlurredSegment()
{
  idBegin=-1;
  idEnd=-1;
}

MaximalBlurredSegment::MaximalBlurredSegment(int idB, int idE, double t, const AlphaThickSegmentComputer2D& seg)
{
  idBegin=idB;
  idEnd=idE;
  thickness=t;
  segment=seg;
}

void MaximalBlurredSegment::setIdBegin(int id)
{
  idBegin=id;
}

void MaximalBlurredSegment::setIdEnd(int id)
{
  idEnd=id;
}

void MaximalBlurredSegment::setSegment(const AlphaThickSegmentComputer2D& seg)
{
  segment=seg;
}
