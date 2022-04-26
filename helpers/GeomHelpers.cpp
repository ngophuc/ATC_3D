#include "GeomHelpers.h"

bool isInConvexPolygon(Z2i::Point testPoint, std::vector<Z2i::Point> polygon) {
  //Check if a triangle or higher n-gon
  if(polygon.size() <= 2) return true;
  
  //n>2 Keep track of cross product sign changes
  int pos = 0;
  int neg = 0;
  
  for (size_t it = 0; it < polygon.size(); it++)
  {
    //If point is in the polygon
    if (polygon.at(it) == testPoint)
      return true;
    
    //Form a segment between the i'th point
    auto x1 = polygon.at(it)[0];
    auto y1 = polygon.at(it)[1];
    
    //And the i+1'th, or if i is the last, with the first point
    int it2 = (it+1)%polygon.size();
    
    auto x2 = polygon.at(it2)[0];
    auto y2 = polygon.at(it2)[1];
    
    auto x = testPoint[0];
    auto y = testPoint[1];
    
    //Compute the cross product
    auto d = (x - x1)*(y2 - y1) - (y - y1)*(x2 - x1);
    
    if (d > 0) pos++;
    if (d < 0) neg++;
    
    //If the sign changes, then point is outside
    if (pos > 0 && neg > 0)
      return false;
  }
  
  //If no change in direction, then on same side of all segments, and thus inside
  return true;
}

void writePointList(const vector<Z3i::Point> &v, const char *filename)
{
  ofstream myfile (filename);
  if (myfile.is_open()) {
    for(vector<Z3i::Point>::const_iterator it=v.begin(); it != v.end(); it++)
      myfile <<(*it)[0]<<" "<<(*it)[1]<<" "<<(*it)[2]<<endl;
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
}

void writeRealPointList(const vector<Z3i::RealPoint> &v, const char *filename)
{
  ofstream myfile (filename);
  if (myfile.is_open()) {
    for(vector<Z3i::RealPoint>::const_iterator it=v.begin(); it != v.end(); it++)
      myfile <<(*it)[0]<<" "<<(*it)[1]<<" "<<(*it)[2]<<endl;
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
}

void writeVectorList(const vector<double> &v, const char *filename)
{
  ofstream myfile (filename);
  if (myfile.is_open()) {
    for(vector<double>::const_iterator it=v.begin(); it != v.end(); it++)
      myfile <<(*it)<<endl;
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
}

vector<Z3i::Point> readPointList(const char* filename)
{
  ifstream myfile (filename);
  vector<Z3i::Point> vecP;
  if (myfile.is_open()) {
    int x, y, z;
    while (myfile >> x >> y >> z) {
      vecP.push_back(Z3i::Point(x,y,z));
    }
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
  return vecP;
}

vector<Z3i::RealPoint> readRealPointList(const char* filename) {
  ifstream myfile (filename);
  vector<Z3i::RealPoint> vecP;
  if (myfile.is_open()) {
    double x, y, z;
    while (myfile >> x >> y >> z) {
      vecP.push_back(Z3i::RealPoint(x,y,z));
    }
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
  return vecP;
}
vector<double> readVectorList(const char *filename) {
  ifstream myfile (filename);
  vector<double> vecV;
  if (myfile.is_open()) {
    double x;
    while (myfile >> x) {
      vecV.push_back(x);
    }
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
  return vecV;
}

Color getColor(int i)
{
  switch (i) {
    case 1: return Color::Red;
    case 2: return Color::Green;
    case 3: return Color::Blue;
    case 4: return Color::Magenta;
    case 5: return Color::Yellow;
    case 6: return Color::Navy;
    case 7: return Color::Cyan;
    case 8: return Color::Lime;
    case 9: return Color::Aqua;
    default: return Color::Black;
  }
}

bool isRepeatedVector(const vector<Z2i::Point>& vec)
{
  for(vector<Z2i::Point>::const_iterator it=vec.begin(); it<vec.end()-1;it++) {
    vector<Z2i::Point>::const_iterator a=find(it+1,vec.end(),*it);
    if(a != vec.end() && a != it+1)
      return true;
  }
  return false;
}

bool isRepeatedPoint(const vector<Z2i::Point>& vec, Z2i::Point p)
{
  vector<Z2i::Point>::const_iterator a=find(vec.begin(),vec.end(),p);
  if(a != vec.end())
    return true;
  return false;
}

bool isValidPoints(const vector<Z2i::Point>& pointsOxy, const vector<Z2i::Point>& pointOxz, const vector<Z2i::Point>& pointOyz, size_t i)
{
  unsigned int lenXY = std::distance ( pointsOxy.begin(), pointsOxy.end() );
  unsigned int lenXZ = std::distance ( pointOxz.begin(), pointOxz.end() );
  unsigned int lenYZ = std::distance ( pointOyz.begin(), pointOyz.end() );
  if ( i == 0 && ( lenYZ >= lenXZ || lenYZ >= lenXY ) )
    return true;
  else if ( i == 1 && ( lenXZ >= lenXY || lenXZ >= lenYZ ) )
    return true;
  else if ( i == 2 && ( lenXY >= lenXZ || lenXY >= lenYZ ) )
    return true;
  return false;
}

bool projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz, vector<std::tuple<bool, bool, bool> >& projection)
{
  bool xy=true, xz=true, yz=true;
  for(vector<Z3i::Point>::const_iterator it=aContour.begin(); it!=aContour.end(); it++) {
    Z2i::Point pxy((*it)[0],(*it)[1]);
    Z2i::Point pxz((*it)[0],(*it)[2]);
    Z2i::Point pyz((*it)[1],(*it)[2]);
    xy=false, xz=false, yz=false;
    pointOxy.push_back(pxy);
    if(!isRepeatedPoint(pointOxy,pxy))
      xy=true;
    pointOxz.push_back(pxz);
    if(!isRepeatedPoint(pointOxz,pxz))
      xz=true;
    pointOyz.push_back(pyz);
    if(!isRepeatedPoint(pointOyz,pyz))
      yz=true;
    projection.push_back(make_tuple(xz,xz,yz));
  }
  
  if((!isRepeatedVector(pointOxy) && !isRepeatedVector(pointOxz)) || (!isRepeatedVector(pointOxy) && !isRepeatedVector(pointOyz)) || (!isRepeatedVector(pointOyz) && !isRepeatedVector(pointOxz)))
    return true;
  return false;
}

bool projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz,
                        bool& validOxy, bool& validOxz, bool& validOyz)
{
  for(vector<Z3i::Point>::const_iterator it=aContour.begin(); it!=aContour.end(); it++) {
    pointOxy.push_back(Z2i::Point((*it)[0],(*it)[1]));
    pointOxz.push_back(Z2i::Point((*it)[0],(*it)[2]));
    pointOyz.push_back(Z2i::Point((*it)[1],(*it)[2]));
  }
  validOxy=true;
  validOxz=true;
  validOyz=true;
  //verify of bijectif projetion
  if(isRepeatedVector(pointOxy))//isValidPoints(pointOxy,pointOxz,pointOyz,0)
    validOxy=false;//cout<<"Oxy is not bijectif"<<endl;
  if(isRepeatedVector(pointOxz))//isValidPoints(pointOxy,pointOxz,pointOyz,1)
    validOxz=false;//cout<<"Oxz is not bijectif"<<endl;
  if(isRepeatedVector(pointOyz))//(isValidPoints(pointOxy,pointOxz,pointOyz,2))
    validOyz=false;//cout<<"Oyz is not bijectif"<<endl;
  if((!isRepeatedVector(pointOxy) && !isRepeatedVector(pointOxz)) || (!isRepeatedVector(pointOxy) && !isRepeatedVector(pointOyz)) || (!isRepeatedVector(pointOyz) && !isRepeatedVector(pointOxz)))
    return true;
  return false;
}

void projectionPoints3D(const vector<Z3i::Point>& aContour, vector<Z2i::Point>& pointOxy, vector<Z2i::Point>& pointOxz, vector<Z2i::Point>& pointOyz)
{
  for(vector<Z3i::Point>::const_iterator it=aContour.begin(); it!=aContour.end(); it++) {
    pointOxy.push_back(Z2i::Point((*it)[0],(*it)[1]));
    pointOxz.push_back(Z2i::Point((*it)[0],(*it)[2]));
    pointOyz.push_back(Z2i::Point((*it)[1],(*it)[2]));
  }
}

vector<double> readMeanindfulThicknessFile(const char* filename)
{
  vector<double> P;
  ifstream myfile (filename);
  string lineTmp;
  if (myfile.is_open()) {
    int count=0;
    //idx noiselvl code x y
    double m=0,x=0,y=0;
    std::getline(myfile, lineTmp);//ignore the first two lines
    std::getline(myfile, lineTmp);
    while(myfile >> x >> y >> m)// X Y noiseLevel
    {
      P.push_back(2.0*m);//*sqrt(2)
      count++;
    }
    myfile.close();
    cout<<count<<" points are read "<<endl;
  }
  else cout << "Unable to open file " <<filename;
  
  return P;
}

void writeFile(const vector<Z2i::Point>& v, const char* filename)
{
  ofstream myfile (filename);
  if (myfile.is_open())
  {
    for(vector<Z2i::Point>::const_iterator it=v.begin(); it != v.end(); it++)
      myfile <<(*it)[0]<<" "<<(*it)[1]<<endl;
    myfile.close();
  }
  else cout << "Unable to open file " <<filename;
}

vector<double> getMeaningfulThickness(const string mtDir, const vector<Z2i::Point> aContour, double maxMT, double stepMT)
{
  cout<<"getMeaningfulThickness"<<endl;
  string contourFile,noiseLevelMTFile;
  noiseLevelMTFile="MeanThickness.txt";
  contourFile="Contour.txt";
  writeFile(aContour,contourFile.c_str());
  std::stringstream instruction;
  instruction << mtDir << "/build/MeaningfulThickness"
  << " -srcPolygon " << contourFile
  << " 0 1 -setSampling "<<maxMT<<" "<<stepMT
  << " -exportNoiseLevel "<< noiseLevelMTFile;
  std::system(instruction.str().c_str());
  return readMeanindfulThicknessFile(noiseLevelMTFile.c_str());
}

vector<double> sortMeaningfulThickness(const vector<double> vecMT)
{
  vector<double> vecThickness;
  vecThickness.push_back(vecMT.front());
  for(vector<double>::const_iterator it = vecMT.begin()+1; it != vecMT.end(); it++)
  {
    double m = (*it);
    if (std::find(vecThickness.begin(), vecThickness.end(),m)==vecThickness.end())
      vecThickness.push_back(m);
  }
  std::sort(vecThickness.begin(),vecThickness.end(),[](double a, double b) {return a < b;});
  return vecThickness;
}

Z2i::Point getStartPoint(const AlphaThickSegmentComputer2D s)
{
  if(s.size()==0)
    return *(s.begin()); //add point from iterator
  return *(s.containerBegin());
}

Z2i::Point getEndPoint(const AlphaThickSegmentComputer2D s)
{
  Z2i::Point p;
  if(s.size()==0)
    for(vector<Z2i::Point>::const_iterator it=s.begin();it != s.end();it++)
      p=*(it);
  else
    for(vector<Z2i::Point>::const_iterator it=s.containerBegin();it != s.containerEnd();it++)
      p=*(it);
  return p;
}

int findElement(const vector<Z2i::Point>& vec, Z2i::Point p)
{
  vector<Z2i::Point>::const_iterator it=find(vec.begin(),vec.end(),p);
  if(it != vec.end())
    return (it-vec.begin());
  return -1;
}

int findElement(const vector<Z2i::Point>& vec, Z2i::Point p, int start)
{
  int it;
  for(it=start; it<vec.size(); it++)
    if(vec.at(it)==p)
      return it;
  for(it=0; it<start; it++)
    if(vec.at(it)==p)
      return (int)vec.size()+it;
  return -1;
}

vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecomposition(const vector<Z2i::Point>& aContour, double thickness)
{
  vector<AlphaThickSegmentComputer2D> fuzzySegmentSet;
  int idBegin=0, idEnd=0;
  for(vector<Z2i::Point>::const_iterator it=aContour.begin();it!=aContour.end();it++) {
    idEnd=idBegin;
    AlphaThickSegmentComputer2D aSegment(thickness);
    aSegment.init(it);
    while(idBegin<aContour.size() && idEnd<aContour.size() && aSegment.extendFront(aContour.at(idEnd%aContour.size()))){idEnd++;}
    if(it==aContour.begin())
      fuzzySegmentSet.push_back(aSegment);
    else if(findElement(aContour,getEndPoint(aSegment))>findElement(aContour,getEndPoint(fuzzySegmentSet.back())))
      fuzzySegmentSet.push_back(aSegment);
    if(getEndPoint(aSegment)==aContour.back() || it==aContour.end() || idBegin>=aContour.size())
      break;
    idBegin++;
  }
  return fuzzySegmentSet;
}

vector<double> labelMeaningfulThickness(const vector<Z2i::Point>& aContour, const vector<double>& vecMT)
{
  if(vecMT.size()==1)
    return vecMT;
  //1. Find vector of thickness element
  vector<double> meaningThicknessElement=sortMeaningfulThickness(vecMT);
  
  //2. Compute different thickness tangent covers (blurred segments)
  vector<vector<AlphaThickSegmentComputer2D> > meaningThicknessTangentCover(meaningThicknessElement.size());
  int index = 0;
  for(vector<double>::const_iterator it = meaningThicknessElement.begin(); it != meaningThicknessElement.end(); it++) {
    double thickness = (*it)*sqrt(2);
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = blurredSegmentCurveDecomposition(aContour,thickness);
    cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl;
    for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
      meaningThicknessTangentCover[index].push_back(*it_bis);
    index++;
  }
  //3. Update thickness of points with tangent covers
  vector<double> vecMTmodified;
  for(vector<double>::const_iterator it = vecMT.begin(); it != vecMT.end(); it++)
    vecMTmodified.push_back(*it);
  for(int it=meaningThicknessTangentCover.size()-1; it>=0; it--) {
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = meaningThicknessTangentCover.at(it);//*it;
    double thickness = meaningThicknessElement.at(it);
    for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++) {
      int idStart = findElement(aContour,getStartPoint(*it_bis));
      int idEnd = findElement(aContour,getEndPoint(*it_bis),idStart);
      if(idStart != -1 && idEnd != -1) {
        double maxThickness = -1;
        for(int i=idStart; i<=idEnd; i++) {
          double thicknessPoint = vecMT.at(i);
          if(thicknessPoint > maxThickness)
            maxThickness = thicknessPoint;
        }
        for(int i=idStart; i<=idEnd; i++) {
          if(maxThickness==thickness)
            vecMTmodified.at(i) = maxThickness;
        }
      }
      else
        cout<<"negatif"<<endl;
    }
  }
  return vecMTmodified;
}

void getParameters (const AlphaThickSegmentComputer2D &mySeg_Pl1, const AlphaThickSegmentComputer2D &mySeg_Pl2,
                    const bool validOxy, const bool validOxz, const bool validOyz,
                    Z3i::RealPoint& direction, Z3i::RealPoint& intercept, Z3i::RealPoint& thickness)
{
  if(validOxy && validOxz)//XY-plane, XZ-plane
  {
    double a1 = -mySeg_Pl1.getNormal()[1];//myXYalgo.b();
    double b1 = mySeg_Pl1.getNormal()[0];//myXYalgo.a();
    double a2 = -mySeg_Pl2.getNormal()[1];//myXZalgo.b();
    double c1 = mySeg_Pl2.getNormal()[0];//myXZalgo.a();
    
    if ( c1 == 0 || ( a1 == 0 && a2 == 0 ) )
      direction = Z3i::RealPoint ( a1, b1, c1 );
    else
    {
      if ( b1 == 0  )
        direction = Z3i::RealPoint ( a2, b1, c1 );
      else
        direction = Z3i::RealPoint ( a1 * a2 , a2 * b1 , a1 * c1 );
    }
    double mu1 = mySeg_Pl1.getMu();//myXYalgo.mu();
    double mu2 = mySeg_Pl2.getMu();//myXZalgo.mu();
    intercept[0] = 0; intercept[1] = -mu1/ a1; intercept[2] = -mu2/a2;
    double omega1 = mySeg_Pl1.getNu();//myXYalgo.omega()-1;
    double omega2 = mySeg_Pl2.getNu();//myXZalgo.omega()-1;
    thickness[0] = 0; thickness[1] = -omega1/a1; thickness[2] = -omega2/a2;
  }
  else if(validOxy && validOyz)//XY-plane, YZ-plane
  {
    double a1 = -mySeg_Pl1.getNormal()[1];//myXYalgo.b();
    double b1 = mySeg_Pl1.getNormal()[0];//myXYalgo.a();
    double b2 = -mySeg_Pl2.getNormal()[1];//myXZalgo.b();
    double c2 = mySeg_Pl2.getNormal()[0];//myXZalgo.a();
    if ( a1 == 0 || ( b2 == 0 && b1 == 0 ) )
      direction = Z3i::RealPoint ( a1, b2, c2 );
    else
    {
      if ( c2 == 0 )
        direction = Z3i::RealPoint ( a1, b1, c2 );
      else
        direction = Z3i::RealPoint ( b2 * a1 , b1 * b2 , b1 * c2 );
    }
    double mu1 = mySeg_Pl1.getMu();//myXYalgo.mu();
    double mu2 = mySeg_Pl2.getMu();//myXZalgo.mu();
    intercept[0] = mu1/b1; intercept[1] = 0; intercept[2] = -mu2/b2;
    double omega1 = mySeg_Pl1.getNu();//myXYalgo.omega()-1;
    double omega2 = mySeg_Pl2.getNu();//myXZalgo.omega()-1;
    thickness[0] = omega1/a1; thickness[1] = 0; thickness[2] = -omega2/b2;
  }
  else //YZ-plane, XZ-plane
  {
    double b2 = -mySeg_Pl2.getNormal()[1];//myXYalgo.b();
    double c2 = mySeg_Pl2.getNormal()[0];//myXYalgo.a();
    double a2 = -mySeg_Pl1.getNormal()[1];//myXZalgo.b();
    double c1 = mySeg_Pl1.getNormal()[0];//myXZalgo.a();
    if ( a2 == 0 || ( c2 == 0 && c1 == 0 ) )
      direction = Z3i::RealPoint ( a2, b2, c2 );
    else
    {
      if ( b2 == 0 )
        direction = Z3i::RealPoint ( a2, b2, c1 );
      else
        direction = Z3i::RealPoint ( c2 * a2, c1 * b2, c1 * c2 );
    }
    double mu1 = mySeg_Pl2.getMu();//myXYalgo.mu();
    double mu2 = mySeg_Pl1.getMu();//myXZalgo.mu();
    intercept[0] = mu2/c1; intercept[1] = mu1/c2; intercept[2] = 0;
    double omega1 = mySeg_Pl2.getNu();//myXYalgo.omega()-1;
    double omega2 = mySeg_Pl1.getNu();//myXZalgo.omega()-1;
    thickness[0] = omega2/c1; thickness[1] = omega1/c2; thickness[2] = 0;
  }
}

void drawAsBoundingBox(MyViewer& display,
                       const vector<Z3i::Point>& aContour,
                       size_t idBegin, size_t idEnd,
                       const AlphaThickSegmentComputer2D& mySeg_Pl1,
                       const AlphaThickSegmentComputer2D& mySeg_Pl2,
                       const bool validOxy, const bool validOxz, const bool validOyz)
{
  //get DSS parameters
  Z3i::RealPoint u; //direction vector
  Z3i::RealPoint mu; //intercept
  Z3i::RealPoint omega; //thickness
  getParameters(mySeg_Pl1, mySeg_Pl2, validOxy, validOxz, validOyz,u, mu, omega);
  //Squared L2 norm of u
  double n = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
  
  //first and last points
  Z3i::Point first = aContour.at(idBegin);
  Z3i::Point last = aContour.at(idEnd);
  Z3i::RealPoint f = Z3i::RealPoint(NumberTraits<int>::castToDouble(first[0]),
                        NumberTraits<int>::castToDouble(first[1]),
                        NumberTraits<int>::castToDouble(first[2]) );
  Z3i::RealPoint l = Z3i::RealPoint(NumberTraits<int>::castToDouble(last[0]),
                        NumberTraits<int>::castToDouble(last[1]),
                        NumberTraits<int>::castToDouble(last[2]) );
  
  if (n != 0) {
    //last coefficient of the normal plane to the DSS direction
    //passing trough f and l
    double df = -u[0]*f[0] -u[1]*f[1] -u[2]*f[2];
    double dl = -u[0]*l[0] -u[1]*l[1] -u[2]*l[2];
    
    //omega masks
    Z3i::RealPoint omega1, omega2;
    if (omega[0] == 0) {
      omega1 = Z3i::RealPoint(0,omega[1],0);
      omega2 = Z3i::RealPoint(0,0,omega[2]);
    } else {
      if (omega[1] == 0) {
        omega1 = Z3i::RealPoint(omega[0],0,0);
        omega2 = Z3i::RealPoint(0,0,omega[2]);
      } else {//omega[2] == 0
        omega1 = Z3i::RealPoint(omega[0],0,0);
        omega2 = Z3i::RealPoint(0,omega[1],0);
      }
    }
    
    double m1 = u[0]*mu[0] + u[1]*mu[1] + u[2]*mu[2];
    double m2 = u[0]*(mu[0]+omega1[0]) + u[1]*(mu[1]+omega1[1]) + u[2]*(mu[2]+omega1[2]);
    double m3 = u[0]*(mu[0]+omega2[0]) + u[1]*(mu[1]+omega2[1]) + u[2]*(mu[2]+omega2[2]);
    double m4 = u[0]*(mu[0]+omega[0]) + u[1]*(mu[1]+omega[1]) + u[2]*(mu[2]+omega[2]);
    
    //4 lines
    Z3i::RealPoint pf = Z3i::RealPoint( mu[0] - ( (m1+df)*u[0] )/n,
                           mu[1] - ( (m1+df)*u[1] )/n,
                           mu[2] - ( (m1+df)*u[2] )/n );
    Z3i::RealPoint pl = Z3i::RealPoint( mu[0] - ( (m1+dl)*u[0] )/n,
                           mu[1] - ( (m1+dl)*u[1] )/n,
                           mu[2] - ( (m1+dl)*u[2] )/n );
    
    display.addLine(pf, pl);
    
    Z3i::RealPoint pf2 = Z3i::RealPoint((mu[0]+omega1[0]) - ( (m2+df)*u[0] )/n,
                            (mu[1]+omega1[1]) - ( (m2+df)*u[1] )/n,
                            (mu[2]+omega1[2]) - ( (m2+df)*u[2] )/n );
    Z3i::RealPoint pl2 = Z3i::RealPoint((mu[0]+omega1[0]) - ( (m2+dl)*u[0] )/n,
                            (mu[1]+omega1[1]) - ( (m2+dl)*u[1] )/n,
                            (mu[2]+omega1[2]) - ( (m2+dl)*u[2] )/n );
    display.addLine(pf2, pl2);
    
    Z3i::RealPoint pf3 = Z3i::RealPoint((mu[0]+omega2[0]) - ( (m3+df)*u[0] )/n,
                            (mu[1]+omega2[1]) - ( (m3+df)*u[1] )/n,
                            (mu[2]+omega2[2]) - ( (m3+df)*u[2] )/n );
    Z3i::RealPoint pl3 = Z3i::RealPoint((mu[0]+omega2[0]) - ( (m3+dl)*u[0] )/n,
                            (mu[1]+omega2[1]) - ( (m3+dl)*u[1] )/n,
                            (mu[2]+omega2[2]) - ( (m3+dl)*u[2] )/n );
    display.addLine(pf3, pl3);
    
    Z3i::RealPoint pf4 = Z3i::RealPoint((mu[0]+omega[0]) - ( (m4+df)*u[0] )/n,
                            (mu[1]+omega[1]) - ( (m4+df)*u[1] )/n,
                            (mu[2]+omega[2]) - ( (m4+df)*u[2] )/n );
    Z3i::RealPoint pl4 = Z3i::RealPoint((mu[0]+omega[0]) - ( (m4+dl)*u[0] )/n,
                            (mu[1]+omega[1]) - ( (m4+dl)*u[1] )/n,
                            (mu[2]+omega[2]) - ( (m4+dl)*u[2] )/n );
    display.addLine(pf4, pl4);
    
    //two end facets
    display.addLine(pf, pf2);
    display.addLine(pf2,pf4);
    display.addLine(pf4, pf3);
    display.addLine(pf3, pf);
    
    display.addLine(pl, pl2);
    display.addLine(pl2, pl4);
    display.addLine(pl4, pl3);
    display.addLine(pl3, pl);
  }
}

void getTangentCoverDifferentNoise(const vector<double>& vecThickness,
                                   const vector<Z2i::Point>& aContour_Pl1,
                                   const vector<Z2i::Point>& aContour_Pl2,
                                   vector<vector<AlphaThickSegmentComputer2D> >& vecTC_Pl1,
                                   vector<vector<AlphaThickSegmentComputer2D> >& vecTC_Pl2)
{
  for(size_t it=0; it<vecThickness.size(); it++) {
    double thickness=vecThickness.at(it);
    vector<AlphaThickSegmentComputer2D> TC_Pl1,TC_Pl2;
    /* Decomposition of curves into MBS with two projection-planes */
    size_t idBegin_Pl1=0, idEnd_Pl1=0;
    size_t idBegin_Pl2=0, idEnd_Pl2=0;
    size_t idEnd2D_prev=0;
    vector<Z2i::Point>::const_iterator it1=aContour_Pl1.begin();
    vector<Z2i::Point>::const_iterator it2=aContour_Pl2.begin();
    
    AlphaThickSegmentComputer2D aSegment_Pl1(thickness), aSegment_Pl2(thickness);
    for(;it1!=aContour_Pl1.end()&&it2!=aContour_Pl2.end();it1++,it2++)
    {
      idEnd_Pl1=idBegin_Pl1;
      idEnd_Pl2=idBegin_Pl1;
      aSegment_Pl1.init(it1);
      aSegment_Pl2.init(it2);
      //while (aSegment_Pl1.end()!=aContour_Pl1.end() && aSegment_Pl1.extendFront()){idEnd_Pl1++;}
      //while (aSegment_Pl2.end()!=aContour_Pl2.end() && aSegment_Pl2.extendFront()){idEnd_Pl2++;}
      while (aSegment_Pl1.end()!=aContour_Pl1.end() && aSegment_Pl1.extendFront() &&
             aSegment_Pl2.end()!=aContour_Pl2.end() && aSegment_Pl2.extendFront()){idEnd_Pl1++;idEnd_Pl2++;}
      idEnd_Pl1 = idEnd_Pl1 < idEnd_Pl2 ? idEnd_Pl1 : idEnd_Pl2;
      idEnd_Pl2 = idEnd_Pl1;
      if(idEnd_Pl1!=idEnd2D_prev)
      {
        TC_Pl1.push_back(aSegment_Pl1);
        TC_Pl2.push_back(aSegment_Pl2);
      }
      idBegin_Pl1++;
      idBegin_Pl2++;
      idEnd2D_prev=idEnd_Pl1;
    }
    /* Decomposition of curves into MBS with two projection-planes */
    cout<<"Thickness="<<thickness<<", TC_Pl1 size="<<TC_Pl1.size()<<endl;
    cout<<"Thickness="<<thickness<<", TC_Pl2 size="<<TC_Pl2.size()<<endl;
    vecTC_Pl1.push_back(TC_Pl1);
    vecTC_Pl2.push_back(TC_Pl2);
  }
  cout<<"vecTC_Pl1 size="<<vecTC_Pl1.size()<<endl;
  cout<<"vecTC_Pl2 size="<<vecTC_Pl2.size()<<endl;
}

void getSegmentsAdaptiveTangentCover(const vector<Z3i::Point>& aContour,
                                     const vector<double>& vecThickness,
                                     const vector<double>& vecMT3D,
                                     const vector<vector<AlphaThickSegmentComputer2D> >& vecTC_Pl1,
                                     const vector<vector<AlphaThickSegmentComputer2D> >&vecTC_Pl2,
                                     vector<vector<AlphaThickSegmentComputer2D>>& vecATC_Pl1,
                                     vector<vector<AlphaThickSegmentComputer2D>>& vecATC_Pl2,
                                     const bool validOxy, const bool validOxz, const bool validOyz)
{
  vector<vector<AlphaThickSegmentComputer2D> > vecATC_Pl1tmp,vecATC_Pl2tmp;
  vector<vector<size_t> > vecBeginTmp,vecEndTmp;
  /* Select seg for ATC */
  vector<vector<AlphaThickSegmentComputer2D>>::const_iterator it1=vecTC_Pl1.begin(),it2=vecTC_Pl2.begin();
  size_t noiseLevel=0;
  for(;it1!=vecTC_Pl1.end() && it2!=vecTC_Pl2.end();it1++,it2++)
  {
    //cout<<"=====> noiseLevel="<<vecThickness.at(noiseLevel)<<endl;
    vector<AlphaThickSegmentComputer2D> ATC_Pl1,ATC_Pl2;
    vector<size_t> vecBegin,vecEnd;
    //for each noise level, keep the segment having at least one point label = MT
    vector<AlphaThickSegmentComputer2D>::const_iterator it1_bis=(*it1).begin(),it2_bis=(*it2).begin();
    for(;it1_bis!=(*it1).end() && it2_bis!=(*it2).end();it1_bis++,it2_bis++)
    {
      Z2i::Point pl1_begin=*((*it1_bis).begin());
      Z2i::Point pl2_begin=*((*it2_bis).begin());
      Z2i::Point pl1_end=*((*it1_bis).end()-1);
      Z2i::Point pl2_end=*((*it2_bis).end()-1);
      //cout<<"pl1_begin="<<pl1_begin<<" and pl1_end="<<pl1_end<<endl;
      //cout<<"pl2_begin="<<pl2_begin<<" and pl2_end="<<pl2_end<<endl;
      
      size_t idBegin, idEnd;
      if(validOxy && validOxz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl1_begin[1],pl2_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl1_end[1],pl2_end[1]));
        idEnd= itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      else if(validOxy && validOyz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl1_begin[1],pl2_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl1_end[1],pl2_end[1]));
        idEnd= itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      else//if(validOxz && validOyz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl2_begin[0],pl1_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl2_end[0],pl1_end[1]));
        idEnd= itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      
      if(idBegin==-1 || idEnd==-1)
        cout<<"ERROR....idBegin="<<idBegin<<" and idEnd="<<idEnd
        <<"Pbegin="<<Z3i::Point(pl1_begin[0],pl1_begin[1],pl2_begin[1])
        <<"and Pend="<<Z3i::Point(pl1_end[0],pl1_end[1],pl2_end[1])
        <<endl;
      else
      {
        cout<<"idBegin="<<idBegin<<" and idEnd="<<idEnd<<endl;
        bool isAccepted=false;
        for(size_t it=idBegin; it<=idEnd && isAccepted==false;it++)
          if(vecMT3D.at(it)==vecThickness.at(noiseLevel))
            isAccepted=true;
        if(isAccepted)
        {
          ATC_Pl1.push_back((*it1_bis));
          ATC_Pl2.push_back((*it2_bis));
          vecBegin.push_back(idBegin);
          vecEnd.push_back(idEnd);
        }
      }
    }
    cout<<"ATC_Pl1 size="<<ATC_Pl1.size()<<endl;
    cout<<"ATC_Pl2 size="<<ATC_Pl2.size()<<endl;
    vecATC_Pl1tmp.push_back(ATC_Pl1);
    vecATC_Pl2tmp.push_back(ATC_Pl2);
    vecBeginTmp.push_back(vecBegin);
    vecEndTmp.push_back(vecEnd);
    noiseLevel++;
  }
  cout<<"vecATC_Pl1tmp size="<<vecATC_Pl1tmp.size()<<endl;
  cout<<"vecATC_Pl1tmp size="<<vecATC_Pl2tmp.size()<<endl;
  /* Select seg for ATC */
  
  /* Remove inclusion */
  //vector<vector<size_t>> vecBegin,vecEnd;
  for(size_t thickness=0; thickness<vecBeginTmp.size()-1; thickness++)
  {
    cout<<"=====> noiseLevel="<<vecThickness.at(thickness)<<endl;
    vector<AlphaThickSegmentComputer2D> ATC_Pl1,ATC_Pl2;
    //vector<size_t> begin,end;
    size_t idBeginCur,idEndCur,idBegin,idEnd;
    for(size_t it=0; it<vecBeginTmp.at(thickness).size();it++)
    {
      idBeginCur=vecBeginTmp.at(thickness).at(it);
      idEndCur=vecEndTmp.at(thickness).at(it);
      bool isInclused=false;
      for(size_t it_bis=0; it_bis<vecBeginTmp.at(thickness+1).size() && isInclused==false;it_bis++)
      {
        idBegin=vecBeginTmp.at(thickness+1).at(it_bis);
        idEnd=vecEndTmp.at(thickness+1).at(it_bis);
        if(idBeginCur>=idBegin && idEndCur<=idEnd)
        {
          cout<<"thickness="<<vecThickness.at(thickness)<<" has seg"<<it<<" to be removed"<<endl;
          isInclused=true;
        }
      }
      if(isInclused==false)
      {
        cout<<"thickness="<<vecThickness.at(thickness)<<" has seg"<<it<<" to be added"<<endl;
        ATC_Pl1.push_back(vecATC_Pl1tmp.at(thickness).at(it));
        ATC_Pl2.push_back(vecATC_Pl2tmp.at(thickness).at(it));
        //begin.push_back(idBeginCur);
        //end.push_back(idEndCur);
      }
    }
    cout<<"ATC_Pl1 size="<<ATC_Pl1.size()<<endl;
    cout<<"ATC_Pl2 size="<<ATC_Pl2.size()<<endl;
    vecATC_Pl1.push_back(ATC_Pl1);
    vecATC_Pl2.push_back(ATC_Pl2);
    //vecBegin.push_back(begin);
    //vecEnd.push_back(end);
  }
  vecATC_Pl1.push_back(vecATC_Pl1tmp.back());
  vecATC_Pl2.push_back(vecATC_Pl2tmp.back());
  cout<<"ATC_Pl1 size="<<vecATC_Pl1.back().size()<<endl;
  cout<<"ATC_Pl2 size="<<vecATC_Pl2.back().size()<<endl;
  /* Remove inclusion */
}

void buildAdaptiveTangentCover(const vector<Z3i::Point>& aContour,
                               const vector<double>& vecThickness,
                               const vector<double>& vecMT3D,
                               const vector<Z2i::Point>& aContour_Pl1,
                               const vector<Z2i::Point>& aContour_Pl2,
                               vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl1,
                               vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl2,
                               const bool validOxy, const bool validOxz, const bool validOyz)
{
  /* TC at different noise levels */
  vector<vector<AlphaThickSegmentComputer2D>> vecTC_Pl1,vecTC_Pl2;
  getTangentCoverDifferentNoise(vecThickness,aContour_Pl1,aContour_Pl2,vecTC_Pl1,vecTC_Pl2);
  /* TC at different noise levels */
  
  /* get segments for ATC */
  getSegmentsAdaptiveTangentCover(aContour,vecThickness,vecMT3D,vecTC_Pl1,vecTC_Pl2,vecATC_Pl1,vecATC_Pl2,validOxy,validOxz,validOyz);
  /* get segments for ATC */
}

void drawAdaptiveTangentCover(MyViewer& display,
                              const vector<Z3i::Point>& aContour,
                              const vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl1,
                              const vector<vector<AlphaThickSegmentComputer2D> >& vecATC_Pl2,
                              const bool validOxy, const bool validOxz, const bool validOyz)
{
  Z3i::RealPoint direction3D,intercept3D,thickness3D;
  Z2i::RealPoint bb11,bb12,bb13,bb14;
  Z2i::RealPoint bb21,bb22,bb23,bb24;
  vector<vector<AlphaThickSegmentComputer2D> >::const_iterator it1=vecATC_Pl1.begin(), it2=vecATC_Pl2.begin();
  size_t noiseLevel=1;
  for(;it1!=vecATC_Pl1.end() && it2!=vecATC_Pl2.end();it1++,it2++)
  {
    vector<AlphaThickSegmentComputer2D>::const_iterator it1_bis=(*it1).begin(),it2_bis=(*it2).begin();
    for(;it1_bis!=(*it1).end() && it2_bis!=(*it2).end();it1_bis++,it2_bis++)
    {
      Z2i::Point pl1_begin=*((*it1_bis).begin());
      Z2i::Point pl2_begin=*((*it2_bis).begin());
      Z2i::Point pl1_end=*((*it1_bis).end()-1);
      Z2i::Point pl2_end=*((*it2_bis).end()-1);
      
      size_t idBegin, idEnd;
      if(validOxy && validOxz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl1_begin[1],pl2_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl1_end[1],pl2_end[1]));
        idEnd=itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      else if(validOxy && validOyz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl1_begin[1],pl2_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl1_end[1],pl2_end[1]));
        idEnd=itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      else//if(validOxz && validOyz)
      {
        vector<Z3i::Point>::const_iterator itBegin = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_begin[0],pl2_begin[0],pl1_begin[1]));
        idBegin=itBegin == aContour.end() ? -1 : itBegin-aContour.begin();
        vector<Z3i::Point>::const_iterator itEnd = std::find(aContour.begin(),aContour.end(),Z3i::Point(pl1_end[0],pl2_end[0],pl1_end[1]));
        idEnd=itEnd == aContour.end() ? -1 : itEnd-aContour.begin();
      }
      
      getParameters((*it1_bis),(*it2_bis),validOxy,validOxz,validOyz,direction3D,intercept3D,thickness3D);
      //cout<<"direction="<<direction3D<<",intercept="<<intercept3D<<",thickness="<<thickness3D<<endl;
      display.setLineColor( getColor(noiseLevel+3));
      drawAsBoundingBox(display,aContour,idBegin,idEnd,(*it1_bis),(*it2_bis),validOxy,validOxz,validOyz);
      
      (*it1_bis).getBoundingBox(bb11,bb12,bb13,bb14);
      (*it2_bis).getBoundingBox(bb21,bb22,bb23,bb24);
      display << CustomColors3D(getColor(noiseLevel+3),  getColor(noiseLevel+3));
      if(validOxy && validOxz) {
        //Oxy
        display.addLine(Z3i::RealPoint(bb11[0],bb11[1],0),Z3i::RealPoint(bb12[0],bb12[1],0));
        display.addLine(Z3i::RealPoint(bb12[0],bb12[1],0),Z3i::RealPoint(bb13[0],bb13[1],0));
        display.addLine(Z3i::RealPoint(bb13[0],bb13[1],0),Z3i::RealPoint(bb14[0],bb14[1],0));
        display.addLine(Z3i::RealPoint(bb14[0],bb14[1],0),Z3i::RealPoint(bb11[0],bb11[1],0));
        //Oxz
        display.addLine(Z3i::RealPoint(bb21[0],0,bb21[1]),Z3i::RealPoint(bb22[0],0,bb22[1]));
        display.addLine(Z3i::RealPoint(bb22[0],0,bb22[1]),Z3i::RealPoint(bb23[0],0,bb23[1]));
        display.addLine(Z3i::RealPoint(bb23[0],0,bb23[1]),Z3i::RealPoint(bb24[0],0,bb24[1]));
        display.addLine(Z3i::RealPoint(bb24[0],0,bb24[1]),Z3i::RealPoint(bb21[0],0,bb21[1]));
      }
      else if(validOxy && validOyz) {
        //Oxy
        display.addLine(Z3i::RealPoint(bb11[0],bb11[1],0),Z3i::RealPoint(bb12[0],bb12[1],0));
        display.addLine(Z3i::RealPoint(bb12[0],bb12[1],0),Z3i::RealPoint(bb13[0],bb13[1],0));
        display.addLine(Z3i::RealPoint(bb13[0],bb13[1],0),Z3i::RealPoint(bb14[0],bb14[1],0));
        display.addLine(Z3i::RealPoint(bb14[0],bb14[1],0),Z3i::RealPoint(bb11[0],bb11[1],0));
        //Oyz
        display.addLine(Z3i::RealPoint(0,bb21[0],bb21[1]),Z3i::RealPoint(0,bb22[0],bb22[1]));
        display.addLine(Z3i::RealPoint(0,bb22[0],bb22[1]),Z3i::RealPoint(0,bb23[0],bb23[1]));
        display.addLine(Z3i::RealPoint(0,bb23[0],bb23[1]),Z3i::RealPoint(0,bb24[0],bb24[1]));
        display.addLine(Z3i::RealPoint(0,bb24[0],bb24[1]),Z3i::RealPoint(0,bb21[0],bb21[1]));
      }
      else {//if(validOxz && validOyz)
        //Oxz
        display.addLine(Z3i::RealPoint(bb11[0],0,bb11[1]),Z3i::RealPoint(bb12[0],0,bb12[1]));
        display.addLine(Z3i::RealPoint(bb12[0],0,bb12[1]),Z3i::RealPoint(bb13[0],0,bb13[1]));
        display.addLine(Z3i::RealPoint(bb13[0],0,bb13[1]),Z3i::RealPoint(bb14[0],0,bb14[1]));
        display.addLine(Z3i::RealPoint(bb14[0],0,bb14[1]),Z3i::RealPoint(bb11[0],0,bb11[1]));
        //Oyz
        display.addLine(Z3i::RealPoint(0,bb21[0],bb21[1]),Z3i::RealPoint(0,bb22[0],bb22[1]));
        display.addLine(Z3i::RealPoint(0,bb22[0],bb22[1]),Z3i::RealPoint(0,bb23[0],bb23[1]));
        display.addLine(Z3i::RealPoint(0,bb23[0],bb23[1]),Z3i::RealPoint(0,bb24[0],bb24[1]));
        display.addLine(Z3i::RealPoint(0,bb24[0],bb24[1]),Z3i::RealPoint(0,bb21[0],bb21[1]));
      }
    }
    noiseLevel++;
  }
}

void drawBoundingBox(MyViewer& display,
                     Z3i::RealPoint pf, Z3i::RealPoint pl,
                     Z3i::RealPoint pf2, Z3i::RealPoint pl2,
                     Z3i::RealPoint pf3, Z3i::RealPoint pl3,
                     Z3i::RealPoint pf4, Z3i::RealPoint pl4)
{
  display.setLineColor( Color::Purple );
  display.addLine(pf, pl);
  //display.setLineColor( Color::Green );
  display.addLine(pf2, pl2);
  //display.setLineColor( Color::Blue );
  display.addLine(pf3, pl3);
  //display.setLineColor( Color::Magenta );
  display.addLine(pf4, pl4);
  
  display.setLineColor( Color::Purple );
  //two end facets
  display.addLine(pf, pf2);
  display.addLine(pf2,pf4);
  display.addLine(pf4, pf3);
  display.addLine(pf3, pf);
  
  display.addLine(pl, pl2);
  display.addLine(pl2, pl4);
  display.addLine(pl4, pl3);
  display.addLine(pl3, pl);
}
