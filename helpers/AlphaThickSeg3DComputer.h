
#if defined(AlphaThickSeg3DComputer_RECURSES)
#error Recursive header files inclusion detected in AlphaThickSeg3DComputer.h
#else // defined(AlphaThickSeg3DComputer_RECURSES)
/** Prevents recursive inclusion of headers. */
#define AlphaThickSeg3DComputer_RECURSES

#if !defined AlphaThickSeg3DComputer_h
/** Prevents repeated inclusion of headers. */
#define AlphaThickSeg3DComputer_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <list>
#include <utility>
#include <array>
#include "DGtal/base/Exceptions.h"
#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/kernel/CInteger.h"
#include "DGtal/geometry/curves/AlphaThickSegmentComputer.h"
#include "DGtal/base/ConstIteratorAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

template <typename TIterator, typename TInteger>
class AlphaThickSeg3DComputer
{
    BOOST_CONCEPT_ASSERT(( concepts::CInteger<TInteger> ) );
    // ----------------------- Types ------------------------------
  public:
     /// Type of integer, devoted to remainders (and intercepts)
    typedef TInteger Integer;

    /// Type which represent quotient of two integers first/second.
    /// typedef std::pair <Integer, Integer> Quotient;
    /// Type of 3d rational point
    /// typedef std::array< Quotient, 3 > PointR3d;

    /// Type of iterator, at least readable and forward
    typedef TIterator ConstIterator;
    ///Self type
    typedef AlphaThickSeg3DComputer< ConstIterator, TInteger > Self;
    ///Reverse type
    typedef AlphaThickSeg3DComputer< ReverseIterator< ConstIterator >,TInteger > Reverse;
    /// Type of 3d digital point
    typedef typename IteratorCirculatorTraits< ConstIterator >::Value Point3d;
    /// Type of 3d digital vector
    typedef typename IteratorCirculatorTraits< ConstIterator >::Value Vector3d;
    /// Type of 3d digital point coordinate
    typedef typename Point3d::Coordinate Coordinate;
    /// Type of 2d digital point
    typedef DGtal::PointVector< 2, Coordinate > Point2d;
    /// Type of 3d double point
    typedef DGtal::PointVector<3,double> PointD3d;
    /// Adapter for iterators
    typedef functors::Projector< SpaceND< 2, Coordinate > > Projector2d;
    /// Iterator over adapter
    typedef ConstIteratorAdapter< ConstIterator, Projector2d, Point2d > IteratorAdapter;
    /// 2D alpha thick segment recognition algorithm
    typedef DGtal::AlphaThickSegmentComputer< Point2d, IteratorAdapter  > AlphaThickSegmentComputer2d;

    // ----------------------- Standard services ------------------------------
  public:


    /**
     * Default constructor.
     * @param maximalThickness maximal thickness of ATS
     */
    AlphaThickSeg3DComputer(const double maximalThickness = 1.0, const int nbPlane = 3);

    /**
     * Constructor with initialisation
     * @param it an iterator
     * @param maximalThickness maximal thickness of ATS
     * @see init
     */
    AlphaThickSeg3DComputer ( const ConstIterator& it, const double maximalThickness = 1.0, const int nbPlane = 3);

    /**
     * Initialisation.
     * @param it an iterator
     */
    void init ( const ConstIterator& it );


    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    AlphaThickSeg3DComputer ( const AlphaThickSeg3DComputer & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    AlphaThickSeg3DComputer & operator= ( const AlphaThickSeg3DComputer & other );

    /**
     * @return a default-constructed instance of Self.
     */
    Self getSelf() const;

    /**
     * @return a default-constructed instance of Reverse.
     */
    Reverse getReverse() const;

    /**
     * Checks whether a point belongs to the ATS or not
     * @param aPoint the point to be checked
     * @return 'true' if yes, 'false' otherwise
     */
    bool isInATS ( const Point3d& aPoint ) const;

    /**
     * Checks whether a point belongs to the ATS or not
     * @param it an iterator on the point to be checked
     * @return 'true' if yes, 'false' otherwise
     */
    bool isInATS ( const ConstIterator & it ) const;

    /**
     * Equality operator.
     * @param other the object to compare with.
     * @return 'true' either if the leaning points perfectly match
     * or if the first leaning points match to the last ones
     * (same ATS scanned in the reverse way)
     * and 'false' otherwise
     */
    bool operator== ( const AlphaThickSeg3DComputer & other ) const;

    /**
     * Difference operator.
     * @param other the object to compare with.
     * @return 'false' if equal
     * 'true' otherwise
     */
    bool operator!= ( const AlphaThickSeg3DComputer & other ) const;

    /**
     * Destructor.
     */
    ~AlphaThickSeg3DComputer(){}

    // ----------------------- Interface --------------------------------------
  public:


    /**
     * Tests whether the current ATS can be extended at the front.
     * Computes the parameters of the extended ATS if yes
     * and adds the point to the current ATS in this case.
     * @return 'true' if yes, 'false' otherwise.
     */
    bool extendFront();


    /**
     * Tests whether the 3d ATS can be extended at the front.
     *
     * @return 'true' if yes, 'false' otherwise
     */
    bool isExtendableFront();

    // ------------------------- Accessors ------------------------------

    /**
     * Computes the parameters
     * (direction, intercept, thickness)
     * of the ATS
     * @param direction direction vector calculated from 2D valid ATS.
     * @param intercept intercept calculated from mu-parameters of 2D valid ATS.
     * @param thickness thickness calculated from omega-parameters of 2D valid ATS.
     */
    void getParameters ( PointD3d& direction, PointD3d& intercept, PointD3d& thickness ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    /**
     *
     * @return begin iterator of the 3d ATS range.
     */
    ConstIterator begin() const;
    /**
     * @return end iterator of the 3d ATS range.
     */
    ConstIterator end() const;

    /**
       @return a const-reference on the arithmetical ATS recognition
       algorithm along the XY plane.
    */
    const AlphaThickSegmentComputer2d & alphaThickSeg2dXY () const;

    /**
       @return a const-reference on the arithmetical ATS recognition
       algorithm along the XZ plane.
    */
    const AlphaThickSegmentComputer2d & alphaThickSeg2dXZ () const;

    /**
       @return a const-reference on the arithmetical ATS recognition
       algorithm along the YZ plane.
    */
    const AlphaThickSegmentComputer2d & alphaThickSeg2dYZ () const;

    /**
       @param i the axis orthogonal to the plane
       i = 0 -> YZ-plane
       i = 1 -> XZ-plane
       i = 2 -> XY-plane
       @return a const-reference on the arithmetical ATS recognition
       algorithm along the plane orthogonal to the \a i-th axis.
    */
    const AlphaThickSegmentComputer2d & alphaThickSeg2d( Dimension i ) const;

    /**
       @param i the axis orthogonal to the plane
       i = 0 -> YZ-plane
       i = 1 -> XZ-plane
       i = 2 -> XY-plane
       @return true if given 2D ATS along orthogonal axis is valid
     */
    bool validalphaThickSeg2d ( Dimension i ) const;


    // ------------------ Display ------------------------------------------
  public:

    /**
     * @return the style name used for drawing this object.
     */
    std::string className() const;

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    // ------------------------- Hidden services ------------------------------
  private:
    /**
     * Tests whether the current 2d-ATS can be extended at the front.
     * Computes the parameters of the extended 2d-ATS if yes
     * and adds the point to the current 2d-ATS in this case.
     * Used internally to simplify extendFront().
     * @param ATS2D reference to 2d-ATSComputer
     * @param blocked reference to status of ATS2D
     * updated if ATS2D cannot be extended at the front.
     * @return 'true' if yes, 'false' otherwise.
     */
    bool extendFront ( AlphaThickSegmentComputer2d & ATS2D, bool & blocked );

    // ------------------------- Protected Datas ------------------------------
  protected:
    /// The maximal thickness of the segment.
    double myMaximalThickness;
    /// Projector for XY-plane.
    Projector2d myProjXY; 
    /// Projector for XZ-plane.
    Projector2d myProjXZ;
    /// Projector for YZ-plane.
    Projector2d myProjYZ;
    /// Number of valide projection for recognizing the segment
    int myNbPlane;
    /// 2d alpha thick segment recognition algorithms for XY-plane.
    AlphaThickSegmentComputer2d myXYalgo;
    double myLengthXY = 0;
    /// 2d alpha thick segment recognition algorithms for XZ-plane.
    AlphaThickSegmentComputer2d myXZalgo;
    double myLengthXZ = 0;
    /// 2d alpha thick segment recognition algorithms for YZ-plane.
    AlphaThickSegmentComputer2d myYZalgo;
    double myLengthYZ = 0;
    //Real length of the segment
    double myLength = 0;
    /**
     * Used internally to store information which 2d alpha thick segment
     * should not be any more extended. This happened when two successive 3D points
     * have same projections onto respective 2d plane.
     */
    bool blockXY, blockXZ, blockYZ;

    /// begin and end iterators
    ConstIterator myBegin, myEnd;
}; // end of class AlphaThickSeg3DComputer

  /**
   * Overloads 'operator<<' for displaying objects of class 'AlphaThickSeg3DComputer'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'AlphaThickSeg3DComputer' to write.
   * @return the output stream after the writing.
   */
  template <typename TIterator, typename TInteger>
  std::ostream&
  operator<< ( std::ostream & out,  const AlphaThickSeg3DComputer<TIterator,TInteger> & object )
  {
    object.selfDisplay( out);
    return out;
  }
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/geometry/curves/AlphaThickSeg3DComputer.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined AlphaThickSeg3DComputer_h

#undef AlphaThickSeg3DComputer_RECURSES
#endif // else defined(AlphaThickSeg3DComputer_RECURSES)
