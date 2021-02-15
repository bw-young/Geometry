/////////////////////////////////////////////////////////////////////
// Objects and methods for managing points of intersection between //
// geometric objects (2d).                                         //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 02/05/2021 - Brennan Young                                      //
// - created                                                       //
// - migrated IPt, isSame, IPtToPoint2d, and reverseForB from      //
//   PolyIntersect.h.                                              //
// 02/08/2021 - Brennan Young                                      //
// - reverseForB renamed to reverseIpt.                            //
// - corrected errors in isCoincident and isSame.                  //
// 02/11/2021 - Brennan Young                                      //
// - overloaded reverseIPt to reverse a single IPt object.         //
// - numerous improvements to the isSame algorithm.                //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_INTERSECTIONPOINT_20210205
#define YOUNG_GEOMETRY_INTERSECTIONPOINT_20210205

#include <vector>   // std::vector
#include "Polyline.hpp" // Point2d, Line2d, Polygon2d, Polyline2d


// Object for tracking an intersection point between two polygons.
class IPt {
public:
    // state constants
    static const char UNSET    = 0;
    static const char UNKNOWN  = 0;
    static const char ENTER    = 1;
    static const char EXIT     = 2;
    static const char EDGE     = 3;
    
    size_t i;       // vertex index in polygon A
    double a;       // fraction of distance from i to i+1
    size_t j;       // vertex index in polygon B
    double b;       // fraction of distance from j to j+1
    char   state;   // indicator of what the point represents
    
    // constructors / destructor
    IPt (size_t I, double A, size_t J, double B, char S=0, char r=0)
        : i(I), a(A), j(J), b(B), state(S) {}
    IPt (const IPt& P)
        : i(P.i), a(P.a), j(P.j), b(P.b), state(P.state) {}
    ~IPt() {}
    
    // operators
    IPt& operator= (const IPt& P)
    {
        if ( this == &P ) return *this;
        i = P.i;
        a = P.a;
        j = P.j;
        b = P.b;
        state = P.state;
        return *this;
    }
    bool operator< (const IPt& P) const
    {
        // sort in order of i, a, j, b (state not considered)
        if ( i < P.i ) return true;
        if ( i > P.i ) return false;
        if ( a < P.a ) return true;
        if ( a > P.a ) return false;
        if ( j < P.j ) return true;
        if ( j > P.j ) return false;
        return b < P.b;
    }
}; // IPt


/////////////////////////////////////////////////////////////////////
// Conversion                                                      //
/////////////////////////////////////////////////////////////////////


// Convert from IPt to Point2d.
// -- Arguments --
// P : IPt to convert.
// A : Container of vertex P.i; vertex must be a Point2d object.
template <class PointArray>
Point2d IPtToPoint2d ( const IPt& P, const PointArray& A )
{
    Line2d L (A[P.i], A[P.i+1]);
    return Point2d (L.a.x + L.r.x * L.len * P.a,
                    L.a.y + L.r.y * L.len * P.a);
}

// Reverse an intersection point -- (i,a) and (j,b) are swapped.
// Does not change the intersection point's state, as this cannot be
// determined with the information given.
IPt reverseIPt ( const IPt& p )
{
    return IPt(p.j, p.b, p.i, p.a, p.state);
}

// Reverse a vector of intersections -- that is, (i,a) and (j,b) are
// swapped, such that the vector of points conforms to the vertex-
// order of B.
std::vector<IPt> reverseIPt ( const std::vector<IPt>& P )
{
    // put points into sorted structure
    std::set<IPt> S;
    for ( size_t i = 0; i < P.size(); ++i )
        S.insert(reverseIPt(P[i]));
    
    // add points to output vector
    std::vector<IPt> out;
    out.insert(out.begin(), S.begin(), S.end());
    
    return out;
}


/////////////////////////////////////////////////////////////////////
// Comparison                                                      //
/////////////////////////////////////////////////////////////////////


// Get previous and next line segments from an intersection point
// involving vertex i in A.
// -- Arguments --
// i    : intersection vertex.
// a    : distance from i toward i+1, where 0 is i, 1 is i+1.
// A    : polygon or polyline object containing i vertex.
// i0   : begining of chain containing i.
// i1   : end of chain containing i (not inclusive).
// loop : true if the chain that contains i loops, false otherwise.
// L0   : output line object for (i-1,i) if a=0, else (i,i+1).
// L1   : output line object for (i,i+1) if a=0, else (i+1,i+2).
// -- Returns --
// 0 if successful, 1 if not (i.e., is crossing a non-looping chain).
// Assigns L0 and L1. If crossing the end of a non-looping chain, L0
// and L1 are both set to (i,i+1).
template <class Poly>
int segmentsBeforeAfter ( size_t i, double a,
    const Poly& A, size_t i0, size_t i1, bool loop,
    Line2d* L0, Line2d* L1 )
{
    bool flag = true;
    
    if ( a < 0.0000001 ) {
        *L1 = Line2d(A[i], A[i+1]);
        if ( i == i0 ) {
            if ( !loop ) {
                *L0 = *L1;
                flag = false;
            }
            else *L0 = Line2d(A[i1-2], A[i]);
        }
        else *L0 = Line2d(A[i-1], A[i]);
    }
    else {
        *L0 = Line2d(A[i], A[i+1]);
        if ( i+2 >= i1 ) {
            if ( !loop ) {
                *L1 = *L0;
                flag = false;
            }
            else *L1 = Line2d(A[i+1],A[i0+1]);
        }
        else *L1 = Line2d(A[i+1],A[i+2]);
    }
    
    return flag ? 0 : 1;
}

// Determine if two intersection points are effectively identical.
// This does not compare their state or the relative orientation of
// the intersecting segments.
// -- Arguments --
// a  : IPt object to compare.
// b  : IPt object to compare.
// i0 : first index of chain that a.i is in.
// i1 : ending index (not inclusive) of chain that a.i is in.
// iL : true if the chain that a.i is in loops.
// j0 : first index of chain that a.j is in.
// j1 : ending index (not inclusive) of chain that a.j is in.
// jL : true if the chain that a.j is in loops.
// -- Returns --
// True if they represent the same point in space, false otherwise.
bool isCoincident ( const IPt& a, const IPt& b,
    size_t i0, size_t i1, bool iL,
    size_t j0, size_t j1, bool jL )
{
    bool aai0 = a.a       < 0.0000001; // in a, hits j at i
    bool aai1 = 1.0 - a.a < 0.0000001; // in a, hits j at i+1
    bool bai0 = b.a       < 0.0000001; // in b, hits j at i
    bool bai1 = 1.0 - b.a < 0.0000001; // in b, hits j at i+1
    bool aaj0 = a.b       < 0.0000001; // in a, hits i at j
    bool aaj1 = 1.0 - a.b < 0.0000001; // in a, hits i at j+1
    bool baj0 = b.b       < 0.0000001; // in b, hits i at j
    bool baj1 = 1.0 - b.b < 0.0000001; // in b, hits i at j+1
    
    bool iSameChain = a.i >= i0 && a.i < i1 && b.i >= i0 && b.i < i1;
    bool jSameChain = a.j >= j0 && a.j < j1 && b.j >= j0 && b.j < j1;
    bool iLoop = iL
        && ((aai1 && a.i+2 == i1 && bai0 && b.i == i0)
        ||  (bai1 && b.i+2 == i1 && aai0 && a.i == i0));
    bool jLoop = jL
        && ((aaj1 && a.j+2 == j1 && baj0 && b.j == j0)
        ||  (baj1 && b.j+2 == j1 && aaj0 && a.j == j0));
    
    int iDiff = ((int) b.i) - ((int) a.i);
    int jDiff = ((int) b.j) - ((int) a.j);
    double aDiff = b.a - a.a;
    double bDiff = b.b - a.b;
    
    // determine if i is the same
    bool iSame =
        (iSameChain
            && ((iDiff == 0 && fabs(aDiff) < 0.0000001)
                || (iDiff == 1 && fabs(aDiff - -1.0) < 0.0000001)
                || (iDiff == -1 && fabs(aDiff - 1.0) < 0.0000001)))
        || (iLoop
            && ((iDiff > 0 && fabs(aDiff - 1.0) < 0.0000001)
                || (iDiff < 0 && fabs(aDiff - -1.0) < 0.0000001)));
    
    // determine if j is the same
    bool jSame =
        (jSameChain
            && ((jDiff == 0 && fabs(bDiff) < 0.0000001)
                || (jDiff == 1 && fabs(bDiff - -1.0) < 0.0000001)
                || (jDiff == -1 && fabs(bDiff - 1.0) < 0.0000001)))
        || (jLoop
            && ((jDiff > 0 && fabs(bDiff - 1.0) < 0.0000001)
                || (jDiff < 0 && fabs(bDiff - -1.0) < 0.0000001)));
    
    return iSame && jSame;
}

// Determine if an intersection represents an approach toward the
// right or left of a line segment. If the intersection lies at a
// vertex (a "joint") in B and the joint turns right (i.e., the
// segment in B after the vertex points right of the first), "right"
// is any approach direction toward the space CW of both segments
// and "left" is any other approach; vice-versa if the joint turns
// left; and joints that don't turn are treated like one segment.
// -- Arguments --
// ip : intersection point being evaluated.
// A  : polygon or polyline object containing i.
// B  : polygon or polyline object containing j.
// -- Returns --
// 1 if (i,i+1) approaches toward the right,
// -1 if toward the left,
// 0 if cannot be determined (e.g., completely parallel).
template <class PolyA, class PolyB>
int jointRight ( const IPt& ip, const PolyA& A, const PolyB& B )
{
    Line2d Li (A[ip.i], A[ip.i+1]);
    
    // middle of segment in B?
    if ( ip.b > 0.0 && ip.b < 1.0 ) {
        Line2d Lj (B[ip.j], B[ip.j+1]);
        return isRight(Li.r, Lj.r);
    }
    
    // intersecting at vertex in B
    else {
        // get position in B
        size_t jj = B.findChain(ip.j);
        size_t j1 = B.chainEnd(jj);
        
        // get line segments in B on either side of intersection
        Line2d Lj0, Lj1;
        int flag = segmentsBeforeAfter(
            ip.j, ip.b, B, B.chain(jj), j1, A.loop(jj), &Lj0, &Lj1);
        
        int R0 = isRight(Li.r, Lj0.r);
        if ( flag > 0 ) return R0;
        
        int R1 = isRight(Li.r, Lj1.r);
        int Rj = isRight(Lj1.r, Lj0.r);
        
        // right turn in B
        if ( Rj > 0 ) {
            if ( (R0 > 0 && R1 > 0)
                    || (R0 == 0 && R1 > 0)
                    || (R1 == 0 && R0 > 0) )
                return 1;
            else return -1;
        }
        
        // left turn in B
        else if ( Rj < 0 ) {
            if ( (R0 < 0 && R1 < 0)
                    || (R0 == 0 && R1 < 0)
                    || (R1 == 0 && R0 < 0) )
                return -1;
            else return 1;
        }
        
        // no turn in B
        else return R0;
    }
    
    return 0;
}

// Determine if two intersection points represent exactly the same
// intersection (they don't if spatially coincident but one
// is entering and the other is exiting).
// -- Arguments --
// a : IPt object to compare.
// b : IPt object to compare.
// A : Polygon or Polyline object containing i vertices.
// B : Polygon object containing j vertices.
// -- Returns --
//  0 = they are not the same.
//  1 = they are the same and one can be discarded.
//  2 = they are the same and not necessary for characterizing enter/
//      exit into a polygon (they mirror vertices/edges in B).
template <class Poly>
int isSame ( const IPt& a, const IPt& b, const Poly& A,
    const Polygon2d& B )
{
    static const double pi = 4.0 * atan(1.0);
    
    // chain limits in A
    size_t ic = A.findChain(a.i);
    size_t i0 = A.chain(ic);
    size_t i1 = ic+1 < A.nchains() ? A.chain(ic+1) : A.size();
    bool iL = A.loop(ic);
    
    // chain limits in B
    size_t jc = B.findChain(a.j);
    size_t j0 = B.chain(jc);
    size_t j1 = jc+1 < B.nchains() ? B.chain(jc+1) : B.size();
    bool jL = B.loop(jc);
    
    // if they are not spatially coincident or in the same chain,
    // no intersection removal is necessary
    if ( !isCoincident(a, b, i0, i1, iL, j0, j1, jL) ) return 0;
    
    // if both intersections represent the same distance from i
    // to i+1, can remove one of them UNLESS intersecting at a
    // vertex in B, where B touches but does not penetrate A
    if ( a.i == b.i ) {
        if ( (a.a > 0 && a.a < 1) && !(a.b > 0 && a.b < 1) ) {
            Line2d ki (A[a.i], A[a.i+1]); // (i,i+1) in a
            Line2d kj (B[a.j], B[a.j+1]); // (j,j+1) in a
            Line2d mj (B[b.j], B[b.j+1]); // (j,j+1) in b
            
            int aR = isRight(kj.r, ki.r);
            int bR = isRight(mj.r, ki.r);
            if ( aR != 0 && bR != 0 && aR != bR ) return 0;
        }
        
        return 1;
    }
    
    Line2d ki (A[a.i], A[a.i+1]); // (i,i+1) in a
    Line2d kj (B[a.j], B[a.j+1]); // (j,j+1) in a
    Line2d mi (A[b.i], A[b.i+1]); // (i,i+1) in b
    Line2d mj (B[b.j], B[b.j+1]); // (j,j+1) in b
    
    int kkR = isRight(ki.r, kj.r);
    int kmR = isRight(ki.r, mj.r);
    int mkR = isRight(mi.r, kj.r);
    int mmR = isRight(mi.r, mj.r);
    
    // determine if same direction
    bool kkD = angle2d(ki, kj) < 0.5 * pi;
    bool kmD = angle2d(ki, mj) < 0.5 * pi;
    bool mkD = angle2d(mi, kj) < 0.5 * pi;
    bool mmD = angle2d(mi, mj) < 0.5 * pi;
    
    // if parallel-to-parallel, may discard both -- don't need to
    // track chain mirroring
    if ( (kkR == 0 && mmR == 0
                && ((kkD && mmD) || (!kkD && !mmD)))
            || (kmR == 0 && mkR == 0
                && ((kmD && mkD) || (!kmD && !mkD))) )
        return 2;
    
    // discard neither if both entering and exiting at intersection,
    // otherwise discard one
    int aR = jointRight(a, A, B);
    int bR = jointRight(b, A, B);
    if ( aR != 0 && bR != 0 && aR != bR ) return 0;
    
    // otherwise discard one intersection, as it only needs to be
    // represented with one
    return 1;
}


#endif // YOUNG_GEOMETRY_INTERSECTIONPOINT_20210205