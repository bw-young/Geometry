/////////////////////////////////////////////////////////////////////
// Objects and methods for managing polygon intersections and      //
// boolean operations.                                             //
//                                                                 //
// The concepts implemented here may be broadly attributed to:     //
//                                                                 //
// Greiner, G. and Hormann, K., 1998, Efficient clipping of        //
//   arbitrary polygons. ACM Transactions on Graphics 17(2), 71-   //
//   83.                                                           //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Operations                       Functions                      //
// -----------------------          -------------------------      //
// AND = intersection = clip        polyAND, intersect, clip2d     //
// OR = union = merge = dissolve                                   //
// NOT                                                             //
// XOR                                                             //
//
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/28/2020 - Brennan Young                                      //
// - created                                                       //
// 01/05/2021 - Brennan Young                                      //
// - IPt now also compares using i and j.                          //
// - added IPt2Point2d.                                            //
// - added shrink_to_fit after creating vectors.                   //
// - added reverseForB.                                            //
// - intersection2d renamed to clip2d.                             //
// - clip2d further developed. updating intersections while        //
//   building At causes the algorithm to fail.                     //
// 01/07/2021 - Brennan Young                                      //
// - clip2d completed, renamed polyAND and intersect.              //
// 01/11/2021 - Brennan Young                                      //
// - added intersect for Point (true/false point lies in polygon). //
// - added intersect between Polyline and Polygon.                 //
// - added intersectionPoints2d for Polyline and Polygon           //
//   comparison.                                                   //
// - corrected error in intersectionsPoints2d that resulted in an  //
//   incorrect comparison of distances.                            //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYINTERSECT_20201228
#define YOUNG_GEOMETRY_POLYINTERSECT_20201228

#include <cstddef> // size_t
#include <set>     // std::set
#include <vector>  // std::vector

#include "Polyline.hpp" // Point2d, Line2d
#include "Polygon.hpp" // Point2d, Line2d, Polygon2d, inPoly()
#include "../BYstdlib/sort.h" // quickSort


/////////////////////////////////////////////////////////////////////
// Object: Intersection Point.                                     //
/////////////////////////////////////////////////////////////////////


// Object for tracking an intersection point between two polygons.
class IPt {
public:
    size_t i; // vertex index in polygon A
    double a; // fraction of distance from i to i+1
    size_t j; // vertex index in polygon B
    double b; // fraction of distance from j to j+1
    
    // constructors / destructor
    IPt (size_t I, double A, size_t J, double B)
        : i(I), a(A), j(J), b(B) {}
    IPt (const IPt& P) : i(P.i), a(P.a), j(P.j), b(P.b) {}
    ~IPt() {}
    
    // operators
    IPt& operator= (const IPt& P)
    {
        if ( this == &P ) return *this;
        i = P.i;
        a = P.a;
        j = P.j;
        b = P.b;
        return *this;
    }
    bool operator< (const IPt& P) const
    {
        // sort by i, then by a, then by j, then by b
        if ( i < P.i ) return true;
        if ( i > P.i ) return false;
        if ( a < P.a ) return true;
        if ( a > P.a ) return false;
        if ( j < P.j ) return true;
        if ( j > P.j ) return false;
        return b < P.b;
    }
}; // IPt

// Convert from IPt to Point2d.
// -- Arguments --
// P : IPt to convert.
// A : Container of vertex i (as a Point2d object).
template <class PointArray>
Point2d IPt2Point2d ( const IPt& P, const PointArray& A )
{
    Line2d L (A[P.i], A[P.i+1]);
    return Point2d (L.a.x + L.r.x * L.len * P.a,
                    L.a.y + L.r.y * L.len * P.a);
}

// Reverse a vector of intersections for traversing B -- that is,
// (i,a) and (j,b) are swapped, such that A and B are swapped.
std::vector<IPt> reverseForB ( const std::vector<IPt>& P )
{
    // put points into sorted structure
    std::set<IPt> S;
    for ( size_t i = 0; i < P.size(); ++i )
        S.insert(IPt(P[i].j, P[i].b, P[i].i, P[i].a));
    
    // add points to output vector
    std::vector<IPt> out;
    out.insert(out.begin(), S.begin(), S.end());
    
    return out;
}


/////////////////////////////////////////////////////////////////////
// Operations                                                      //
/////////////////////////////////////////////////////////////////////


// Get all points where polygons A and B intersect. e.g., for use in
// clipping, union, or xor operations. Complexity: O(m*n)
// -- Arguments --
// A : Polygon A.
// B : Polygon B.
// -- Returns --
// A vector of IPt objects.
std::vector<IPt> intersectionPoints2d ( const Polygon2d& A,
    const Polygon2d& B )
{
    std::vector<IPt> intersections;
    intersections.reserve(A.size());
    
    size_t m, j, i0=0, i1=0, j0=0, j1=0;
    Line2d Li, Lj;
    bool dd=false, d0=false, R=false, R0=false, Rp=false;
    double di, dj, dot;
    for ( size_t i = 0; i+1 < A.size(); ++i ) {
        // start of polygon
        if ( i == i1 ) {
            i0 = i1;
            m = -1;
        }
        
        // end of polygon
        else if ( A[i] == A[i0] ) {
            i1 = i+1; // start of next polygon
            continue;
        }
        
        // get line object to represent segment i -> i+1
        Li = Line2d (A[i], A[i+1]);
        
        for ( j = 0; j+1 < B.size(); ++j ) {
            // start of polygon
            if ( j == j1 ) j0 = j1;
            
            // end of polygon
            else if ( B[j] == B[j0] ) {
                j1 = j+1; // start of next polygon
                continue;
            }
            
            // get line object to represent segment j -> j+1
            Lj = Line2d (B[j], B[j+1]);
            
            // get distance to intersection
            di = raySegment_intersect2d(Li.a, Li.r, Lj.a, Lj.b);
            
            // skip if doesn't intersect segment
            if ( di < 0 || di > Li.len ) continue;
            di = di / Li.len;
            
            // determine if same as first or previous intersection
            d0 = A[i+1] == A[i0]
                && m < intersections.size()
                && fabs(intersections[m].a - di) < 0.0000001;
            dd = intersections.size() > 0
                && fabs(intersections.back().a - di) < 0.0000001;
            
            // get dot product of right-orthogonal vector (+ = clockwise)
            dot = Li.r.x * (-Li.dy) + Li.r.y * Li.dx;
            R = dot > 0.0;
            
            if ( i == i0 ) R0 = Rp = R;
            
            // record intersection only if not duplicating a vertex
            // unless ray touches but doesn't cross polygon boundary
            // (that is, both exits and enters the polygon)
            if ( !((d0 && R0 == R) || (dd && Rp == R)) ) {
                dj = raySegment_intersect2d(Lj.a, Lj.r, Li.a, Li.b);
                intersections.push_back(
                    IPt(i, di, j, dj/Lj.len));
                if ( i == i0 ) {
                    m = intersections.size() - 1;
                    R0 = R;
                }
            }
            
            Rp = R;
        }
    }
    
    intersections.shrink_to_fit(); // free excess reserved memory
    bystd::quickSort(&intersections); // sort IPt objects
    
    return intersections;
}

// Get all points where polygon A and polyline B intersect. e.g., for
// use in clipping, union, or xor operations. Complexity: O(m*n)
// -- Arguments --
// A : Polygon A.
// B : Polyline B.
// -- Returns --
// A vector of IPt objects.
std::vector<IPt> intersectionPoints2d ( const Polygon2d& A,
    const Polyline2d& B )
{
    std::vector<IPt> intersections;
    intersections.reserve(A.size());
    
    size_t m, j, i0=0, i1=0, J=1;
    Line2d Li, Lj;
    bool dd=false, d0=false, R=false, R0=false, Rp=false;
    double di, dj, dot;
    for ( size_t i = 0; i+1 < A.size(); ++i ) {
        // start of polygon
        if ( i == i1 ) {
            i0 = i1;
            m = -1;
        }
        
        // end of polygon
        else if ( A[i] == A[i0] ) {
            i1 = i+1; // start of next polygon
            continue;
        }
        
        // get line object to represent segment i -> i+1
        Li = Line2d (A[i], A[i+1]);
        
        for ( j = 0; j+1 < B.size(); ++j ) {
            // end of chain
            if ( J < B.s.size() && j+1 == B.s[J] ) {
                ++J;
                continue;
            }
            
            // get line object to represent segment j -> j+1
            Lj = Line2d (B[j], B[j+1]);
            
            // get distance to intersection
            di = raySegment_intersect2d(Li.a, Li.r, Lj.a, Lj.b);
            
            // skip if doesn't intersect segment
            if ( di < 0 || di > Li.len ) continue;
            di = di / Li.len;
            
            // determine if same as first or previous intersection
            d0 = A[i+1] == A[i0]
                && m < intersections.size()
                && fabs(intersections[m].a - di) < 0.0000001;
            dd = intersections.size() > 0
                && fabs(intersections.back().a - di) < 0.0000001;
            
            // get dot product of right-orthogonal vector (+ = clockwise)
            dot = Li.r.x * (-Li.dy) + Li.r.y * Li.dx;
            R = dot > 0.0;
            
            if ( i == i0 ) R0 = Rp = R;
            
            // record intersection only if not duplicating a vertex
            // unless ray touches but doesn't cross polygon boundary
            // (that is, both exits and enters the polygon)
            if ( !((d0 && R0 == R) || (dd && Rp == R)) ) {
                dj = raySegment_intersect2d(Lj.a, Lj.r, Li.a, Li.b);
                intersections.push_back(
                    IPt(i, di, j, dj/Lj.len));
                if ( i == i0 ) {
                    m = intersections.size() - 1;
                    R0 = R;
                }
            }
            
            Rp = R;
        }
    }
    
    intersections.shrink_to_fit(); // free excess reserved memory
    bystd::quickSort(&intersections); // sort IPt objects
    
    return intersections;
}


/////////////////////////////////////////////////////////////////////
// AND

// Return the polygon where two polygons overlap.
// -- Arguments -- 
// A : polygon to clip.
// B : clipping polygon.
// -- Returns --
// The clipped polygon. If the two objects do not intersect, returns
// a polygon with no vertices.
Polygon2d polyAND ( const Polygon2d& A, const Polygon2d& B )
{
    // get all edge-pair intersections
    std::vector<IPt> intersections = intersectionPoints2d(A, B);
    
    // if no intersection points, return inner polygon if overlaps
    if ( intersections.size() == 0 ) {
        if ( A.size() > 0 && inPoly(B, A[0]) ) return A;
        if ( B.size() > 0 && inPoly(A, B[0]) ) return B;
        return Polygon2d();
    }
    
    // remove points that are simultaneously entry and exit points
    // (vertices that touch the edge of the other polygon)
    for ( size_t i = 0; i+1 < intersections.size(); ++i ) {
        double d = sqDist(IPt2Point2d(intersections[i], A),
            IPt2Point2d(intersections[i+1], A));
        if ( d < 0.0000001 ) {
            intersections.erase(intersections.begin()+i);
            intersections.erase(intersections.begin()+i);
            --i;
        }
    }
    
    // trace and mark entry and exit points to the other polygon
    std::vector<Point2d> At, Bt; // polygons with intersection points
    std::vector<size_t> An, Bn; // index of neighbor in other polygon
                                // at intersection, otherwise -1
    std::vector<bool> Ae, Be; // true if inside or entering other
                              // polygon
    At.reserve(A.size() + intersections.size());
    An.reserve(A.size() + intersections.size());
    Ae.reserve(A.size() + intersections.size());
    Bt.reserve(B.size() + intersections.size());
    Bn.reserve(B.size() + intersections.size());
    Be.reserve(B.size() + intersections.size());
    
    int w = inPoly(B, A[0]);
    bool state = w > 0 && w % 2 != 0; // true if in B
    size_t k = 0;
    for ( size_t i = 0; i < A.size(); ++i ) {
        // add if vertex is NOT an intersection point
        if ( !(k < intersections.size() && i == intersections[k].i
                && intersections[k].a < 0.0000001) )
        {
            At.push_back(A[i]);
            An.push_back(-1);
            Ae.push_back(state);
        }
        
        // add intersection points
        for ( ; k < intersections.size() && i == intersections[k].i;
                ++k )
        {
            state = !state;
            At.push_back(IPt2Point2d(intersections[k], A));
            An.push_back(intersections[k].j);
            Ae.push_back(state);
            intersections[k].i = At.size() - 1;
        }
    }
    
    intersections = reverseForB(intersections);
    w = inPoly(A, B[0]);
    state = w > 0 && w % 2 != 0; // true if in A
    k = 0;
    for ( size_t i = 0; i < B.size(); ++i ) {
        // add if vertex is NOT an intersection point
        if ( !(k < intersections.size() && i == intersections[k].i
                && intersections[k].a < 0.0000001) )
        {
            Bt.push_back(B[i]);
            Bn.push_back(-1);
            Be.push_back(state);
        }
        
        // add intersection points
        for ( ; k < intersections.size() && i == intersections[k].i;
                ++k )
        {
            state = !state;
            An[intersections[k].j] = Bt.size();
            Bt.push_back(IPt2Point2d(intersections[k], B));
            Bn.push_back(intersections[k].j);
            Be.push_back(state);
        }
    }
    intersections = std::vector<IPt> ();
    
    // create the output polygon
    Polygon2d out;
    out.v.reserve(At.size());
    
    bool flag1, flagB;
    size_t K = -1; // beginning of current polygon
    k = 0; // current position
    for ( size_t i = 0; i < At.size(); ++i ) {
        // skip to next unprocessed intersection
        if ( !(An[i] < Bt.size()) ) continue;
        
        // create a new polygon
        out.v.push_back(At[i]);
        K = k = i;
        flag1 = flagB = false; // first pt in polygon; current in At
        
        // traverse polygons until the polygon closes
        while ( (!flagB && An[k] < Bt.size())
                || (flagB && Bn[k] < At.size()) ) {
            
            // remove from list of unprocessed intersections
            if ( !flagB ) An[k] = -1;
            else          Bn[k] = -1;
            
            // check if polygon is closed
            if ( flag1 && !flagB && k == K ) break;
            flag1 = true;
            
            // if current is entering other polygon
            if ( (!flagB && Ae[k]) || (flagB && Be[k]) ) {
                // add verticecs until next intersection
                while ( (!flagB && !(An[k] < Bt.size()))
                        || (flagB && !(Bn[k] < At.size())) ) {
                    ++k;
                    out.v.push_back(flagB ? Bt[k] : At[k]);
                }
            }
            // if current is exiting other polygon
            else {
                // add vertices until next intersection
                while ( (!flagB && !(An[k] < Bt.size()))
                        || (flagB && !(Bn[k] < At.size())) ) {
                    --k;
                    out.v.push_back(flagB ? Bt[k] : At[k]);
                }
            }
            
            // jump to neighbor
            size_t n = k;
            k = flagB ? Bn[k] : An[k];
            flagB = !flagB;
            
            // remove from list of unprocessed intersections
            if ( flagB ) An[n] = -1;
            else         Bn[n] = -1;
        }
    }
    
    out.v.shrink_to_fit(); // free excess reserved memory
    
    return out;
}

Polygon2d intersect ( const Polygon2d& A, const Polygon2d& B )
{
    return polyAND(A, B);
}

Polygon2d clip2d ( const Polygon2d& A, const Polygon2d& B )
{
    return polyAND(A, B);
}

// Return the polyline where a polygon and a polyline overlap.
Polyline2d intersect ( const Polygon2d& poly,
    const Polyline2d& line )
{
    if ( poly.size() + line.size() == 0 ) return Polyline2d();
    
    // get intersection points
    std::vector<IPt> intersections =
        intersectionPoints2d(poly, line);
    intersections = reverseForB(intersections);
    
    // determine if the first vertex in the line is inside the
    // polygon
    int w = inPoly(poly, line[0]);
    bool state = w > 0 && w % 2 != 0; // true if in polygon
    bool newChain = true;
    
    // build list of vertices
    Polyline2d out;
    out.v.reserve(line.size() + poly.size());
    out.s.reserve(line.s.size() + intersections.size());
    
    size_t I = 0; // polyline chain
    size_t k = 0;
    for ( size_t i = 0; i < line.size(); ++i ) {
        if ( I < line.s.size() && i == line.s[I] ) {
            newChain = true;
            ++I;
        }
        
        // add if vertex is NOT an intersection point
        if ( state &&
                !(k < intersections.size() && i == intersections[k].i
                && intersections[k].a < 0.0000001) )
        {
            out.v.push_back(line[i]);
            if ( newChain ) out.s.push_back(out.v.size()-1);
            newChain = false;
        }
        
        // add intersection points
        for ( ; k < intersections.size() && i == intersections[k].i;
                ++k )
        {
            state = !state;
            if ( state ) newChain = true;
            out.v.push_back(IPt2Point2d(intersections[k], line));
            if ( newChain ) out.s.push_back(out.v.size()-1);
            newChain = false;
        }
    }
    
    out.v.shrink_to_fit();
    out.s.shrink_to_fit();
    return out;
}

// Return if a point intersects a polygon.
bool intersect ( const Polygon2d& poly, const Point2d& p )
{
    int w = inPoly(poly, p); // winding number
    return w > 0 && w % 2 != 0;
}


/////////////////////////////////////////////////////////////////////
// OR

//


#endif // YOUNG_GEOMETRY_POLYINTERSECT_20201228