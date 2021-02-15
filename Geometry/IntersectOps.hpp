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
// - added IPtToPoint2d.                                           //
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
// 01/21/2021 - Brennan Young                                      //
// - corrected error where intersectionPoints2d(Polgon,Polygon)    //
//   would not detect sub-polygon terminations in polygon B.       //
// 01/25/2021 - Brennan Young                                      //
// - made numerous adjustments to mitigate error in intersection(  //
//   Polygon, Polyline).                                           //
// 01/26/2021 - Brennan Young                                      //
// - fixed errors in intersection(Polygon, Polyline).              //
// 02/01/2021 - Brennan Young                                      //
// - added 'state' member to IPt to track whether the intersection //
//   represents entering, exiting, or tracing a polygon.           //
// 02/04/2021 - Brennan Young                                      //
// - added intersectionPoints2dState.                              //
// - intersectionPoints2d(Polygon, Polyline) now takes the line    //
//   first.                                                        //
// - intersectionPoints2d(polyline, poly) now relies on            //
//   intersectionPoints2dState to resolve duplicate/unnecessary    //
//   intersection points.                                          //
// 02/05/2021 - Brennan Young                                      //
// - added isSame.                                                 //
// - migrated IPt, isSame(), IPtToPoint2d(), and reverseForB() to  //
//   IPt.hpp.                                                      //
// - added intersection(Point2d, Line2d).                          //
// - renamed from PolyIntersect.hpp to IntersectOps.hpp.           //
// - update to utilize Polygon's new structure.                    //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_POLYINTERSECT_20201228
#define YOUNG_GEOMETRY_POLYINTERSECT_20201228

#include <cstddef>            // size_t
#include <set>                // std::set
#include <vector>             // std::vector

#include "Polyline.hpp"       // Point2d, Line2d, Polygon2d,
                              //   Polyline2d, inPoly()
#include "VectorOps.hpp"      // dot(), isRight(), angle()
#include "IPt.hpp"            // IPt, IPtToPoint2d(), reverseForB()
#include "../BYstdlib/sort.h" // quickSort()


/////////////////////////////////////////////////////////////////////
// Intersect - Point-Line                                          //
/////////////////////////////////////////////////////////////////////


// Return if a point intersects a line segment.
// -- Arguments --
// p : point to evaluate.
// L : line to evaluate.
// t : tolerance distance (within and including this distance).
// -- Returns --
// True if p is within t of L, false otherwise.
bool intersect ( const Point2d& p, const Line2d& L, double t )
{ return sqDistToLine2d(p, L) <= t*t; }
bool intersect ( const Line2d& L, const Point2d& p, double t )
{ return intersect(p, L, t); }


/////////////////////////////////////////////////////////////////////
// Intersect - Point-Polygon                                       //
/////////////////////////////////////////////////////////////////////


// Return if a point intersects a polygon.
// -- Arguments --
// p : point to evaluate.
// G : polygon to evaluate.
// -- Returns --
// True if p is within G, false otherwise.
bool intersect ( const Point2d& p, const Polygon2d& G )
{
    int w = inPoly(G, p); // winding number
    return w > 0 && w % 2 != 0;
}
bool intersect ( const Polygon2d& G, const Point2d& p )
{ return intersect(p, G); }


/////////////////////////////////////////////////////////////////////
// Intersect - Polygon-Polygon and Polyline-Polygon                //
/////////////////////////////////////////////////////////////////////


// Validate whether the chain in A begins inside B, accounting for
// the ambiguous case where the chain begins on the edge of B.
// -- Arguments --
// p : vector of intersection points. Must be sorted by i, a, j, b.
// k : intersection point to check (first in intersection in chain).
// A : polygon or polyline object containing all i in p.
// B : polygon object being intersected.
// -- Returns --
// True if the chain containing p[k].i should be treated as beginning
// inside of B; false otherwise.
template <class Poly>
bool validateInOut ( const std::vector<IPt>& p, size_t k,
    const Poly& A, const Polygon2d& B )
{
    size_t ii = A.findChain(p[k].i); // chain containing p[k].i
    size_t i0 = A.chain(ii);         // beginning of chain
    size_t i1 = A.chainEnd(ii);      // end of chain
    
    // get first determinable right- or left-crossing intersection
    int R = 0;
    for ( size_t kk = k; R == 0 && kk < p.size() && p[kk].i < i1;
            ++kk )
        R = jointRight(p[kk], A, B);
    if ( R == 0 ) return intersect(A[i0], B); // cannot be determined
    
    // get segment perpendicular to (i0,i0+1) and passing through i0
    Line2d L (A[i0], A[i0+1]);
    // a is to the left of (i0,i0+1)
    Point2d a (A[i0].x - 0.001 * L.r.y,
               A[i0].y + 0.001 * L.r.x);
    // b is to the right of (i0,i0+1)
    Point2d b (A[i0].x + 0.001 * L.r.y,
               A[i0].y - 0.001 * L.r.x);
    L = Line2d(a, b);
    std::vector<double> d = segmentPoly_intersect2d(B, L);
    if ( d.size() == 0 ) return intersect(A[i0], B); // not at edge
    
    // if leave B to the right, check if inside to the left, or
    // vice-versa
    if ( R > 0 ) return intersect(a, B);
    else return intersect(b, B);
}

// Determine the state of all given intersection points and removes
// points unless they represent an entering/exiting of the polygon.
// -- Arguments --
// ip : vector of intersection points, where (i,a) refers to A and
//      (j,b) refers to B.
// A  : an organized container of vertices (Polyline or Polygon).
// B  : the polygon being intersected with.
// -- Returns --
// A vector of intersection points.
template <class Poly>
std::vector<IPt> intersectionPoints2dState (
    const std::vector<IPt>& ip, const Poly& A,
    const Polygon2d& B )
{
    // prepare output
    std::vector<IPt> p = ip;
    
    // compare each point with every other point and remove
    // all points that don't result in a change of enter/exit
    // state
    for ( size_t k = 0; k+1 < p.size(); ++k ) {
        bool flag = false;
        
        for ( size_t m = k+1; m < p.size(); ++m ) {
            // determine if they represent the same point
            // and whether to discard
            int same = isSame(p[k], p[m], A, B);
            
            // discard one
            if ( same == 1 ) {
                p.erase(p.begin() + m);
                --m;
                continue;
            }
            
            // discard both
            else if ( same == 2 ) {
                p.erase(p.begin() + m);
                p.erase(p.begin() + k);
                --k;
                flag = true;
                break;
            }
        }
        
        if ( flag ) continue;
    }
    
    // discard points if no state-change between the last
    // intersection point and this (check for situations where A
    // enters the edge of B at one intersection point and then leaves
    // the edge at a further intersection point, but does not change
    // whether A is in B)
    for ( size_t k = 0; k+1 < p.size(); ++k ) {
        size_t ii0 = A.findChain(p[k].i);
        size_t ii1 = A.findChain(p[k+1].i);
        size_t jj0 = B.findChain(p[k].j);
        size_t jj1 = B.findChain(p[k+1].j);
        int R0 = jointRight(p[k],   A, B); // k on right?
        int R1 = jointRight(p[k+1], A, B); // k+1 on right?
        bool same = isCoincident(p[k], p[k+1],
            A.chain(ii0), A.chainEnd(ii0), A.loop(ii0),
            B.chain(jj0), B.chainEnd(jj0), B.loop(jj0));
        
        // in same chains but not spatially coincident
        if ( ii0 == ii1 && jj0 == jj1 && !same
                && (R0 == 0 || R1 == 0 || R0 == R1) ) {
            p.erase(p.begin() + k+1);
            --k;
            continue;
        }
    }
    
    // determine the state change associated with each intersection
    // (entering or exiting)
    bool AinB = false;
    bool newChain = true;
    size_t ii = 0; // current chain in A
    size_t i1 = A.chainEnd(ii);
    for ( size_t k = 0; k < ip.size(); ++k ) {
        // get position in A
        size_t i = ip[k].i;
        while ( i >= i1 && ii+1 < A.nchains() ) {
            newChain = true;
            ++ii;
            i1 = A.chainEnd(ii);
        }
        
        // determine if chain in A begins inside B
        if ( newChain ) {
            AinB = validateInOut(p, k, A, B);
            
            // mark whether the intersection represents entering B
            p[k].state = AinB ? IPt::EXIT : IPt::ENTER;
            AinB = !AinB;
        }
        
        // continuing in chain -- toggle in/out
        else {
            p[k].state = AinB ? IPt::EXIT : IPt::ENTER;
            AinB = !AinB;
        }
        
        newChain = false;
    }
    
    p.shrink_to_fit();
    return p;
}

// Get all points where a polygon or polyline intersects a polygon.
// e.g., for use in clipping, union, or xor operations.
// -- Arguments --
// A : Polygon or Polyline A.
// B : Polygon B.
// -- Returns --
// A vector of IPt objects.
template <class Poly>
std::vector<IPt> intersectionPoints2d ( const Poly& A,
    const Polygon2d& B )
{
    std::vector<IPt> ip;
    ip.reserve(A.size());
    
    // for each chain in A
    for ( size_t k = 0; k < A.nchains(); ++k ) {
        // for each vertex in chain in A
        size_t i1 = (k+1 < A.nchains() ? A.chain(k+1) : A.size());
        for ( size_t i = A.chain(k); i+1 < i1; ++i ) {
            
            Line2d AL (A[i], A[i+1]); // (i,i+1)
            
            // skip if (i,i+1) extent doesn't overlap with B
            Extent Aext (AL.a.x, AL.a.y, AL.b.x, AL.b.y, 1.0);
            if ( !Aext.overlap(B.extent()) ) continue;
            
            // for each chain in B
            for ( size_t m = 0; m < B.nchains(); ++m ) {
                // for each vertex in chain in B
                size_t j1 =
                    (m+1 < B.nchains() ? B.chain(m+1) : B.size());
                for ( size_t j = B.chain(m); j+1 < j1; ++j ) {
                    Line2d BL (B[j], B[j+1]); // (j,j+1)
                    
                    // distance to intersection from i toward i+1
                    double Ad = lineIntersect2d(AL, BL);
                    if ( Ad < 0.0 ) continue; // no intersection
                    double Bd = lineIntersect2d(BL, AL);
                    
                    // relationship between segments
                    // 1=(i,i+1) points to right of (j,j+1)
                    // -1=left
                    // 0=parallel
                    int R = isRight(AL.r, BL.r);
                    
                    // add this intersection point
                    ip.push_back(
                        IPt(i, Ad / AL.len, j, Bd / BL.len, 0, R));
                }
            }
        }
    }
    
    ip.shrink_to_fit(); // free excess memory
    bystd::quickSort(&ip); // sort IPt objects
    
    // characterize points' state and clear unnecessary points
    ip = intersectionPoints2dState(ip, A, B);
    return ip;
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
    std::vector<IPt> ip = intersectionPoints2d(A, B);
    
    // remove points that are simultaneously entry and exit points
    // (vertices that touch the edge of the other polygon)
    for ( size_t i = 0; i+1 < ip.size(); ++i ) {
        double d = sqDist(IPtToPoint2d(ip[i], A),
            IPtToPoint2d(ip[i+1], A));
        if ( d < 0.0000001
                && ((ip[i].state == IPt::ENTER
                    && ip[i+1].state == IPt::EXIT)
                || (ip[i].state == IPt::EXIT
                    && ip[i+1].state == IPt::ENTER)) ) {
            ip.erase(ip.begin()+i);
            ip.erase(ip.begin()+i);
            --i;
        }
    }
    
    // if no intersection points, return inner polygon if overlaps
    if ( ip.size() == 0 ) {
        if ( A.size() > 0 && inPoly(B, A[0]) ) return A;
        if ( B.size() > 0 && inPoly(A, B[0]) ) return B;
        return Polygon2d();
    }
    
    // trace and mark entry and exit points to the other polygon
    std::vector<Point2d> At, Bt; // polygons with intersection points
    std::vector<size_t> An, Bn; // index of neighbor in other polygon
                                // at intersection, otherwise -1
    std::vector<bool> Ae, Be; // true if inside or entering other
                              // polygon
    At.reserve(A.size() + ip.size());
    An.reserve(A.size() + ip.size());
    Ae.reserve(A.size() + ip.size());
    Bt.reserve(B.size() + ip.size());
    Bn.reserve(B.size() + ip.size());
    Be.reserve(B.size() + ip.size());
    
    int w = inPoly(B, A[0]);
    bool state = w > 0 && w % 2 != 0; // true if in B
    size_t k = 0;
    for ( size_t i = 0; i < A.size(); ++i ) {
        // add if vertex is NOT an intersection point
        if ( !(k < ip.size() && i == ip[k].i
                && ip[k].a < 0.0000001) )
        {
            At.push_back(A[i]);
            An.push_back(-1);
            Ae.push_back(state);
        }
        
        // add intersection points
        for ( ; k < ip.size() && i == ip[k].i;
                ++k )
        {
            state = !state;
            At.push_back(IPtToPoint2d(ip[k], A));
            An.push_back(ip[k].j);
            Ae.push_back(state);
            ip[k].i = At.size() - 1;
        }
    }
    
    ip = reverseIPt(ip);
    w = inPoly(A, B[0]);
    state = w > 0 && w % 2 != 0; // true if in A
    k = 0;
    for ( size_t i = 0; i < B.size(); ++i ) {
        // add if vertex is NOT an intersection point
        if ( !(k < ip.size() && i == ip[k].i
                && ip[k].a < 0.0000001) )
        {
            Bt.push_back(B[i]);
            Bn.push_back(-1);
            Be.push_back(state);
        }
        
        // add intersection points
        for ( ; k < ip.size() && i == ip[k].i;
                ++k )
        {
            state = !state;
            An[ip[k].j] = Bt.size();
            Bt.push_back(IPtToPoint2d(ip[k], B));
            Bn.push_back(ip[k].j);
            Be.push_back(state);
        }
    }
    ip = std::vector<IPt> ();
    
    // create the output polygon
    Polygon2d out;
    std::vector<Point2d> v;
    v.reserve(At.size());
    
    bool flag1, flagB;
    size_t K = -1; // beginning of current polygon
    k = 0; // current position
    for ( size_t i = 0; i < At.size(); ++i ) {
        // skip to next unprocessed intersection
        if ( !(An[i] < Bt.size()) ) continue;
        
        // create a new polygon
        v.push_back(At[i]);
        out.addChain(out.size(), v);
        v.clear();
        
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
                    v.push_back(flagB ? Bt[k] : At[k]);
                }
            }
            // if current is exiting other polygon
            else {
                // add vertices until next intersection
                while ( (!flagB && !(An[k] < Bt.size()))
                        || (flagB && !(Bn[k] < At.size())) ) {
                    --k;
                    v.push_back(flagB ? Bt[k] : At[k]);
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

// Return the polyline where a polyline and a polygon overlap.
// -- Arguments --
// A : polyline to clip.
// B : clipping polygon.
// -- Returns --
// A Polyline object containing chains of vertices where A intersects
// B.
Polyline2d intersect ( const Polyline2d& A, const Polygon2d& B )
{
    if ( B.size() + A.size() == 0 ) return Polyline2d();
    
    // get intersection points
    std::vector<IPt> ip = intersectionPoints2d(A, B);
    size_t k = 0; // intersection
    
    // remove points that are simultaneously entry and exit points
    // (vertices that touch the edge of the other polygon)
    for ( size_t i = 0; i+1 < ip.size(); ++i ) {
        double d = sqDist(IPtToPoint2d(ip[i], A),
            IPtToPoint2d(ip[i+1], A));
        if ( d < 0.0000001
                && ((ip[i].state == IPt::ENTER
                    && ip[i+1].state == IPt::EXIT)
                || (ip[i].state == IPt::EXIT
                    && ip[i+1].state == IPt::ENTER)) ) {
            ip.erase(ip.begin()+i);
            ip.erase(ip.begin()+i);
            --i;
        }
    }
    
    // if no intersection points, return inner polyline if overlaps
    if ( ip.size() == 0 ) {
        if ( A.size() > 0 && inPoly(B, A[0]) ) return A;
        return Polyline2d();
    }
    
    // build list of vertices
    Polyline2d out;
    std::vector<Point2d> chain;
    chain.reserve(A.size() + ip.size());
    
    // for each chain in A
    for ( size_t ii = 0; ii < A.nchains(); ++ii ) {
        // get end of chain
        size_t i1 = A.chainEnd(ii);
        
        // get first intersection point in chain
        for ( ; k < ip.size() && ip[k].i < A.chain(ii); ++k ) {}
        
        // is first vertex in chain inside polygon?
        bool state = false;
        if ( k < ip.size() ) state = ip[k].state == IPt::EXIT;
        else state = intersect(A[A.chain(ii)], B);
        size_t m = out.nchains(); // 1st output chain from this chain
        
        // for each vertex in chain in A
        for ( size_t i = A.chain(ii); i < i1; ++i ) {
            // add vertex if in polygon and NOT an intersection point
            if ( state
                    && !(k < ip.size()
                    && i == ip[k].i
                    && ip[k].a < 0.0000001)
                    && !(k-1 < ip.size()
                    && i-1 == ip[k-1].i
                    && fabs(ip[k-1].a - 1) < 0.0000001) ) {
                chain.push_back(A[i]);
                
                // check for end of chain
                if ( i+1 == i1 ) {
                    // check if entire loop is not inside polygon
                    if ( m < out.nchains() && A.loop(ii)
                            && chain.back() == out[out.chain(m)] ) {
                        chain.pop_back();
                        out.append_front(m, chain, false);
                    }
                    else out.addChain(out.nchains(), chain, true);
                    
                    // reset chain
                    chain.clear();
                    chain.reserve(
                        A.size() - i + ip.size() - k);
                }
            }
            
            // add intersection points
            for ( ; k < ip.size() && i == ip[k].i; ++k ) {
                state = !state;
                
                // entering polygon
                if ( state && chain.size() > 0 ) {
                    // add currently collected chain to output
                    // (if closing an input loop, attach to the chain
                    // in output that represents the start of the
                    // looping input chain)
                    if ( i+1 == i1 && A.loop(ii)
                            && m < out.nchains() )
                        out.append_front(m, chain, false);
                    else out.addChain(out.nchains(), chain, false);
                    
                    // reset chain
                    chain.clear();
                    chain.reserve(
                        A.size() - i + ip.size() - k);
                }
                
                chain.push_back(IPtToPoint2d(ip[k], A));
            }
        }
    }
    
    if ( chain.size() > 0 )
        out.addChain(out.nchains(), chain, false);
    
    return out;
}


/////////////////////////////////////////////////////////////////////
// OR

//


#endif // YOUNG_GEOMETRY_POLYINTERSECT_20201228