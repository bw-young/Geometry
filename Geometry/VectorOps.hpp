/////////////////////////////////////////////////////////////////////
// Methods for vector operations using Geometry objects.           //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 02/04/2021 - Brennan Young                                      //
// - created                                                       //
// 02/09/2021 - Brennan Young                                      //
// - migrate angle(Line) to Line.hpp.                              //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_GEOMETRY_VECTOROPS_20210204
#define YOUNG_GEOMETRY_VECTOROPS_20210204

#include "Point.hpp" // Point2d, Point3d


// Get the dot product of two vectors, represented by Point objects.
// -- Arguments --
// a : Point object representing a vector.
// b : Point object representing a vector.
// -- Returns --
// A floating point scalar. Divide by the magnitude of b to get the
// length of a projected onto b.
double dot ( const Point2d& a, const Point2d& b )
{ return a.x * b.x + a.y * b.y; }
double dot ( const Point3d& a, const Point3d& b )
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

// Determine if vector a points to the right (clockwise) of b.
// -- Arguments --
// a : Point object representing a vector.
// b : Point object representing a vector.
// -- Returns --
// 1 if a is clockwise of b, -1 if counterclockwise, 0 if parallel
// (or in opposite direction).
int isCW ( const Point2d& a, const Point2d& b )
{
    double d = dot(a, Point2d(b.y, -b.x));
    if ( d > 0.0 ) return 1;
    if ( d < 0.0 ) return -1;
    return 0;
}
int isRight ( const Point2d& a, const Point2d& b )
{ return isCW(a, b); }


#endif // YOUNG_GEOMETRY_VECTOROPS_20210204