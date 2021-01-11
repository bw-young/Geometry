# Geometry
Objects and operations for representing geometric objects and problems.

This is an ongoing work in progress and is liable to change at any time without warning. For example, a planned changes includes generalization of the objects to n dimensions, which will deprecate the present Point2d, Point3d, etc., etc. framework currently in place. Another change that may occur is the generalization of vertex chains for both polylines and polygons, perhaps even to representing polylines and polygons as the same data structure.

Likewise, this documentation will expand and become more complex as needed.

The simplest way to include these objects and functions is:
```cpp
#include "Geometry.hpp"
```

## Objects

### Point2d
Header file: Point.hpp

```cpp
class Point2d {double x, y;};

Point2d p (x,y);  // construct with (x,y); default (0,0)
Point2d p2 (p);   // copy constructor

p2 = p;           // assignment copies the object's contents
p < p2;           // checks y then x
p > p2;           // checks y then x
p == p2;          // true if x and y are exactly the same.
p != p2;          // false if x and y are exactly the same.
```

### Point3d
Header file: Point.hpp

```cpp
class Point3d {double x, y, z;};

Point3d p (x,y,z);  // construct with (x,y,z); default (0,0,0)
Point3d p2 (p);     // copy constructor

p2 = p;             // assignment copies the object's contents
p < p2;             // checks z, y, then x
p > p2;             // checks z, y, then x
p == p2;            // true if x, y, and z are exactly the same.
p != p2;            // false if x, y, and z are exactly the same.
```

### Line
Header file: Line.hpp

```cpp
class Line2d{
  Point2d a;      // segment start point
  Point2d b;      // segment end point
  Point2d r;      // segment direction vector a->b
  double dx;      // b.x - a.x
  double dy;      // b.y - a.y
  double len;     // segment length
}; // Line2d

Line2d s (p1,p2); // construct with (p1,p2) and computes all other members of s; default (Point2d(),Point2d())
Line2d s2 (s);    // copy constructor
s2 = s;           // assignment copies the object's contents
s < s2;           // true if both points in s are less than either point in s2
s == s2;          // true if points in s are equal to the same points in s2
s != s2;          // false if points in s are equal to the same points in s2
```

### Polyline
Header file: Polyline.hpp

```cpp
class Polyline2d{
  std::vector<Point2d> v; // vertices
  std::vector<size_t> s;  // indexes in v that start a new chain
  Point2d bMin;           // minimum bounding coordinates
  Point2d bMax;           // maximum bounding coordinates
}; // Polyline2d

Polyline2d (V);           // construct with a vector of Point2d objects and determine the bounds; default is an empty vector
Polyline2d (L);           // copy constructor
L2 = L;                   // assignment copies the object's contents
L[i];                     // access element i in v, vector of Point2d objects
L < L2;                   // checks bounding coordinates, then number of chains, then number of vertices, then each vertex
L == L2;                  // true if chains and vertices are the same
L != L2;                  // false if chains and vertices are the same
L.size();                 // get number of vertices
L.chains();               // get number of chains
L.computeBounds();        // determine minimum and maximum coordinates, and pushes 0 to L.s if L.s is empty but L.v is not.
```

### Polygon
Header file: Polygon.hpp

```cpp
class Polygon2d{
  std::vector<Point2d> v; // vertices; each polygon begins and ends with equivalent vertices
  Point2d bMin;           // minimum bounding coordinates
  Point2d bMax;           // maximum bounding coordinates
}; // Polygon2d

Polygon2d (V);            // construct with a vector of Point2d objects; default is an empty vector
Polygon2d (P);            // copy constructor
P2 = P;                   // assignment copies the object's contents
P[i];                     // access element i in v, vector of Point2d objects
P < P2;                   // checks bounding coordinates, then number of vertices, then each vertex
P == P2;                  // true if vertices are the same
P != P2;                  // false if vertices are the same
P.size();                 // get number of vertices
P.computeBounds();        // determine minimum and maximum coordinates
```

### IPt
Header file: PolyIntersect.hpp

```cpp
class IPt{
  size_t i;   // vertex index in polygon A
  double a;   // fraction of distance from i to i+1
  size_t j;   // vertex index in polygon B
  double b;   // fraction of distance from j to j+1
}; // IPt

IPt(I,A,J,B); // construct an intersection point with vertex-distance pairs (I,A) and (J,B)
IPt(ipt);     // copy constructor
ipt2 = ipt;   // assignment copies the object's contents
ipt < ipt2;   // checks i, a, j, then b
```

## Conversion Operations
```cpp
// convert from IPt to Point2d; PointArray could be a std::vector<Point2d>, Polyline2d, or Polygon2d object, or any other object that has Point2d objects accessible with the brackets [] operator.
template <class PointArray> Point2d IPt2Point2d (IPt ipt, PointArray A);

// swap (i,a) and (j,b) and reorder
std::vector<IPt> reverseForB (std::vector<IPt> ipt);
```

## Distance Operations
```cpp
// square euclidean distance between points a and b
template <class Point> double sqDist2d (Point a, Point b);
double sqDist(Point2d a, Point2d b);
  
// euclidean distance between points a and b
template <class Point> double dist2d (Point a, Point b);
double dist(Point2d a, Point2d b);

// euclidean distance to ray r0->r1 intersection with line segment p0->p1, negative otherwise
template <class Point> double raySegment_intersect2d (Point r0, Point r1, Point p0, Point p1);

// square euclidean distance between point a and line cd
template <class Point, class Line> double sqDistToLine2d (Point a, Line cd);
double sqDistToLine (Point2d a, Line2d cd);

// euclidean distance between point a and line cd
template <class Point, class Line> double distToLine2d (Point a, Line cd);
double distToLine (Point2d a, Line2d cd);

// minimum square euclidean distance between lines ab and cd
template <class Line> double sqDistToLine2d (Line ab, Line cd);
double sqDistToLine (Line2d ab, Line2d cd);

// minimum euclidean distance between lines ab and cd
template <class Line> double distToLine2d (Line ab, Line cd);
double distToLine (Line2d ab, Line2d cd);

// minimum square euclidean distance between the given point and the nearest edge of the polygon
double sqDistToEdge2d (Polygon2d P, Point2d p);

// minimum euclidean distance between the given point and the nearest edge of the polygon
double distToEdge2d (Polygon2d P, Point2d p);

// minimum square euclidean distance between the given line segment and the nearest edge of the polygon
double sqDistToEdge2d (Polygon2d P, Line2d s);

// minimum euclidean distance between the given line segment and the nearest edge of the polygon
double distToEdge2d (Polygon2d P, Line2d s);
```

## Topological Operations
```cpp
// get the winding number for p in P. p is in P if w > 0 && w % 2 != 0, where w is the winding number
template <class Point> int inPoly (Polygon2d P, Point p);

// get a vector of distances from r0 to each point where the ray r0->r1 intersects the polygon; note that argument rmax presently does nothing
template <class Point> std::vector<double> rayPoly_intersect2d (Point r0, Point r1, Polygon2d P, double rmax);

// get a vector of distances from s.a to each point where the line segment s intersects the polygon
template <class Point> std::vector<double> segmentPoly_intersect2d (Polygon2d P, Line2d s);

// get all points where the polygons intersects another geometry object (e.g., for use in clipping, union, or xor operations), ordered for polygon A
std::vector<IPt> intersectionPoints2d (Polygon2d A, Polygon2d B);
std::vector<IPt> intersectionPoints2d (Polygon2d A, Polyline2d B)
```

## 2D Boolean Operations
```cpp
// AND = intersection = clip
Polygon2d polyAND (Polygon2d A, Polygon2d B);
Polygon2d intersect (Polygon2d A, Polygon2d B);   // calls polyAND
Polygon2d clip2d (Polygon2d A, Polygon2d B);      // calls PolyAND
Polyline2d intersect (Polygon2d P, Polyline2d L); // get polyline clipped by polygon
bool intersect (Polygon2d P, Point2d p);          // true if p is inside P (calls and interprets inPoly)

// OR = union = merge = dissolve
// <!> UNDER CONSTRUCTION <!>

// NOT
// <!> UNDER CONSTRUCTION <!>

// XOR
// <!> UNDER CONSTRUCTION <!>
```
