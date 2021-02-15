# Geometry
Objects and operations for representing geometric objects and problems.

This is an ongoing work in progress and is liable to change at any time without warning. For example, a planned changes includes generalization of the objects to n dimensions, which will deprecate the present Point2d, Point3d, etc., etc. framework currently in place. Another change that may occur is the generalization of vertex chains for both polylines and polygons, perhaps even to representing polylines and polygons as the same data structure.

Likewise, this documentation will expand and become more complex as needed.

## Installation
The simplest way to include these objects and functions is:
```cpp
#include "Geometry.hpp"
```

A proper DLL or SDK is low-priority for now, but may come in future updates.

## Contents
| Objects |
| --- |
| [Extent](#extent) |
| [Point2d](#point2d) |
| [Point3d](#point2d) |
| [Line2d](#line2d) |
| [Polyline2d](#polyline2d) |
| [Polygon2d](#polygon2d) |
| [IPt](#ipt) |

| Operations |
| --- |
| [Conversion](#conversion-operations) |
| [Distance](#scale-operations) |
| [Topological](#topological-operations) |
| [Boolean](#boolean-operations) |

## Objects

### Extent
Header file: Extent.hpp

Purpose: Represent the bounding box of an n-dimensional (hyper-)volume.

```cpp
class Extent {
  size_t n;                      // number of dimensions
  double *a;                     // minimum bounds
  double *b;                     // maximum bounds
};

Extent e (N);                    // initialize to N dimensions, minimum N=1; default N=1
Extent e (x0,x1,b);              // initialize to 1-dimensional extent from min(x0,x1)-b to max(x0,x1)+b; default b=0
Extent e (x0,y0,x1,y1,b);        // initialize to 2-dimensonal extent; default b=0
Extent e (x0,y0,z0,x1,y1,z1,b);  // initialize to 3-dimensional extent; default b=0
Extent e (e2);                   // copy constructor

e2 = e;                          // assignment copies the object's contents
e < e2;                          // checks dimensions, then minimum bound, then maximum bound
e > e2;                          // checks dimensions, then minimum bound, then maximum bound
e == e2;                         // true if number of dimensions and all values are exactly the same
e != e2;                         // false if number of dimensions and all values are exactly the same

e.size();                        // get number of dimensions
e.min(i); e.max(i);              // access minimum/maximum bound in dimension i
e.xmin(); e.xmax();              // access minimum/maximum bound in dimension 0
e.ymin(); e.ymax();              // access minimum/maximum bound in dimension 1
e.zmin(); e.zmax();              // access minimum/maximum bound in dimension 2
e.len(i);                        // get the length of dimension i: b[i] - a[i]
e.xlen();                        // get the length of dimension 0: b[0] - a[0]
e.ylen();                        // get the length of dimension 1: b[1] - a[1]
e.zlen();                        // get the length of dimension 2: b[2] - a[2]
e.volume();                      // get the (hyper-)volume of the bounding box

e.add(e2); e3 = e.add(e2);       // enlarge the bounding box to include another bounding box
e.clip(e2);                      // reduce to the volume occupied by both extents
e3 = e.intersect(e2);            // get the overlap extent
e + e2; e += e2;                 // enlarge the bounding box to include another bounding box
e3 = e - e2; e -= e2;            // get the overlap extent
e.overlap(e2);                   // true if extents overlap (checks for intersection volume > 0)
```

### Point2d
Header file: Point.hpp

Purpose: Represent a 2d vertex.

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

Purpose: Represent a 3d vertex.

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

### Line2d
Header file: Line.hpp

Purpose: Represent a 2d line segment.

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

### Polyline2d
Header file: Polyline.hpp

Purpose: Represent 2d vertex chains.

```cpp
class Polyline2d{
  struct Chain {size_t start, bool loop}; // chain break structure
  std::vector<Point2d> v; // vertices
  std::vector<Chain> s;   // chain breaks
  Extent ext;             // minimum bounding box
}; // Polyline2d

Polyline2d L (V,c);       // construct with a vector of Point2d objects that is looping if c=true, and determine the bounds; default is an empty vector
Polyline2d L (P);         // construct from a Polygon2d object
Polyline2d L (L);         // copy constructor
L2 = L;                   // assignment copies the object's contents
L[i];                     // access element i in v, vector of Point2d objects
L < L2;                   // checks extent, then number of chains, then number of vertices, then each vertex
L == L2;                  // true if chains and vertices are the same
L != L2;                  // false if chains and vertices are the same
L + L2; L += L2;          // append polyline
L.size();                 // get number of vertices
L.nchains();              // get number of chains
L.findChain(j);           // get the index for the chain that contains vertex j
L.chain(i);               // get the index for the first vertex in chain i
L.chainEnd(i);            // get 1 + the index of the last vertex in chain i
L.loop(i);                // true if chain i loops
L.extent();               // get ext
L.vertices();             // access the vertex vector

L.addChainBreak(j,c);     // add a chain break at vertex j that loops if c=true; if a chain already begins at j, updates instead
L.removeChainBreak(i);    // remove chain break i; if either side of the break is looping, the new chain will only be looping if the last and first vertices are the same
L.addChain(i,L,c);        // adds a chain in the ith chain position from a vector of Point2d objects that is looping if c=true
L.extractChain(i);        // get a Polyline2d object representing chain i
L.removeChain(i);         // remove the chain break and all vertices associated with chain i
L.append_front(i,L,c);    // append the vector of Point2d objects, which loops if c=true, to the beginning of chain i
L.append_back(i,L,c);     // append the vector of Point2d objects, which loops if c=true, to the end of chain i
L.insert(j,p);            // insert point p before vertex j
L.insert(j,beg,end);      // insert point objects as vertices before vertex j, from iterator beg to but not including iterator end
L.push_back(p);           // append point p as a vertex to the end of the vertex chain
L.erase(j);               // delete vertex j
L.computeBounds();        // determine minimum and maximum coordinates, and pushes 0 to L.s if L.s is empty but L.v is not.
```

### Polygon2d
Header file: Polygon.hpp

Purpose: Represent a 2d polygon.

```cpp
class Polygon2d{
  std::vector<Point2d> v; // vertices; each polygon begins and ends with equivalent vertices
  std::vector<size_t> s;  // beginning of polygons (vertex chains)
  Extent ext;             // minimum bounding box
}; // Polygon2d

Polygon2d (V);            // construct with a vector of Point2d objects; default is an empty vector
Polygon2d (P);            // copy constructor
P2 = P;                   // assignment copies the object's contents
P[i];                     // access element i in v
P < P2;                   // checks bounding coordinates, then number of vertices, then each vertex
P == P2;                  // true if vertices are the same
P != P2;                  // false if vertices are the same
P.size();                 // get number of vertices
P.nchains();              // number of polygons (vertex chains)
P.findChain(j);           // get the index for the chain that contains vertex j
L.chain(i);               // get the index for the first vertex in chain i
L.chainEnd(i);            // get 1 + the index of the last vertex in chain i
L.loop(i);                // true if chain i loops (always true for polygons)
L.extent();               // get ext
L.vertices();             // access the vertex vector

L.addChain(i,L);          // adds a chain in the ith chain position from a vector of Point2d objects
L.extractChain(i);        // get a Polygon2d object representing chain i
L.removeChain(i);         // remove the chain and all vertices associated with chain i
L.insert(j,p);            // insert point p before vertex j
L.insert(j,beg,end);      // insert point objects as vertices before vertex j, from iterator beg to but not including iterator end
L.push_back(p);           // append point p as a vertex to the end of the vertex chain; if the last vertex is the same as the one beginning the chain, creates a new chain beginning with this vertex
L.erase(j);               // delete vertex j
P.computeBounds();        // determine ext
```

### IPt
Header file: PolyIntersect.hpp

Purpose: Represent the intersection point between two objects, where that intersection can be represented with a vertex and a distance to the next vertex.

```cpp
class IPt{
  static const char UNSET, UNKNOWN, ENTER, EXIT, EDGE; // public state constants
  
  size_t i;     // vertex index in polygon A
  double a;     // fraction of distance from i to i+1
  size_t j;     // vertex index in polygon B
  double b;     // fraction of distance from j to j+1
  char state;   // indicator of what the point represents
}; // IPt

IPt::UNSET;     // get state constant
IPt::UNKNOWN;   // get state constant
IPt::ENTER;     // get state constant
IPt::EXIT;      // get state constant
IPt::EDGE;      // get state constant
IPt(I,A,J,B,S); // construct an intersection point with vertex-distance pairs (I,A) and (J,B) with state S, default=IPt(0,0,0,0,IPt::UNSET)
IPt(ipt);     // copy constructor
ipt2 = ipt;   // assignment copies the object's contents
ipt < ipt2;   // checks i, a, j, then b

// Return true if two intersection point objects are spatially coincident and within the same chain; chain extents are [i0,i1) (looping if iL=true) for i vertices, and [j0,j1) for j vertices (looping if jL=true)
bool isCoincident (IPt a, IPt b, size_t i0, size_t i1, bool iL, size_t j0, size_t j1, bool jL);

// Returns 1 if line (ipt.i,ipt.i+1) in A points to the right of (ipt.j) in B, -1 if left, or 0 if neither; "right" is any approach direction toward the space CW of both segments if interescting B at a vertex where the segments on either side of the vertex turns right; vice-versa for left
template <class PolyA, class PolyB> // PolyA and PolyB are either polyline or polygon objects
int jointRight (IPt ipt, PolyA A, Poly B);

// Return 0 if two intersection points are not the same, 1 if they are the same and one can be discarded, and 2 if they are the same and neither is necessary for characterizing entry/exit into a polygon; determine if intersection points ipt and ipt2 represent exactly the same intersection of A (contains i vertices) and B (contains j vertices); the points aren't the same if spatially coincident but one is entering and the other is exiting
template <class Poly> // Poly is either a polyline of polygon object
int isSame (ipt, ipt2, Poly A, Polygon2d B)
```

## Conversion Operations
```cpp
// swap (i,a) and (j,b)
IPt reverseIPt (IPt ipt);

// swap (i,a) and (j,b) and reorder
std::vector<IPt> reverseIPt (std::vector<IPt> ipt);

// convert from IPt to Point2d
template <class PointArray> // should be either a polyline or polygon object, or any other container of Point2d objects accessible with the brackets [] operator
Point2d IPt2Point2d (IPt ipt, PointArray A);

// get the previous and next line segments (L0 and L1) from an intersection point (i,a) involving vertex i in A confined within vertex chain [i0,i1); returns 0 if successful, 1 if not (i.e., crossing the end of a non-looping vertex chain)
template<Poly> // Poly should be either a polyline or a polygon object
int segmentsBeforeAfter (size_t i, double a, Poly A, size_t i0, size_t i1, bool loop, Line2d* L0, Line2d* L1)
```

## Scale Operations
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

// maximum square euclidean distance between lines ab and cd
template <class Line> double sqMaxDistToLine2d (Line ab, Line cd);
double sqMaxDistToLine (Line2d ab, Line2d cd);

// maximum euclidean distance between lines ab and cd
template <class Line> double maxDistToLine2d (Line ab, Line cd);
double maxDistToLine(Line2d ab, Line2d cd);

// minimum square euclidean distance between the given point and the nearest edge of the polygon
double sqDistToEdge2d (Polygon2d P, Point2d p);

// minimum euclidean distance between the given point and the nearest edge of the polygon
double distToEdge2d (Polygon2d P, Point2d p);

// minimum square euclidean distance between the given line segment and the nearest edge of the polygon
double sqDistToEdge2d (Polygon2d P, Line2d s);

// minimum euclidean distance between the given line segment and the nearest edge of the polygon
double distToEdge2d (Polygon2d P, Line2d s);

// angle between lines ab and cd, treated as vectors
double angle2d (Line2d ab, Line2d cd);
```

## Topological Operations
```cpp
// get the winding number for p in P. p is in P if w > 0 && w % 2 != 0, where w is the winding number
template <class Point> int inPoly (Polygon2d P, Point p);

// get a vector of distances from r0 to each point where the ray r0->r1 intersects the polygon; note that argument rmax presently does nothing
template <class Point> std::vector<double> rayPoly_intersect2d (Point r0, Point r1, Polygon2d P, double rmax);

// get a vector of distances from s.a to each point where the line segment s intersects the polygon
template <class Point> std::vector<double> segmentPoly_intersect2d (Polygon2d P, Line2d s);

// Return true if the chain in A that contains ipt[k].i should be treated as beginning inside of B; false otherwise
template <class Poly> // Poly should be either a polyline or polygon object
bool validateInOut (std::vector<IPt> p, size_t k, Poly A, Polygon2d B)

// Return a vector of IPt objects with assigned states and culled of all points
// that do not represent a new entrance or exit into/from B
template <class Poly> // Poly should be either a polyline or polygon object
std::vector<IPt> intersectionPoints2dState (std::vector<IPt> ipt, Poly A, Polygon2s B)

// Get all points where a polygon or polyline intersects a polygon (e.g., for use in boolean operations)
template <class Poly>  // Poly should be either a polyline or polygon object
std::vector<IPt> intersectionPoints2d (Poly A, Polygon2d B)
```

## Boolean Operations
```cpp
// AND = intersection = clip
bool intersect (Point2d p, Line2d ab, double t);  // true if p lies within tolerance distance t of ab
bool intersect (Line2d ab, Point2d p, double t);

bool intersect (Point2d p, Polygon2d P);          // true if p is inside P with tolerance distance t (calls and interprets inPoly)
bool intersect (Polygon2d P, Point2d p);

Polygon2d polyAND (Polygon2d A, Polygon2d B);
Polygon2d intersect (Polygon2d A, Polygon2d B);   // calls polyAND
Polygon2d clip2d (Polygon2d A, Polygon2d B);      // calls PolyAND
Polyline2d intersect (Polyline2d L, Polygon2d P); // get polyline clipped by polygon

// OR = union = merge = dissolve
// <!> UNDER CONSTRUCTION <!>

// NOT
// <!> UNDER CONSTRUCTION <!>

// XOR
// <!> UNDER CONSTRUCTION <!>
```
