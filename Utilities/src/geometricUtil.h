#ifndef GEOMETRICUTIL_H
#define GEOMETRICUTIL_H

#include "metroUtil.h"
#include <math.h>

#define PI 3.14159265

/**
  Structure representing a geometic point
*/
struct Point{
    double x;
    double y;
    double z;
};

/**
  @param x - the x coord of the point.
  @param y - the y coord of the point.
  @param z - the z coord of the point.
  @return - a point at x,y,z
*/
Point createPoint(double X, double Y, double Z);

struct Vector{
    Point origin;
    Point end;
};

Vector createVector(Point startPoint, Point endPoint);

Vector createVector(Point startPoint, double deltaX, double deltaY, double deltaZ);

/**
  Structure representing a geometric plane.
  A plane is described by three points.
*/
struct Plane{
    Atom atom1;
    Atom atom2;
    Atom atom3;
};

void printPoint(Point p);

/**
  @param p1 - the first point in the plane.
  @param p2 - the second point in the plane.
  @param p3 - the third point in the plane.
  @return - a plane represented by points p1, p2 and p3
*/
Plane createPlane(Atom p1, Atom p2, Atom p3);

/**
  Returns a specific atom from a vector of atoms.
  @param atoms - vector of atoms from which to get the Atom
  @param atomID - the ID of the atom to find.
  @return - the atom with atomID or Atom with parameters = -1 if not found.
*/
Atom getAtom(vector<Atom> atoms, unsigned long atomID);

/**
  Returns a bond (if it exists) containing a1 and a2.
  @param bonds - the vector of bonds to be tested.
  @param a1 - the first atom needed in the bond
  @param a2 - the second atom needed in the bond
  @return - the bond containing a1 and a2
          - Bond(-1,-1,false) if a suitable bond has not been found
 */
Bond getBond(vector<Bond> bonds, unsigned long a1, unsigned long a2);

/**
  Returns a vector of all of the atoms bonded to atom id
  @param bonds - the vector of bonds to be checked.
  @param atomID - the atom that must be involved in the bonds
  @return - a vector of atomIDs that are bonded to  atomID
*/
vector<unsigned long> getAllBonds(vector<Bond> bonds, unsigned long atomID);

/**
  Returns the intersection of the two vectors.
  @param v1 - the first vector  set.
  @param v2 - the seconds vector set.
  @return - vector containing the elements of the intersection of the two sets.
*/
vector<unsigned long> getIntersection(vector<unsigned long> v1, vector<unsigned long> v2);

/**
  Determines wheter the atom is a member of the vector.
  @param atoms - the list of atom ids to check.
  @param toCheck - the atomID to be checked for in the vector.
  @return - true if toCheck is a member of the vector.
            false if toCheck is not a member of the vector.
*/
bool isMember(vector<unsigned long> atoms, unsigned long toCheck);

/**
  @param degrees - the measure of an angle in degrees
  @return - the value of the same angle in radians.
*/
double degreesToRadians(double degrees);

/**
  @param radians - the measure of an angle in radians
  @return - the value of the same angle in degrees
*/
double radiansToDegrees(double radians);

/**
  Returns the id in the bond that is not atomID
  @param bond - the bond to be checked.
  @param atomID - the atomID that you want the opposite of.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the bond
*/
unsigned long getOppositeAtom(Bond bond, unsigned long atomID);

/**
  Returns the id in the angle that is not atomID
  @param angle - the angle to be checked.
  @param atomID - the atomID that you want the opposite of.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the angle
*/
unsigned long getOppositeAtom(Angle angle, unsigned long atomID);

/**
  Returns the id of the atom in the dihedral that is not atomID
  @param dihedral - the dihedral to be checked.
  @param atomID - the atomID that do not want to return.
  @return - the atom id of the atom opposite atomID
            -1 if atomID is not found in the dihedral 
*/
unsigned long getOppositeAtom(Dihedral dihedral, unsigned long atomID);

/**
  Returns the id of the third atom in the angle between atom1 and atom2.
  Returns -1 if no common angle is found.
  @param bonds - vector of bonds to be checked.
  @param atom1 - the first atom in the angle.
  @param atom2 - the second atom in the angle.
  @return - the id of the third atom in the angle or -1 if none exists.
*/
unsigned long getCommonAtom(vector<Bond> bonds, unsigned long atom1,
        unsigned long atom2);

/**
  @param - the first atom.
  @param - the second atom.
  @param - the distance between atom1 and atom2.
*/
double getDistance(Atom atom1, Atom atom2);

/**
  Determines if an atom is on the xz plane
  @atom - the atom to examine.
  @return - true if the atom is in the xz plane
*/
bool inXZPlane(Atom atom);

/**
  Returns the angle made by the three atoms.
  @param atom1 - edge of the angle. Bonded to atom2
  @param atom2 - the corner of the angle, Bonded to atom1 and atom3.
  @param atom3 - the other edge of the angle.  Bonded to atom2.
  @return - the angle in degrees created by the three atoms.
*/
double getAngle(Atom atom1, Atom atom2, Atom atom3);

/**
  Returns the acute angle between the two planes.
  @param p1 - the first plane that creates the angle.
  @param p2 - the second plane that creates the angle.
  @return - the acute angle in degrees between the planes.
          - returns 90 degrees if the planes are perdenicular.
*/
double getAngle(Plane p1, Plane p2);

/**
  Returns the normal vector of a plane.
  The vector has a start point of 0,0,0 and an end point of @return.
  @param p - the plane used to find the normal.
  @return - the end point of the normal vector
*/
Point getNormal(Plane p);

/**
  Returns a point that creates a normal vector.
  @param atom1 - the first atom in line (start).
  @param atom2 - the middle atom in the line.
  @param atom3 - the end of the line.
  @return - point that can be used to create an normal vector with the line)
*/
Point getNormal(Atom atom1, Atom atom2, Atom atom3);

/**
  returns the distance between atom0 and the line generated by atom1, atom2 and atom3
*/
double getDistance(Atom atom0, Atom atom1, Atom atom2, Atom atom3);

/**
  Translates the atom by x, y, z
  @param atom - the atom to be translated.
  @param x - the distance in the x direction to be translated.
  @param y - the distance in the y direction to be translated.
  @param z - the distance in the z direction to be translated.
  @return - a copy of atom that has been translated accordingly.
*/
Atom translateAtom(Atom atom, double x, double y, double z);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the x axis.
*/
Atom rotateAboutX(Atom atom, double theta);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the y axis.
*/
Atom rotateAboutY(Atom atom, double theta);

/**
  @param atom - the atom to be rotated
  @param theta - the distance in degrees to be rotated.
  @return - atom rotated theta degrees about the z axis.
*/
Atom rotateAboutZ(Atom atom, double theta);


/**
  Rotates atom1 about atom2 in the plane defined by atom1, atom2 and atom3
  theta degrees.
  @param atom1 - the atom to be rotated.
  @param atom2 - the atom about which atom1 will be rotated.
  @param atom3 - the atom used to complete the plane.
  @param theta - the number of degrees to be rotated.
  @return - atom that has been rotated theta degrees in the plane defined
  by the locations of atom1, atom2 and atom3.
  If these points are collinear then the atom will be rotated about the 
  vector created by atom2 and (atom2.x + 1, atom2.y, atom2.z)
*/
Atom rotateAtomInPlane(Atom atom1, Atom atom2, Atom atom3, double theta);

/**
  Rotates atom1 about the vector defined by atom3 and atom2.
  @param atom1 - the atom to be rotated.
  @param atom2 - the atom defining the start point of the vector.
  @param atom3 - the atom defining the end point of the vector.
*/
Atom rotateAtomAboutVector(Atom atom1, Atom atom2, Atom atom3, double theta);

#endif //GEOMETRICUTIL_H
