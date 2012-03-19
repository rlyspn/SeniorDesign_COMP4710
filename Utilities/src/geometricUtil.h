#ifndef GEOMETRICUTIL_H
#define GEOMETRICUTIL_H

#include "metroUtil.h"

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

/**
  Structure representing a geometric plane.
  A plane is described by three points.
*/
struct Plane{
    Point point1;
    Point point2;
    Point point3;
};

/**
  @param p1 - the first point in the plane.
  @param p2 - the second point in the plane.
  @param p3 - the third point in the plane.
  @return - a plane represented by points p1, p2 and p3
*/
Plane createPlane(Point p1, Point p2, Point p3);

/**
  Returns a specific atom from a vector of atoms.
  @param atoms - vector of atoms from which to get the Atom
  @param atomID - the ID of the atom to find.
  @return - the atom with atomID or Atom with parameters = -1 if not found.
*/
Atom getAtom(vector<Atom> atoms, unsigned long atomID);

/**
  Returns the id in the bond that is not atomID
  @param bond - the bond to be checked.
  @param atomID - the atomID that you want the opposite of.
*/
unsigned long getOppositeAtom(Bond bond, unsigned long atomID);

/**
  Returns the id in the angle that is not atomID
  @param angle - the angle to be checked.
  @param atomID - the atomID that you want the opposite of.
*/
unsigned long getOppositeAtom(Angle angle, unsigned long atomID);

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

#endif //GEOMETRICUTIL_H
