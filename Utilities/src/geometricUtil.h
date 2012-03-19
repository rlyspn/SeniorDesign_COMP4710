#ifndef GEOMETRICUTIL_H
#define GEOMETRICUTIL_H

#include "metroUtil.h"

/**
  Returns a specific atom from a vector of atoms.
  @param atoms - vector of atoms from which to get the Atom
  @param atomID - the ID of the atom to find.
  @return - the atom with atomID or NULL if not found.
*/
Atom getAtom(vector<Atom> atoms, unsigned long atomID);

/**
  Returns the id in the bond that is not atomID
  @param bond - the bond to be checked.
  @param atomID - the atomID that you want the opposite of.
*/
unsigned long getOppositeAtom(Bond bond, unsigned long atomID);

#endif //GEOMETRICUTIL_H
