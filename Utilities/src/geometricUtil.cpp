#include "geometricUtil.h"

Atom getAtom(vector<Atom> atoms, unsigned long atomID){
    for(int i = 0; i < atoms.size(); i++){
        if(atoms[i].id == atomID)
            return atoms[i];
    }
}

unsigned long getOppositeAtom(Bond bond, unsigned long atomID){
    if(bond.atom1 == atomID)
        return bond.atom2;
    else
        return bond.atom1;
}
