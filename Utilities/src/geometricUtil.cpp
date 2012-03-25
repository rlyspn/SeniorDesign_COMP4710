#include "geometricUtil.h"

Point createPoint(double X, double Y, double Z){
    Point p;
    p.x = X;
    p.y = Y;
    p.z = Z;
    return p;
}

Plane createPlane(Atom a1, Atom a2, Atom a3){
    Plane p;
    p.atom1 = a1;
    p.atom2 = a2;
    p.atom3 = a3;
    return p;
}

Atom getAtom(vector<Atom> atoms, unsigned long atomID){
    for(int i = 0; i < atoms.size(); i++){
        if(atoms[i].id == atomID)
            return atoms[i];
    }
    return createAtom(-1,-1,-1,-1);
}

Bond getBond(vector<Bond> bonds, unsigned long a1, unsigned long a2){
    for(int i = 0; i < bonds.size(); i++){
        if(getOppositeAtom(bonds[i], a1) == a2)
            return bonds[i];
    }
    return createBond(-1, -1, -1.0, false);
}

vector<unsigned long> getAllBonds(vector<Bond> bonds, unsigned long atomID){
    vector<unsigned long> toReturn;

    for(int i = 0; i < bonds.size(); i++){
        unsigned long oppositeAtom = getOppositeAtom(bonds[i], atomID);
        if(oppositeAtom != -1)
            toReturn.push_back(oppositeAtom);
    }
    return toReturn;
}

vector<unsigned long> getIntersection(vector<unsigned long> v1, vector<unsigned long> v2){
    vector<unsigned long> intersection;

    //not effecient but I will be working with small data sets.
    for(int i = 0; i < v1.size(); i++){
        for(int j = 0; j < v2.size(); j++){
            if(v1[i] == v2[j])
                intersection.push_back(v1[i]);
        }
    }
    return intersection;
}

bool isMember(vector<unsigned long> atoms, unsigned long toCheck){
    for(int i = 0; i < atoms.size(); i++){
        if(atoms[i] == toCheck)
            return true;
    }
    return false;
}

double degreesToRadians(double degrees){
    return degrees * PI / 180.0;
}

double radiansToDegrees(double radians){
    return radians * 180 / PI;
}

unsigned long getOppositeAtom(Bond bond, unsigned long atomID){
    if(bond.atom1 == atomID)
        return bond.atom2;
    else if(bond.atom2 == atomID)
        return bond.atom1;
    else
        return -1;
}

unsigned long getOppositeAtom(Angle angle, unsigned long atomID){
    if(angle.atom1 == atomID)
        return angle.atom2;
    else if(angle.atom2 == atomID)
        return angle.atom1;
    else
        return -1;
}

unsigned long getOppositeAtom(Dihedral dihedral, unsigned long atomID){
    if(dihedral.atom1 == atomID)
        return dihedral.atom2;
    else if(dihedral.atom2 == atomID)
        return dihedral.atom1;
    else
        return -1;
}

unsigned long getCommonAtom(vector<Bond> bonds, unsigned long atom1,
        unsigned long atom2){
    vector<unsigned long> atom1Bonds; // atom ids bonded to atom1
    vector<unsigned long> atom2Bonds; // atom ids bonded to atom2

    for(int i = 0; i < bonds.size(); i++){
        unsigned long opp1 = getOppositeAtom(bonds[i], atom1);
        unsigned long opp2 = getOppositeAtom(bonds[i], atom2);
        if(opp1 != -1)
            atom1Bonds.push_back(opp1);
        if(opp2 != -1)
            atom2Bonds.push_back(opp2);
    }

    for(int i = 0; i < atom1Bonds.size(); i++){
        unsigned long currentAtom1 = atom1Bonds[i];
        for(int j = 0; j < atom2Bonds.size(); j++){
            if(currentAtom1 == atom2Bonds[j])
                return currentAtom1;
        }
    }
    
    return -1;
}

double getDistance(Atom atom1, Atom atom2){
    return sqrt( pow(atom1.x - atom2.x, 2) +
                 pow(atom1.y - atom2.y, 2) +
                 pow(atom1.z - atom2.z, 2));
}

double getAngle(Atom atom1, Atom atom2, Atom atom3){
    // uses law of cosines, (d3_1)^2 = (d1_2)^2 - 2(d1_2)(d2_3)*cos(theta)
    // theta = arccos(((d3_1)^2 - (d1_2)^2)/(2(d1_2)(d2_3)))
    double d1_2 = getDistance(atom1, atom2); // distance atom1 to atom2
    double d2_3 = getDistance(atom2, atom3); // distance atom2 to atom3
    double d3_1 = getDistance(atom3, atom1); // distance atom3 to atom1
    
    double numerator = pow(d3_1, 2) - pow(d1_2, 2) - pow(d2_3, 2);
    double denominator = -2 * d1_2 * d2_3; 

    double radians = acos(numerator / denominator);

    double degrees = radiansToDegrees(radians);

    if(degrees > 180){
        degrees -= 180;
    }

    return degrees;
}

Point getNormal(Plane p){
    //normal vector = (point2 - point1) CROSS (point3 - point1)
    double atomAX = p.atom2.x - p.atom1.x;
    double atomAY = p.atom2.y - p.atom1.y;
    double atomAZ = p.atom2.z - p.atom1.z;
    Point atomA = createPoint(atomAX, atomAY, atomAZ);

    double atomBX = p.atom3.x - p.atom1.x;
    double atomBY = p.atom3.y - p.atom1.y;
    double atomBZ = p.atom3.z - p.atom1.z;
    Point atomB = createPoint(atomBX, atomBY, atomBZ);
    
    //calculate the cross product of atomA and atomB
    double normX = atomA.y * atomB.z - atomA.z * atomB.y;
    double normY = atomA.z * atomB.x - atomA.x * atomB.z;
    double normZ = atomA.x * atomB.y - atomA.y * atomB.x;
    return createPoint(normX, normY, normZ); 
}

double getAngle(Plane p1, Plane p2){

    //the normal vectors for each plane defined by a poin
    Point normal1 = getNormal(p1);
    Point normal2 = getNormal(p2);

    Atom a1 = createAtom(-1, normal1.x, normal1.y, normal1.z);
    Atom a2 = createAtom(-1, 0, 0, 0);
    Atom a3 = createAtom(-1, normal2.x, normal1.y, normal1.z);

    return getAngle(a1, a2, a3);

}

Atom translateAtom(Atom atom, double x, double y, double z){
    atom.x += x;
    atom.y += y;
    atom.z += z;

    return atom;
}

Atom rotateAboutX(Atom atom, double theta){
    double thetaRadians = degreesToRadians(theta);
    Atom returnAtom = atom;
    returnAtom.y = atom.y * cos(thetaRadians) + atom.z * sin(thetaRadians);
    returnAtom.z = atom.z * cos(thetaRadians) - atom.y * sin(thetaRadians);
    return returnAtom;
}

Atom rotateAboutY(Atom atom, double theta){
    double thetaRadians = degreesToRadians(theta);
    Atom returnAtom = atom;
    returnAtom.z = atom.z * cos(thetaRadians) + atom.x * sin(thetaRadians);
    returnAtom.x = atom.x * cos(thetaRadians) - atom.z * sin(thetaRadians);
    return returnAtom; 
}

Atom rotateAboutZ(Atom atom, double theta){
    double thetaRadians = degreesToRadians(theta);
    Atom returnAtom = atom;
    returnAtom.x = atom.x * cos(thetaRadians) + atom.y * sin(thetaRadians);
    returnAtom.y = atom.y * cos(thetaRadians) - atom.x * sin(thetaRadians);
    return returnAtom;
}

Atom rotateAtomInPlane(Atom atom1, Atom atom2, Atom atom3, double theta){
    //Create a plane.
    Plane atomPlane = createPlane(atom1, atom2, atom3);
    //Find the normal vector.
    Point normal;
    if(getAngle(atom1, atom2, atom3) == 0){
        normal = createPoint(1, 0, 0);
    }
    else{
        normal = getNormal(atomPlane);
    }
    Atom vectorEnd = createAtom(-1, atom2.x + normal.x, atom2.y + normal.y,
            atom2.z + normal.z);
    //Rotate about that normal vector 
    return rotateAtomAboutVector(atom1, atom2, vectorEnd, theta);
}

Atom rotateAtomAboutVector(Atom atom1, Atom atom2, Atom atom3, double theta){
    printf("Rotating atom %f degrees.\n", theta);
    theta *= -1;
    printf("Rotating atom %f degrees.\n", theta);
    //Translate all atoms so that atom2 is at the origin.
    //The rotation axis needs to pass through the origin
    atom1 = translateAtom(atom1, -atom2.x, -atom2.y, -atom3.z);
    atom2 = translateAtom(atom2, -atom2.x, -atom2.y, -atom3.z);
    atom3 = translateAtom(atom3, -atom2.x, -atom2.y, -atom3.z);

    //find the angle between the vector and xz plane
    Atom xzVector = createAtom(-1, 1, 0, 0);
    double xzAngle = getAngle(atom3, atom2, xzVector);  
    printf("xzAngle = %f\n", xzAngle);
    //rotate about z axis so that vector is parallel to xz plane
    atom1 = rotateAboutZ(atom1, xzAngle);
    //atom2 should not change because atom2 is at the origin
    //atom2 = rotateAboutZ(atom2, xzAngle);
    atom3 = rotateAboutZ(atom3, xzAngle);

    //find the angle between the vector and the z axis
    Atom zAxis = createAtom(-1, 0, 0, 1);
    double zAngle = getAngle(atom3, atom2, zAxis);
    printf("zAngle = %f\n", zAngle);
    //rotate about y axis so that the vector is parallel to z axis
    atom1 = rotateAboutY(atom1, zAngle);
    atom2 = rotateAboutY(atom2, zAngle);
    atom3 = rotateAboutY(atom3, zAngle);

    //rotate atom1 theta about the z axis.
    atom1 = rotateAboutZ(atom1, theta);

    //invert rotation about y axis
    atom1 = rotateAboutY(atom1, -zAngle);
    atom2 = rotateAboutY(atom2, -zAngle);
    atom3 = rotateAboutY(atom3, -zAngle);

    //invert rotation about z axis
    atom1 = rotateAboutZ(atom1, -xzAngle);
    //atom2 = rotateAboutZ(atom2, -xzAngle);
    atom3 = rotateAboutZ(atom3, -xzAngle);
    
    //invert translation to origin
    atom1 = translateAtom(atom1, atom2.x, atom2.y, atom3.z);
    atom2 = translateAtom(atom2, atom2.x, atom2.y, atom3.z);
    atom3 = translateAtom(atom3, atom2.x, atom2.y, atom3.z);

    //the inversions for atoms 2 and 3 are not neccesary b/c of copy by value.

    return atom1;
}
