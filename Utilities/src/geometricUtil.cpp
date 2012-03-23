#include "geometricUtil.h"

Point createPoint(double X, double Y, double Z){
    Point p;
    p.x = X;
    p.y = Y;
    p.z = Z;
    return p;
}

Plane createPlane(Point p1, Point p2, Point p3){
    Plane p;
    p.point1 = p1;
    p.point2 = p2;
    p.point3 = p3;
    return p;
}

Atom getAtom(vector<Atom> atoms, unsigned long atomID){
    for(int i = 0; i < atoms.size(); i++){
        if(atoms[i].id == atomID)
            return atoms[i];
    }
    return createAtom(-1,-1,-1,-1);
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
    
    double numerator = pow(d3_1, 2) - pow(d1_2, 2);
    double denominator = 2 * d1_2 * d2_3; 

    double radians = acos(numerator / denominator);

    return radiansToDegrees(radians);

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

Atom rotateAtom(Atom atom1, Atom atom2, Atom atom3, double theta){
    //Translate atom2 to the origin
    //translate (-atom2.x, -atom2.y, -atom2.z)
    atom1 = translateAtom(atom1, -atom2.x, -atom2.y, -atom2.z);

    //rotate thetaX about x axis into xz plane
         // temporary atom that is directly below atom1 on the xz plane.
         Atom atomA = createAtom(-1, atom1.x, 0, atom1.z);
         double thetaX = getAngle(atom1, atom2, atomA);       
    atom1 = rotateAboutX(atom1, thetaX);
    
    //rotate about y axis thetaY to align with z axis
        // temporary atom that is on the z axis
        Atom atomB = createAtom(-1, atom1.x, 0, 0);
        double thetaY = getAngle(atom1, atom2, atomB);
    atom1 = rotateAboutY(atom1, thetaY);

    //rotate theta about the z axis
    atom1 = rotateAboutZ(atom1, theta);
    
    //rotate -thetaY about y axis
    atom1 = rotateAboutY(atom1, -thetaY);

    //rotate -thetaX out of xz plane
    atom1 = rotateAboutX(atom1, -thetaX);
    
    //translate (atom2.x, atom2.y, atom2.z)
    atom1 = translateAtom(atom1, atom2.x, atom2.y, atom2.z);
 
    return atom1;
}

