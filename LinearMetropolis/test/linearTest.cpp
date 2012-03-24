#include "../../Utilities/src/metroUtil.h"
#include "baseTests.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

void testMakePeriodic(){
    double box = 12.0;
    for(int x = box; x <= box; x++){
        //double baseX = make_periodic(x, box);
        double baseX = make_periodic(x, box);
        double newX = makePeriodic(x, box);
        assert(baseX == newX);
    }

    printf("makePeriodicTest succeeded\n");
}

void testWrapBox(){
    double box = 12.0;
    for(int x = -box; x <= box; x++){
        double baseX = wrap_into_box(x, box);
        double newX = wrapBox(x, box);
        assert(baseX == newX);
    }
    printf("testWrapIntoBox succeeded\n");
}

void testCalculateEnergy(){
    double x1 = 1.0;
    double y1 = 1.0;
    double z1 = 1.0;

    double x2 = 2.0;
    double y2 = 2.0;
    double z2 = 2.0;

    double **coords = new double*[2];
    
    coords[0] = new double[3];
    coords[1] = new double[3];

    coords[0][0] = x1;
    coords[0][1] = y1;
    coords[0][2] = z1;
    coords[1][0] = x2;
    coords[1][1] = y2;
    coords[1][2] = z2;
       
    double sigma = 3.624;
    double epsilon = 0.317;
    double box_size[3] = {5.0, 5.0, 5.0};
    double boxSize = 5.0;
    
    Atom *atoms = new Atom[2];
    atoms[0] = createAtom(0, x1, y1, z1, sigma, epsilon);
    atoms[1] = createAtom(1, x2, y2, z2, sigma, epsilon);
    
    double baseEnergy = calculate_energy(coords, 2, box_size, sigma, epsilon);
    double testEnergy = calculateEnergy(atoms, 2, boxSize);
    
    assert(baseEnergy == testEnergy);
    
    printf("testCalculateEnergy succeeded\n");  
}

int main(){
    testMakePeriodic();
    testWrapBox();
    testCalculateEnergy();
}
