#include "../src/metroUtil.h"
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

}

int main(){
    testMakePeriodic();
    testWrapBox();
}
