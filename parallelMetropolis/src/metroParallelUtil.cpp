#include "metroParallelUtil.h"

int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrt(discriminant)) / 2;
    return qv;
}

int getYFromIndex(int x, int idx){
    return idx - x;
}
