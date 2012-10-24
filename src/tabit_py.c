
#include "tadbit.h"

extern "C" {
    tadbit* tacbit_new(){ return new tadbit(); }
    void tadbit_bar(tadbit* foo){ foo->bar(); }
}
