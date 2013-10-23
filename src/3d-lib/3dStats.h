/* @(#)3dStats.h
 */

#include <sstream>
#include <map>
#include <set>
using namespace std;

#ifndef _3DSTATS_H
#define _3DSTATS_H 1


extern void avgCoord(map<string, float**>::iterator it1, map<string, float**>::iterator it2,
		     int size, set<string> modelList, bool add_first, float** avg);
extern void rmsdRMSD(float** xyzA, float** xyzB, int size, float thres, 
		     int &eqv, float &rms, float &drms);
extern void consistency(float** xyzA, float** xyzB, int size, float thres, 
			int * &cons_list);


extern float findCenrtroid (map<string, float**>::iterator it1, float** avg, int size);
extern float** populateMap(int size, float** xyz);

#endif /* _3DSTATS_H */

