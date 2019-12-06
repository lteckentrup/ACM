//#ifndef ACM
//#define ACM
#include <iostream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Define pi
const double pi = atan(1)*4;

// Parameters from appendix of paper
struct Point {
    static const double a1 = 2.155;
    static const double a2 = 0.0142;
    static const double a3 = 217.9;
    static const double a4 = 0.980;
    static const double a5 = 0.155;
    static const double a6 = 2.653;
    static const double a7 = 4.309;
    static const double a8 = 0.060;
    static const double a9 = 1.062;
    static const double a10 = 0.0006;
} params;

double t1 = 4.4 * (pow(10., -6.));
double t2 = 0.47;
double t3 = 0.31;
double t4 = 0.43;
double t5 = 2.7 * (pow(10., -3.));
double t6 = 2.06 * (pow(10., -6.));
double t7 = 2.48 * (pow(10., -3.));
double t8 = 2.2 * (pow(10., -2.));
double t9 = 2.65 * (pow(10., -6.));

// Carbon pool
double cf = 58.;              // carbon stored in foliage
double cw = 770.;             // carbon stored in wood
double cr = 102.;             // carbon stored in fine roots
double clit = 40.;            // carbon stored in fresh foliar and fine
                              // root litter
double csomwd = 9897.;        // Soil organic matter plus woody debris
 
double GPP(double, double, double, double, double, double, double, double,
           double, struct Point);

//#endif
