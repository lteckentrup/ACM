// ACM model form Williams et al., 2005: An improved analysis of forest carbon
// dynamics using data assimilation

#include <iostream>     // std::cout
#include <algorithm>    // std::min
#include <math.h>
#include <iostream>
#include <fstream>

// Define pi
const double pi = atan(1)*4;

// Parameters from appendix of paper 
const double a1 = 2.155;
const double a2 = 0.0142;
const double a3 = 217.9;
const double a4 = 0.980;
const double a5 = 0.155;
const double a6 = 2.653;
const double a7 = 4.309;
const double a8 = 0.060;
const double a9 = 1.062;
const double a10 = 0.0006;

// Parameter estimates from paper
const double t1 = 4.4 * (pow(10, -6));
const double t2 = 0.47;
const double t3 = 0.31;
const double t4 = 0.43;
const double t5 = 2.7 * (pow(10, -3));
const double t6 = 2.06 * (pow(10, -6));
const double t7 = 2.48 * (pow(10, -3));
const double t8 = 2.2 * (pow(10, -2));
const double t9 = 2.65 * (pow(10, -6));

// Carbon pool
const double cf = 58;               // carbon stored in foliage
const double cw = 770;              // carbon stored in wood
const double cr = 102;              // carbon stored in fine roots
const double clit = 40;             // carbon stored in fresh foliar and fine  
                                    // root litter
const double csomwd = 9897;         // Soil organic matter plus woody debris

// Required input values
const double doy = 172;             // Day of the year
const double lat_deg = 30;          // Latitude (degrees)
const double water_pot = 0.7524;    // Max. soil-leaf water potential 
                                    // difference (MPa)
const double Tmin = 249.742;        // Daily minimum Temperature
const double Tmax = 299.95;         // Daily maximum Temperature 
const double Rtot = 1.0276;         // Total plant-soil hydraulic resistance
const double N = 2.7;               // Average foliar N (gNm-2 leaf area)
const double I = 26.0118;           // Irradiance (MJ m-2 day-1)
const double ca = 355;              // Atmospheric CO2 concentration

// Variables that are calculated
double q;                           // ?
double solar_declination;           // Solar declination (radians)
double L;                           // Leaf area index
double e0;                          // ?
double lat_rad;                     // Latitude (radians)
double s;                           // day length (hours)
double gc;                          // ?
double p;                           // ?
double ci;                          // CO2 concentration in leaf 
double G;                           // GPP (g C m-2 day-1)
            
// Define function to calculate GPP
void GPP(double q, double solar_declination, double L, double e0, 
double lat_rad, double s, double gc, double p, double ci, double G){ 

    q = a3 + a4;
    solar_declination = -0.408 * cos(((360 * (doy + 10))/(365)) * (pi/180));   
    L  = cf/111;
    e0  = (a7 * pow(L, 2))/((pow(L,2))+a9);
    lat_rad = lat_deg * (pi/180);
    if ((tan(lat_rad)*tan(solar_declination))>=1) {
        s = 24;
    } else {
        s = 24 * acos((-1 * tan(lat_rad) * tan(solar_declination))) / pi;
    }
    gc = (pow(water_pot,a10))/((0.5 *(Tmax-Tmin)) +(a6 * Rtot));  
    p = ((a1 * N * L)/gc) * exp(Tmax * a8);
    ci = 0.5 * (ca + q - p + sqrt(pow((ca + q - p),2) - ((4 * ca * q) - \
         (p *a3))));
    G = (e0 * I * gc * (ca-ci))/((e0 * I) + (gc * (ca - ci)))*((a2*s)+a5);

    // Print values 
    std::cout << solar_declination << '\n';
    std::cout << L << '\n';
    std::cout << s << '\n';
    std::cout << q << '\n';
    std::cout << G << '\n';

    // Write output to ASCII
    std::ofstream outfile ("acm_output_doy_172_lat_30.txt");
    outfile  << solar_declination << std::endl;
    outfile  << L << std::endl;
    outfile  << s << std::endl;
    outfile  << q << std::endl;
    outfile  << G << std::endl;
    outfile.close();

    return;  
    } 

// Calculate GPP
int main () {
    GPP(q, solar_declination, L, e0, lat_rad, s, gc, p, ci, G);
      
    return 0;
    }