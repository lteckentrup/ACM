// ACM model form Williams et al., 2005: An improved analysis of forest carbon
// dynamics using data assimilation

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
const double t1 = 4.4 * (pow(10., -6.));
const double t2 = 0.47;
const double t3 = 0.31;
const double t4 = 0.43;
const double t5 = 2.7 * (pow(10., -3.));
const double t6 = 2.06 * (pow(10., -6.));
const double t7 = 2.48 * (pow(10., -3.));
const double t8 = 2.2 * (pow(10., -2.));
const double t9 = 2.65 * (pow(10., -6.));

// Carbon pool
const double cf = 58.;              // carbon stored in foliage
const double cw = 770.;             // carbon stored in wood
const double cr = 102.;             // carbon stored in fine roots
const double clit = 40.;            // carbon stored in fresh foliar and fine
                                    // root litter
const double csomwd = 9897.;        // Soil organic matter plus woody debris

double GPP(double, double, double, double, double, double, double, double, 
           double);

// Calculate GPP
int main () {
    
    // Required input values
    const double lat_deg = 30.;         // Latitude (degrees)
    const double ca = 355.;             // Atmospheric CO2 concentration
    const double water_pot = 0.7524;    // Max. soil-leaf water potential
                                        // difference (MPa)
    const double Tmin = 249.742;        // Daily minimum Temperature
    const double Tmax = 299.95;         // Daily maximum Temperature
    const double Rtot = 1.0276;         // Total plant-soil hydraulic resistance
    const double N = 2.7;               // Average foliar N (gNm-2 leaf area)
    const double I = 26.0118;           // Irradiance (MJ m-2 day-1)
    double G;
    
    std::vector<int> day;
    for( int i = 1; i <= 365; i++ )
        day.push_back( i );

    // Write output to ASCII: Open file
    std::ofstream outfile ("acm_output_lat_30.txt");
    
    for(std::size_t i=0; i<day.size(); ++i){
        G=GPP(lat_deg, ca, water_pot, Tmin, Tmax, Rtot, N, I, day[i]);
        std::cout << G << '\n';
        
        // Write G into file
        outfile  << G << std::endl;
    }  
    
    // Close ASCII file
    outfile.close();

    return 0;
}

// Define function to calculate GPP
double GPP(double lat_deg, double ca, double water_pot, double Tmin, double Tmax,
           double Rtot, double N, double I, double doy){

    // Variables that need to be calculated

    double solar_declination;           // Solar declination (radians)
    double lai;                         // Leaf area index
    double e0;                          // ?
    double lat_rad;                     // Latitude (radians)
    double s;                           // day length (hours)
    double gc;                          // stomatal conductance
    double p;                           // ?
    double ci;                          // Internal CO2 concentration
    double G;                            // GPP (g C m-2 day-1)
    double q;

    q = a3 + a4;
    solar_declination = -0.408 * cos(((360. * (doy + 10.))/(365.)) * (pi/180.));
    lai  = cf/111.;
    e0  = (a7 * pow(lai, 2.))/((pow(lai, 2.))+a9);
    lat_rad = lat_deg * (pi/180.);
    if ((tan(lat_rad)*tan(solar_declination))>=1) {
        s = 24.;
    } else {
        s = 24. * acos((-1. * tan(lat_rad) * tan(solar_declination))) / pi;
    }
    gc = (pow(water_pot,a10))/((0.5 *(Tmax-Tmin)) +(a6 * Rtot));
    p = ((a1 * N * lai)/gc) * exp(Tmax * a8);
    ci = 0.5 * (ca + q - p + sqrt(pow((ca + q - p), 2.) - ((4. * ca * q) - \
         (p *a3))));

    // GPP (g C m-2 day-1)
    G = (e0 * I * gc * (ca-ci))/((e0 * I) + (gc * (ca - ci)))*((a2*s)+a5);

    //return;
    return (G); // ASK MARTIN
    }
