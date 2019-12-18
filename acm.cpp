// ACM model form Williams et al., 2005: An improved analysis of forest carbon
// dynamics using data assimilation

#include <iostream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sstream>
#include <cstdlib>
#include "acm.h"

using namespace std;

double GPP(double, double, double, double, double, double, double, double,
           double, struct Point);

vector<string> split (const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

int main () {

    // Required input values
    double lat_deg = 30.;         // Latitude (degrees)
    double ca;                    // Atmospheric CO2 concentration
    double water_pot;             // Max. soil-leaf water potential
                                  // difference (MPa)
    double Tmin_celsius;          // Daily minimum Temperature (Celsius)
    double Tmax_celsius;          // Daily maximum Temperature (Celsius)
    double Rtot;                  // Total plant-soil hydraulic resistance
    double N;                     // Average foliar N (gNm-2 leaf area)
    double I;                     // Irradiance (MJ m-2 day-1)
    double G;                     // GPP (g C m-2 day-1)

    double doy;
    double temp;
   
    struct Point params;
    
    int i, j, crap, cnt;
    string line;
    
    // Define strings for variables from ASCII file
    string DOY;
    string TEMP;
    string TMAX;
    string TMIN;
    string RAD;
    string PSID;
    string CA;
    string RTOT;
    string NIT;
    string LAI;
        
    // Open file     
    ifstream myfile ("dalec_drivers.OREGON.no_obs.dat");     
    
    // Figure out how big the file is...
    if (myfile.is_open()){
        cnt = 0;
        while(true){
            
            // Get one line from the file
            getline(myfile, line);
            
            // possible as ifstream offers a bool cast operator... 
            if(!myfile) 
                break;
            cnt++;
        }
        
        // Rewind the input file
        // Clear error flags
        myfile.clear(); 
        
        // Set the file get pointer back to the beginning 
        myfile.seekg(0, std::ios::beg); 
        
    } else {
        // If the file is not open output
        cout << "Unable to open file"; 
    }
    
    // Read the file.

    if (myfile.is_open()){
        i = 0;
        std::ofstream outfile ("acm_gpp_lat_30.txt");  
        while(i < cnt){
             
            // Get one line from the file
            getline(myfile, line); 
            vector<string> v = split (line, ' ');

            for (j = 0; j < v.size(); j++){
                // Grab different variables
                DOY = v[0];
                TEMP = v[1];
                TMAX = v[2];
                TMIN = v[3];
                RAD = v[4];
                PSID = v[5];
                CA = v[6];
                RTOT = v[7];
                NIT = v[8];
              
                // Convert string to double
                doy = atof(DOY.c_str());
                temp = atof(TEMP.c_str());
                Tmax_celsius = atof(TMAX.c_str());
                Tmin_celsius = atof(TMIN.c_str());
                I = atof(RAD.c_str());
                water_pot = atof(PSID.c_str());
                ca = atof(CA.c_str());
                Rtot = atof(RTOT.c_str());
                N = atof(NIT.c_str());
      
                // Calculate GPP
                G=GPP(lat_deg, ca, water_pot, Tmin_celsius, Tmax_celsius, Rtot, 
                      N, I, doy, params);
                                       
            }   
            
            // Write GPP to ASCII file
            outfile  << G << std::endl;   
            i++;
            
            // Print GPP
            std::cout<<G<<'\n';            
            
        }
        
        // Close input file
        myfile.close(); 
        
        // Close output file
        outfile.close();        
    } else {
        // If the file is not open output
        cout << "Unable to open file"; 
    }

    return(0);
}

// Define function for solar stuff FIXME function is not working yet
double sun(double doy, double lat_deg){

    double solar_declination;           // Solar declination (radians)
    double lat_rad;                     // Latitude (radians)
    double s;                           // Day length (hours)

    // Calculate solar declination    
    solar_declination = -0.408 * cos(((360. * (doy + 10.))/(365.)) * (pi/180.));
    
    // Convert latitude from degrees to radians
    lat_rad = lat_deg * (pi/180.);
    
    // Calculate day length
    if ((tan(lat_rad)*tan(solar_declination))>=1) {
        s = 24.;
    } else {
        s = 24. * acos((-1. * tan(lat_rad) * tan(solar_declination))) / pi;
    }
    
    return(s);
}

// Define function to calculate GPP
double GPP(double lat_deg, double ca, double water_pot, double Tmin_celsius, 
           double Tmax_celsius, double Rtot, double N, double I, double doy, 
           struct Point params){

    // Variables that need to be calculated

    double lai;                         // Leaf area index
    double e0;                          // Canopy-level quantum yield
    double s;                           // Day length (hours)
    double gc;                          // Stomatal conductance
    double p;                           // ?
    double ci;                          // Internal CO2 concentration
    double G;                           // GPP (g C m-2 day-1)
    double q;
    double Tmin_kelvin;                 // Daily minimum Temperature (Kelvin)
    double Tmax_kelvin;                 // Daily maximum Temperature (Kelvin)
    double water_pot_abs;               // Absolute value of Max. soil-leaf 
                                        // water potential difference (MPa)
    
    // Convert temperature from Celsius to Kelvin
    Tmax_kelvin = Tmax_celsius + 273.15;
    Tmin_kelvin = Tmin_celsius + 273.15;    
    
    // Calculate absolute value for water potential
    if (water_pot<0){
        water_pot_abs = water_pot * (-1);      
    }
    else {
        water_pot_abs = water_pot;
        }
    
    // Call function to calculate sun stuff FIXME is not working yet
    s = sun(doy, lat_deg);   
    
    // Calculate q    
    q = params.a3 + params.a4;
   
    // Calculate leaf area index LAI
    lai  = cf/111.;
    
    // Calculate canopy-level quantum yield eo
    e0  = (params.a7 * pow(lai, 2.))/((pow(lai, 2.))+params.a9);
    
    // Calculate stomatal conductance
    gc = (pow(water_pot_abs,params.a10))/((0.5 *(Tmax_kelvin-Tmin_kelvin)) \
         + (params.a6 * Rtot));
    
    // Calculate p
    p = ((params.a1 * N * lai)/gc) * exp(Tmax_kelvin * params.a8);
    
    // Calculate Internal CO2 concentration
    ci = 0.5 * (ca + q - p + sqrt(pow((ca + q - p), 2.) - ((4. * ca * q) - \
         (p *params.a3))));

    // Calculate GPP
    G = (e0 * I * gc * (ca-ci))/((e0 * I) + (gc * (ca - ci)))*((params.a2*s) + \
         params.a5);
        
    return (G);
}
