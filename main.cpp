#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>

// For RNG
//#include <stdlib.h> // srand, rand
//#include <time.h> // time
#include <random> // for generator and uniform distribution, need C++11, so remember flag -std=c++11

//.................
// Global Constants
//'''''''''''''''''
const double pi = 4.*atan(1.); // Portable definition of pi

const double tau = 1;
const double L = 0.5;
const double alpha = 0.25;
const double DeltaU = 1;

const double T = 1;
const double kB = 1;

const double systemSize = 1;
const size_t n = 201; //system size/resolution must be odd to include 0
double t = 1.77;
const size_t m = 200000; // time steps
const double dt = 1;
const double D = (kB*T)/DeltaU;


//.................
// Function Headers
//'''''''''''''''''
double potential(double x, double t); // tau, L, alpha, Delta_U)
double force(double x, double t);
double gaussianRNG();
std::vector<double> gaussianRNGArray(size_t m); //takes in number of time steps
double langevinEuler(double x, double t);
std::vector<double> make1DBox(const size_t n);
void createPotentialOutput(); // For plotting (testing) the potential, uses make1Dbox()
void createForceOutput(); // For plotting (testing) the force, uses make1Dbox()
void create2DRandomOutput(); // For plotting (testing) the our RNG, uses make1Dbox()
void create3DRandomOutput(); // For plotting (testing) the our RNG, uses make1Dbox()

// ???????????????????????????????????????????????????????????????????????????????????????????????????????
// CAN I AVOID HAVING THESE GLOBALLY DEFINED TO ALLOW IT TO GAIN ACCESS TO MY gaussianRNG()-FUNCTION??????
// ???????????????????????????????????????????????????????????????????????????????????????????????????????
// Activate RNG and create uniform distribution U(0,1):
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

//......
// Main
//''''''
int main(int argc, char** argv){


   std::ofstream outStream("Output/gaussianRandom.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'gaussianRandom.txt' opening failed. \n ";
      exit(1);
   }
   // Creating gaussian random numbers
   std::vector<double> gaussian_numbers;
   gaussian_numbers.resize(m);
   gaussian_numbers = gaussianRNGArray(m);

   // Creating gaussian plot for comparison
   std::vector<double> gaussian_curve;
   gaussian_curve.resize(m);
   std::vector<double> gaussian_xrange;
   gaussian_xrange.resize(m);
   double curve_interval = 6;
   for(int i=0; i<m; ++i){
      gaussian_xrange[i] = i*curve_interval/m-curve_interval/2.0;
      gaussian_curve[i] = 1/sqrt(2*pi)*exp( -(pow(gaussian_xrange[i],2.0))/2.0 );
   }

   for(int i=0; i<gaussian_numbers.size(); ++i){
      outStream << gaussian_numbers[i] << " " << gaussian_xrange[i] << " " << gaussian_curve[i] << std::endl;
   }
   outStream.close();

   return 0;
}


//.........................
// Function Implementations
//'''''''''''''''''''''''''
double potential(double x, double t){
   if(  t - floor(t/tau)*tau < 0.75*tau ){ // OFF is more likely than ON, so this if comes first
      return 0;
   }
   else{
      double loc_x = x-floor(x/L)*L;
      if( loc_x < alpha*L)
         return (DeltaU*loc_x)/(alpha*L);
      else
         return ( DeltaU*(L-loc_x) )/( L*(1-alpha) );
   }
}

double force(double x, double t){
   if(  t - floor(t/tau)*tau < 0.75*tau ){ // OFF is more likely than ON, so this if comes first
      return 0;
   }
   else{
      double loc_x = x-floor(x/L)*L;
      if( loc_x < alpha*L)
         return (DeltaU)/(alpha*L);
      else
         return -DeltaU/( L*(1-alpha) );
   }
}

std::vector<double> gaussianRNGArray(size_t m){ //takes in number of time steps m=even
   double xi_1, xi_2;
   std::vector<double> gaussian_numbers;
   gaussian_numbers.resize(m);
   for (int i=0; i<m; i+=2){
      xi_1 = distribution(generator);
      xi_2 = distribution(generator);
      gaussian_numbers[ i ] = sqrt( -2*log(xi_1) )*cos(2*pi*xi_2);
      gaussian_numbers[i+1] = sqrt( -2*log(xi_1) )*sin(2*pi*xi_2);
   }
   return gaussian_numbers;
}

double langevinEuler(double x, double t, double xi){
   x = x - force(x,t)*dt + sqrt(2*D*dt)*xi;
   return x;
}

double gaussianRNG(){
   return distribution(generator);
}

std::vector<double> make1DBox(const size_t n){
   if(n%2 == 0){
      std::cout << "Error in make1DBox():" << std::endl
                << "System size/resolution n must be an odd number to include 0!" << std::endl;
      exit(1);
   }
   std::vector<double> box;
   box.resize(n);
   for (int i=0; i<n; ++i){
      box[i] = i*(systemSize/(n-1))- systemSize/2;
   }

   return box;
}

// Uses make1Dbox()
void createPotentialOutput(){
   // Output to plotPotential.py
   std::vector<double> box;
   box = make1DBox(n);

   std::ofstream outStream("Output/potential.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'potential.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<box.size(); ++i){
      outStream << box[i] << " " << potential(box[i],t) << std::endl;
   }
   outStream.close();
}

// Uses make1Dbox()
void createForceOutput(){
   // Output to plotForce.py
   std::vector<double> box;
   box = make1DBox(n);

   std::ofstream outStream("Output/force.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'force.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<box.size(); ++i){
      outStream << box[i] << " " << force(box[i],t) << std::endl;
   }
   outStream.close();
}

// Uses gaussianRNG()
void create2DRandomOutput(){
   // Output to plotRandom.py
   std::ofstream outStream("Output/random2d.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'random.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<6000; ++i){
      outStream << gaussianRNG() << " " << gaussianRNG() << std::endl;
   }
   outStream.close();
}

void create3DRandomOutput(){
   // Output to plotRandom.py
   std::ofstream outStream("Output/random3d.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'random3d.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<6000; ++i){
      outStream << gaussianRNG() << " "
                << gaussianRNG() << " "
                << gaussianRNG() << std::endl;
   }
   outStream.close();
}
