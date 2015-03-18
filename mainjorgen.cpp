#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random> // for generator and uniform distribution, need C++11, so remember flag -std=c++11
#include <time.h> // clock_t, clock, CLOCKS_PER_SEC

//Trying to parellelize using OpenMP:
// compile:  mpic++ mpisend.cpp -o runmpisend -std=c++11 -fopenmp
// run:   ./runme
// if this doesn't work try:   ./runme -omp OMP_NUM_THREADS=2
//#include <omp.h>
//#include <mpi.h>


//...................
// Input Parameters
//'''''''''''''''''''
// Parameter 1: bool, rachet potential ON(true)/OFF(false) 
// Parameter 2: tau, time period of the flashing rachet potential, [s]
// Parameter 3: DeltaU, the potential strength/amplitude
//
//
//


//.................
// Global Constants
//'''''''''''''''''
const double pi = 4.*atan(1.); // Portable definition of pi

const double tau = 1;
const double L = 2*pow(10,-6); // m
const double alpha = 0.2;
const double ri = 12*pow(10,-9); // m, particle size
const double eta = pow(10,-3); // Pa*s = kg/m*s
const double kBT = 26*pow(10,-3); //meV
const double DeltaU = 0.10*kBT;// eV        //*1.6*pow(10,-19); // J
//const double DeltaU = 80; // eV         //*1.6*pow(10,-19); // J

const double gammai = 6*pi*ri*eta;
const double omega = DeltaU/(L*L*gammai); 
const double D = (kBT)/( DeltaU/( 1.6*pow(10,-19) ) ); // (DeltaU converted back to eV)

const double systemSize = 5;
const size_t n = 20001; //system size/resolution/#-of-lattice-points, must be odd to include 0
const size_t timeSteps = 1000; // time steps
//const double dt = 0.001;
const double dt = 0.001;


//.................
// Function Headers
//'''''''''''''''''
double potential(double x, double t); // tau, L, alpha, Delta_U)
double force(double x, double t);
std::vector<double> createGaussianRNGArray(size_t m); //takes in number of time steps
double langevinEuler(double x, double t, double xi);
std::vector<double> make1DBox(const size_t n);
void outputPotential(double t); // For plotting (testing) the potential, uses make1Dbox()
void outputForce(double t); // For plotting (testing) the force, uses make1Dbox()
void outputGaussianDist(size_t m); // For plotting (testing) the force, uses make1Dbox()
void outputTimeCriteria();

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

   // Start clock, runtime
   //double tStart, tStop;              // With OMP
   //tStart = omp_get_wtime();          // With OMP
   clock_t tics;                       // In General(without OMP)
   tics = clock();                     // In General(without OMP)

   double t = 0;
   double potentialON = 0.74*tau;

   outputPotential(potentialON);
   //// Plotting Trajectory
   //std::ofstream outStream("Output/trajectory.txt");
   //if (outStream.fail()){
      //std::cout << "Outputfile 'trajectory.txt' opening failed. \n ";
      //exit(1);
   //}

   //// Making Gaussian Distributed Random numbers
   //std::vector<double> xiArray;
   //xiArray = createGaussianRNGArray(timeSteps);

   //// Initial position
   //double x=0;

   //// Iterate using Euler scheme
   //outStream << x << " " << t << std::endl; // Output values converted to their physical units
////#pragma omp parallel for schedule (static)
   //for (int m = 1; m < timeSteps+1; ++m){
      //t = double(m)*dt;
      ////x = langevinEuler(x, t, xiArray[m]);
      //x = langevinEuler(x, potentialON, xiArray[m]);
      //outStream << x << " " << t << std::endl; // Output values converted to their physical units
   //}
   //outStream.close();


   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/distribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      exit(1);
   }
   size_t numberOfParticles =100000;
   double x;

   // Declaring vector for gaussian distributed random numbers
   std::vector<double> xiArray;

   for(int i=0; i<numberOfParticles; ++i){
      // Initialize gaussian distributed random kicks for particle i:
      xiArray = createGaussianRNGArray(timeSteps);
      // Find final position using Euler scheme:
      for (int m = 0; m < timeSteps; ++m){
         t = double(m)*dt;
         x = langevinEuler(x, potentialON, xiArray[m]); // REMEMBER TO CHANGE potentialON to t  !!!!!!!!!!
      }
      outStream << x << std::endl;
   }
   outStream.close();


   //// Plotting final distribution of potential energy of the particles
   //std::ofstream outStream("Output/energydistribution.txt");
   //if (outStream.fail()){
      //std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      //exit(1);
   //}
   //size_t numberOfParticles =1000;
   //double x;

   //// Declaring vector for gaussian distributed random numbers
   //std::vector<double> xiArray;

   //// Making Boltzmann distribution for comparison:
   //std::vector<double> boltz_dist;
   //boltz_dist.resize(numberOfParticles);
   //std::vector<double> rangeU;
   //rangeU.resize(numberOfParticles);
   //for(int i=0; i<numberOfParticles; ++i){
      //rangeU[i] = double(i)/( double(numberOfParticles)-1.0 );
      //boltz_dist[i] = exp(-rangeU[i])/
                      //( kBT*(exp(DeltaU/kBT)-1) );   //REMEMBER TO CHANGE potentialON to t!!!!!
   //}
   
   
   //for(int i=0; i<numberOfParticles; ++i){
      //// Initialize gaussian distributed random kicks for particle i:
      //xiArray = createGaussianRNGArray(timeSteps);
      //// Find final position using Euler scheme:
      //for (int m = 0; m < timeSteps; ++m){
         //t = double(m)*dt;
         //x = langevinEuler(x, potentialON, xiArray[m]); // REMEMBER TO CHANGE potentialON to t  !!!!!!!!!!
      //}
      //outStream << potential(x, potentialON) << " " << rangeU[i] << " " << boltz_dist[i] << std::endl;
      ////std::cout << potential(x, potentialON) << " " << rangeU[i] << " " << boltz_dist[i] << std::endl;
      ////outStream << x << std::endl;
   //}
   //outStream.close();


   // Get final time and Print Runtime
   std::cout << "runtime: " << (double)(clock()-tics)/CLOCKS_PER_SEC << " seconds" << std::endl; // without OMP
   //tStop = omp_get_wtime();                                                // With OMP
   //std::cout << "Runtime: " << tStop-tStart << " seconds" << std::endl;    // With OMP
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
      double loc_x = x-floor(x);
      if( loc_x < alpha)
         return loc_x/alpha;
      else
         return (1-loc_x)/(1.0-alpha);
   }
}

double force(double x, double t){
   if(  t - floor(t/tau)*tau < 0.75*tau ){ // OFF is more likely than ON, so this if comes first
      return 0;
   }
   else{
      double loc_x = x-floor(x);
      if( loc_x < alpha)
         return 1/alpha;
      else
         return -1/(1.0-alpha);
   }
}
//// FOLLOWING FUNCTIONS ARE NOT SCALED:
//double potential(double x, double t){
   //if(  t - floor(t/tau)*tau < 0.75*tau ){ // OFF is more likely than ON, so this if comes first
      //return 0;
   //}
   //else{
      //double loc_x = x-floor(x/L)*L;
      //if( loc_x < alpha*L)
         //return (DeltaU*loc_x)/(alpha*L);
      //else
         //return ( DeltaU*(L-loc_x) )/( L*(1.0-alpha) );
   //}
//}
//double force(double x, double t){
   //if(  t - floor(t/tau)*tau < 0.75*tau ){ // OFF is more likely than ON, so this if comes first
      //return 0;
   //}
   //else{
      //double loc_x = x-floor(x/L)*L;
      //if( loc_x < alpha*L)
         //return (DeltaU)/(alpha*L);
      //else
         //return -DeltaU/( L*(1.0-alpha) );
   //}
//}

std::vector<double> createGaussianRNGArray(size_t m){ //takes in number of time steps m=even
   double xi_1, xi_2;
   std::vector<double> gaussian_numbers;
   gaussian_numbers.resize(m);
//#pragma omp parallel for schedule (static)
   for (int i=0; i<m; i+=2){
      xi_1 = distribution(generator);
      xi_2 = distribution(generator);
      gaussian_numbers[ i ] = sqrt( -2*log(xi_1) )*cos(2*pi*xi_2);
      gaussian_numbers[i+1] = sqrt( -2*log(xi_1) )*sin(2*pi*xi_2);
   }
   return gaussian_numbers;
}

double langevinEuler(double x, double t, double xi){
//   std::cout << "x = " << x << ", t = " << t << ", xi = " << xi << ", force = " << force(x,t) << ", dt = " << dt << ", force*dt = " << force(x,t)*dt <<  std::endl;
   x = x - force(x,t)*dt + sqrt(2*D*dt)*xi;
   return x;
}

// This function is mainly for the Plotting Functions below
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
void outputPotential(double t){
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
void outputForce(double t){
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

// Uses createGaussianRNGARRAY()
void outputGaussianDist(size_t m){
   std::ofstream outStream("Output/gaussianRandom.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'gaussianRandom.txt' opening failed. \n ";
      exit(1);
   }
   // Creating gaussian random numbers
   std::vector<double> gaussian_numbers;
   gaussian_numbers.resize(m);
   gaussian_numbers = createGaussianRNGArray(m);

   // Creating gaussian plot for comparison
   std::vector<double> gaussian_curve;
   gaussian_curve.resize(m);
   std::vector<double> gaussian_xrange;
   gaussian_xrange.resize(m);
   double curve_interval = 6;
   for(int i=0; i<m; ++i){
      gaussian_xrange[i] = i*curve_interval/(m-1)-curve_interval/2.0;
      gaussian_curve[i] = 1/sqrt(2*pi)*exp( -(pow(gaussian_xrange[i],2.0))/2.0 );
   }

   for(int i=0; i<gaussian_numbers.size(); ++i){
      outStream << gaussian_numbers[i] << " " << gaussian_xrange[i] << " " << gaussian_curve[i] << std::endl;
   }
   outStream.close();
}

// WORK iN PROGRESS! NEED TO OUTPUT THE BOLTZMANN DISTRUBUTION TOGETHER WITH THE OUTPUT DATA:
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void outputBoltzmannDist(size_t numberOfParticles){
   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/distribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      exit(1);
   }

   // Boltzmann distribution

   // Simulation output
   std::vector<double> xiArray;
   double x;
   double t;
   for(int i=0; i<numberOfParticles; ++i){
      // Initialize gaussian distributed random kicks for particle i:
      xiArray = createGaussianRNGArray(timeSteps);
      // Find final position using Euler scheme:
      for (int m = 0; m < timeSteps; ++m){
         t = double(m)*dt;
         x = langevinEuler(x, t, xiArray[m]);
      }
      outStream << x << std::endl;
   }

   outStream.close();
}

void outputTimeCriteria(){
   std::ofstream outStream("Output/timecriteria.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'timecriteria.txt' opening failed. \n ";
      exit(1);
   }
   size_t iterations=1000;
   double idt;
   double LHS;
   for(int i = 0; i < iterations; ++i){
      //idt = i*(alpha/(iterations-1));
      idt = (double)i*0.0015/(double)(iterations-1);
      LHS = idt/(alpha*alpha) + (4*sqrt(2*D*idt))/alpha;
      outStream << idt << " " << LHS << std::endl;
   }
   outStream.close();
}
