#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random> // for generator and uniform distribution, need C++11, so remember flag -std=c++11
#include <time.h> // clock_t, clock, CLOCKS_PER_SEC

// If trying to parellelize using OpenMP:
// compile:  g++ file.cpp -o runfile -std=c++11 -fopenmp
// run:   ./runme
// if this doesn't work try:   ./runme -omp OMP_NUM_THREADS=2
//#include <omp.h>
//#pragma omp parallel for schedule(static) num_threads(2)

  //...................................................//
 //                 Input Parameters                  //
//'''''''''''''''''''''''''''''''''''''''''''''''''''//'
// Parameter 1: bool, rachet potential ON(true)/OFF(false) 
// Parameter 2: tau, time period of the flashing rachet potential, [s]
// Parameter 3: D (double): diffusion constant D = kBT/DeltaU

  //.............................................//
 //            Global Constants                 //
//'''''''''''''''''''''''''''''''''''''''''''''//
const double pi = 4.*atan(1.); // Portable definition of pi

const double L = 2*pow(10,-6); // m
const double alpha = 0.2;
const double ri = 12*pow(10,-9); // m, particle size
const double eta = pow(10,-3); // Pa*s = kg/m*s
const double kBT = 26*pow(10,-3); //26 meV
const double DeltaU = 80; // eV         //*1.6*pow(10,-19); // J

const double gammai = 6*pi*ri*eta;
const double omega = ( DeltaU*1.6*pow(10,-19) )/(L*L*gammai); 

double D = kBT/DeltaU;
bool flashON = true;

const double systemSize = 5;
const size_t n = 20001; //system size/resolution/#-of-lattice-points, must be odd to include 0
const double dt = 0.01;
//const size_t timeSteps = 10000; // time steps
//const double timeSteps = 0.5*alpha*alpha*(1/gammai)*(1/kBT)*(10/dt);
double timeSteps;

  //................................................//
 //             Function Declarations              //
//''''''''''''''''''''''''''''''''''''''''''''''''//'
double potential(const double x, const double t, const double tau);
double force(const double x, const double t, const double tau);
std::vector<double> createGaussianRNGArray(size_t m); //takes in number of time steps m
void langevinEuler(double& x, const double t, const double tau, const double xi);
void fillGaussianRNGArray(std::vector<double>& timeStepSizedVector);
void runSimulationForParticleDistribution(std::vector<double> & particleDistribution, 
                                             const size_t timeSteps, const double tau);
std::vector<double> make1DBox(const size_t n);
void outputPotentialAndForce(const double tau); // For plotting (testing) the potential, uses make1Dbox()
void outputFlash(const double tau);
void outputGaussianDist(size_t m); // For plotting (testing) the force, uses make1Dbox()
void outputDistributions(size_t numberOfParticles, const double tau);
void outputTrajectory(double x0, const double tau);
void outputDriftVelocity(size_t numberOfParticles, const size_t tau_steps, const double tau_max);



// Activate RNG and create uniform distribution U(0,1):
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

  //................................................//
 //                   Main                         //
//''''''''''''''''''''''''''''''''''''''''''''''''//'
int main(int argc, char** argv){

   if (argc < 3) {
      std::cout << "The program needs " << 3 << " imput parameters:" << std::endl
                << "  Parameter 1: bool, rachet potential flashing ON(true)/OFF(false)" << std::endl 
                << "  Parameter 2: tau, time period of the flashing rachet potential, [s]" << std::endl
                << "  Parameter 3: D (double): diffusion constant D = kBT/DeltaU" << std::endl
                << "   (if U=80eV and kBT = 26 meV, then D = 0.000325 )" << std::endl;
      return -1;
   }
   flashON = atoi(argv[1]);
   double tau = atof(argv[2]);
   D = atof(argv[3]);

   timeSteps = size_t( 10*alpha*alpha/(2*D*dt) );
   //timeSteps = size_t(10*alpha*alpha/(2*D*dt)*0.1) + 1;

   std::cout << "flashON=" << flashON << std::endl
             << "tau=" << tau << std::endl
             << "D=" << D << std::endl
             << "timeSteps=" << timeSteps << std::endl;

   // Start clock, runtime
   //double tStart, tStop;              // With OMP
   //tStart = omp_get_wtime();          // With OMP
   clock_t tics;                       // In General(without OMP)
   tics = clock();                     // In General(without OMP)


   //std::cout << "Starting outputPotentialAndForce()..." << std::endl;
   //outputPotentialAndForce(tau);
   //std::cout << "Starting outputFlash()..." << std::endl;
   //outputFlash(tau);
   //std::cout << "Starting outputTrajectory()..." << std::endl;
   //outputTrajectory(0.0, tau);
   //std::cout << "Starting outputGaussianDist()..." << std::endl;
   //outputGaussianDist(10000);
   //std::cout << "Starting outputDistributions()..." << std::endl;
   //outputDistributions(10000, tau);
   //std::cout << "Starting outputDriftVelocity()..." << std::endl;
   //outputDriftVelocity(1000, 100, 15);


   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/avdriftvelocity.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'avdriftvelocity.txt' opening failed. \n ";
      exit(1);
   }

   size_t numberOfParticles = 10000;

   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   double velocitySum;
   double velocityAvg;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   std::vector<double> particleDistribution;
   particleDistribution.resize(numberOfParticles);


   size_t tau_steps = 50;
   double tau_min = 0.5*alpha*alpha/(2*D);
   tau = tau_min;
   double tau_max = 2*alpha*alpha/(2*D);
   std::cout << "omega=" << omega << std::endl;
   std::cout << "tau_steps=" << tau_steps << std::endl;
   std::cout << "tau_min=" << tau_min << std::endl;
   std::cout << "tau_max=" << tau_max << std::endl;
   for(int k=0; k<tau_steps; ++k){
      // Run simulation to obtain the final particle distribution:
      runSimulationForParticleDistribution(particleDistribution, timeSteps, tau);
      // Calculate average drift velocity:
      velocitySum=0;
      for (int i=0; i<particleDistribution.size(); ++i){
         velocitySum += fabs(particleDistribution[i])/(timeSteps*dt);
      }
      velocityAvg = velocitySum/numberOfParticles;
      outStream << tau << " " << velocityAvg << std::endl;
      tau += (tau_max-tau_min)/double(tau_steps);
   }
   outStream.close();


   
   // Get final time and Print Runtime:
   std::cout << "runtime: " << (double)(clock()-tics)/CLOCKS_PER_SEC << " seconds" << std::endl; // without OMP
   //tStop = omp_get_wtime();                                                // With OMP
   //std::cout << "Runtime: " << tStop-tStart << " seconds" << std::endl;    // With OMP

   return 0;
}



  //................................................//
 //            Function Implementations            //
//''''''''''''''''''''''''''''''''''''''''''''''''//'
double potential(const double x, const double t, const double tau){
   if (t - floor(t/tau)*tau < 0.75*tau && flashON == true) {
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

double force(const double x, const double t, const double tau){
   if (t - floor(t/tau)*tau < 0.75*tau && flashON == true) { 
      return 0;
   }
   else {
      double loc_x = x-floor(x);
      if( loc_x < alpha)
         return 1/alpha;
      else
         return -1/(1.0-alpha);
   }
}

inline void langevinEuler(double& x, const double t, const double tau, const double xi) {
   x = x - force(x,t,tau)*dt + sqrt(2*D*dt)*xi;
}

void fillGaussianRNGArray(std::vector<double>& timeStepSizedVector){ //takes in number of time steps m=even
   if (timeStepSizedVector.size()%2 != 0){
      std::cout << "ERROR! fillGaussianRNGAarray() needs an array of even number dimention!" << std::endl;
      exit(1);
   }
   double xi_1, xi_2;
   for (int i=0; i<timeStepSizedVector.size(); i+=2){
      xi_1 = distribution(generator);
      xi_2 = distribution(generator);
      timeStepSizedVector[ i ] = sqrt( -2*log(xi_1) )*cos(2*pi*xi_2);
      timeStepSizedVector[i+1] = sqrt( -2*log(xi_1) )*sin(2*pi*xi_2);
   }
}

void runSimulationForParticleDistribution(std::vector<double> & particleDistribution,
                                                    const size_t timeSteps, const double tau) {
   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   for(int i=0; i<particleDistribution.size(); ++i){
      // Initialize position and gaussian distributed random kicks for particle i:
      x = 0;
      t = 0;
      fillGaussianRNGArray(xiArray);
      // Find final position using Euler scheme:
      for (int m = 0; m < timeSteps; ++m){
         t = double(m)*dt;
         langevinEuler(x, t, tau, xiArray[m]);
      }
      // Stores the final position in the argument vector:
      particleDistribution[i] = x;
   }
}

void outputFlash(const double tau){
   std::ofstream outStream("Output/flash.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'flash.txt' opening failed. \n ";
      exit(1);
   }
      double t;
      for (int m=0; m<timeSteps; m++){
         t = dt*double(m);
         outStream << t << " " << potential(alpha, t, tau) << std::endl;
      }
   outStream.close();
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
void outputPotentialAndForce(const double tau){
   // Output to plotPotential.py
   std::vector<double> box;
   box = make1DBox(n);

   std::ofstream outStream("Output/potential.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'potential.txt' opening failed. \n ";
      exit(1);
   }
      for(int i=0; i<box.size(); ++i)
         outStream << box[i] << " " << potential(box[i], 0.75*tau, tau) << std::endl;
   outStream.close();

   outStream.open("Output/force.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'force.txt' opening failed. \n ";
      exit(1);
   }
      for(int i=0; i<box.size(); ++i)
         outStream << box[i] << " " << force(box[i],0.75*tau, tau) << std::endl;
   outStream.close();
}

// Uses createGaussianRNGARRAY()
void outputGaussianDist(size_t m){
   // Creating gaussian random numbers
   std::vector<double> gaussian_numbers;
   gaussian_numbers.resize(m);
   fillGaussianRNGArray(gaussian_numbers);

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

   // Output the gaussian distributed numbers together with the theoretical gaussian curve:
   std::ofstream outStream("Output/gaussianRandom.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'gaussianRandom.txt' opening failed. \n ";
      exit(1);
   }
      for(int i=0; i<gaussian_numbers.size(); ++i){
         outStream << gaussian_numbers[i] << " " << gaussian_xrange[i] << " " << gaussian_curve[i] << std::endl;
   }
   outStream.close();
}

void outputTrajectory(double x0, const double tau){
   // Initial position and time
   double x=x0;
   double t=0;

   // Making Gaussian Distributed Random numbers:
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);
   fillGaussianRNGArray(xiArray);

   std::ofstream outStream("Output/trajectory.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'trajectory.txt' opening failed. \n ";
      exit(1);
   }
      // Iterate using Euler scheme
      outStream << x << " " << t << std::endl; // Output values converted to their physical units
      for (int m = 1; m < timeSteps+1; ++m){
         t = double(m)*dt;
         langevinEuler(x, t, tau, xiArray[m]);
         outStream << x << " " << t << std::endl; // Output values converted to their physical units
      }
   outStream.close();
}

void outputDistributions(size_t numberOfParticles, const double tau){
   /* Declaring vector for random numbers: */
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   /* Declaring vector for distribution: */
   std::vector<double> distribution;
   distribution.resize(numberOfParticles);

   /* Run simulation and get final distribution: */
   runSimulationForParticleDistribution(distribution, timeSteps, tau);

   /* Output Distribution: */
   std::ofstream outStream("Output/distribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      exit(1);
   }
      for (int i = 0; i<distribution.size(); ++i)
         outStream << distribution[i] << std::endl;
   outStream.close();

   /* Output Energy Distribution: */
   outStream.open("Output/energydistribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'energydistribution.txt' opening failed. \n ";
      exit(1);
   }
      // Making Boltzmann distribution for comparison:
      std::vector<double> boltz_dist;
      boltz_dist.resize(numberOfParticles);
      std::vector<double> rangeU;
      rangeU.resize(numberOfParticles);
      for(int i=0; i<numberOfParticles; ++i){
         rangeU[i] = double(i)/( double(numberOfParticles)-1.0 );
         boltz_dist[i] = exp( -rangeU[i]/D )/( D*(1-exp(-1.0/D)) );
      }
      // Output data:
      for (int i = 0; i<distribution.size(); ++i)
         outStream << potential(distribution[i], 0.75*tau, tau) << " " 
                   << rangeU[i] << " " << boltz_dist[i] << std::endl;
   outStream.close();
}

void outputDriftVelocity(size_t numberOfParticles, const size_t tau_steps, const double tau_max=15.0) {
   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/avdriftvelocity.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'avdriftvelocity.txt' opening failed. \n ";
      exit(1);
   }
   //size_t numberOfParticles = 1000;

   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   double velocitySum;
   double velocityAvg;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   double tau = 0;
   //tau_steps = 1000;
   //tau_max = 15;
   std::cout << "omega=" << omega << std::endl;
   std::cout << "tau_steps=" << tau_steps << std::endl;
   for(int k=0; k<tau_steps; ++k){
//      std::cout << "tau = " << tau << std::endl;
      // Update tau:
      // Reset Sum:
      velocitySum=0;
      for(int i=0; i<numberOfParticles; ++i){
         // Initialize position and gaussian distributed random kicks for particle i:
         x = 0;
         fillGaussianRNGArray(xiArray);
         // Find final position using Euler scheme:
         for (int m = 0; m < timeSteps; ++m){
            t = double(m)*dt;
            langevinEuler(x, t, tau, xiArray[m]);
         }
         velocitySum += fabs(x)/(timeSteps*dt);
      }
      velocityAvg = velocitySum/numberOfParticles;
      outStream << tau << " " << velocityAvg << std::endl;
      tau += tau_max/double(tau_steps);
   }
   outStream.close();
}
