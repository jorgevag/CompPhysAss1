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
const double DeltaU = 10*kBT;// eV        //*1.6*pow(10,-19); // J
//const double DeltaU = 80; // eV         //*1.6*pow(10,-19); // J

const double gammai = 6*pi*ri*eta;
const double omega = ( DeltaU*1.6*pow(10,-19) )/(L*L*gammai); 
//const double D = (kBT)/( DeltaU/( 1.6*pow(10,-19) ) ); // (DeltaU converted back to eV)
double D = kBT/DeltaU; // (DeltaU converted back to eV)
double tau = 1;
bool flashON = true;

const double systemSize = 5;
const size_t n = 20001; //system size/resolution/#-of-lattice-points, must be odd to include 0
const double dt = 0.01;
//const size_t timeSteps = 10000; // time steps
//const double timeSteps = 0.5*alpha*alpha*(1/gammai)*(1/kBT)*(10/dt);
double timeSteps;

  //................................................//
 //             Function Implementations           //
//''''''''''''''''''''''''''''''''''''''''''''''''//'
double potential(const double x, const double t, const double tau);
double force(const double x, const double t, const double tau);
std::vector<double> createGaussianRNGArray(size_t m); //takes in number of time steps m
void langevinEuler(double& x, const double t, const double tau, const double xi);

void fillGaussianRNGArray(std::vector<double>& timeStepSizedVector);
void runSimulationForParticleDistribution
      (std::vector<double> & particleDistribution, const size_t timeSteps, const double tau);

std::vector<double> make1DBox(const size_t n);
void outputPotential(const double tau); // For plotting (testing) the potential, uses make1Dbox()
void outputForce(const double tau); // For plotting (testing) the force, uses make1Dbox()
void outputGaussianDist(size_t m); // For plotting (testing) the force, uses make1Dbox()
void outputTrajectory(double x0, const double tau);
void outputParticleDist(size_t numberOfParticles, const double tau);
void outputEnergyDist(size_t numberOfParticles, const double tau);
void outputDriftVelocity(size_t numberOfParticles, const size_t tau_steps, const double tau_max);
void outputTimeCriteria();

// ???????????????????????????????????????????????????????????????????????????????????????????????????????
// CAN I AVOID HAVING THESE GLOBALLY DEFINED TO ALLOW IT TO GAIN ACCESS TO MY gaussianRNG()-FUNCTION??????
// ???????????????????????????????????????????????????????????????????????????????????????????????????????
// Activate RNG and create uniform distribution U(0,1):
std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

double potentialON = 0.75*tau;

  //................................................//
 //                   Main                         //
//''''''''''''''''''''''''''''''''''''''''''''''''//'
int main(int argc, char** argv){

   if (argc < 3) {
      std::cout << " Need " << 3 << " imput parameters!" << std::endl
                << std::endl
                << "Parameter 1: bool, rachet potential flashing ON(true)/OFF(false)" << std::endl 
                << "Parameter 2: tau, time period of the flashing rachet potential, [s]" << std::endl
                << "Parameter 3: D (double): diffusion constant D = kBT/DeltaU" << std::endl
                << "  (if U=80eV and kBT = 26 meV, then D = 0.000325 )" << std::endl;
      return -1;
   }
   flashON = atoi(argv[1]);
   tau = atof(argv[2]);
   D = atof(argv[3]);


   timeSteps = 10*alpha*alpha/(2*D*dt);
   //const bool flashON = true;
   //const double tau = 1;
   //const double D = kBT/DeltaU;
   std::cout << "flashON=" << flashON << std::endl
             << "tau=" << tau << std::endl
             << "D=" << D << std::endl
             << "timeSteps=" << timeSteps << std::endl;

   // Start clock, runtime
   //double tStart, tStop;              // With OMP
   //tStart = omp_get_wtime();          // With OMP
   clock_t tics;                       // In General(without OMP)
   tics = clock();                     // In General(without OMP)


   //outputPotential(tau);
   //outputForce(tau);
   //outputTrajectory(0.5, tau);
   //outputGaussianDist(1000);
   //outputParticleDist(10000, tau);
   //outputEnergyDist(100000, tau);
   
   //outputDriftVelocity(1000, 100, 15);


   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/avdriftvelocity.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'avdriftvelocity.txt' opening failed. \n ";
      exit(1);
   }

   size_t numberOfParticles = 100;

   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   double velocitySum;
   double velocityAvg;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   std::vector<double> particleDistribution;
   particleDistribution.resize(numberOfParticles);


   double tau = 0;
   size_t tau_steps = 10;
   double tau_max = alpha*alpha/(2*D);
   std::cout << "omega=" << omega << std::endl;
   std::cout << "tau_steps=" << tau_steps << std::endl;
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
      tau += tau_max/double(tau_steps);
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

inline 
void langevinEuler(double& x, const double t, const double tau, const double xi) {
   x = x - force(x,t,tau)*dt + sqrt(2*D*dt)*xi;
}

void fillGaussianRNGArray(std::vector<double>& timeStepSizedVector){ //takes in number of time steps m=even
   double xi_1, xi_2;
   for (int i=0; i<timeStepSizedVector.size(); i+=2){
      xi_1 = distribution(generator);
      xi_2 = distribution(generator);
      timeStepSizedVector[ i ] = sqrt( -2*log(xi_1) )*cos(2*pi*xi_2);
      timeStepSizedVector[i+1] = sqrt( -2*log(xi_1) )*sin(2*pi*xi_2);
   }
}

inline
void runSimulationForParticleDistribution
         (std::vector<double> & particleDistribution, const size_t timeSteps, const double tau) {

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
void outputPotential(const double tau){
   // Output to plotPotential.py
   std::vector<double> box;
   box = make1DBox(n);

   std::ofstream outStream("Output/potential.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'potential.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<box.size(); ++i){
      outStream << box[i] << " " << potential(box[i], 0.75*tau, tau) << std::endl;
   }
   outStream.close();
}

// Uses make1Dbox()
void outputForce(const double tau){
   // Output to plotForce.py
   std::vector<double> box;
   box = make1DBox(n);

   std::ofstream outStream("Output/force.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'force.txt' opening failed. \n ";
      exit(1);
   }
   for(int i=0; i<box.size(); ++i){
      outStream << box[i] << " " << force(box[i],0.75*tau, tau) << std::endl;
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

void outputTrajectory(double x0, const double tau){
   // Plotting Trajectory
   std::ofstream outStream("Output/trajectory.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'trajectory.txt' opening failed. \n ";
      exit(1);
   }

   // Initial position and time
   double x=x0;
   double t=0;
   // Making Gaussian Distributed Random numbers
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);
   xiArray = createGaussianRNGArray(timeSteps);


   // Iterate using Euler scheme
   outStream << x << " " << t << std::endl; // Output values converted to their physical units
//#pragma omp parallel for schedule (static)
   for (int m = 1; m < timeSteps+1; ++m){
      t = double(m)*dt;
      //x = langevinEuler(x, t, xiArray[m]);
      langevinEuler(x, t, tau, xiArray[m]);
      outStream << x << " " << t << std::endl; // Output values converted to their physical units
   }
   outStream.close();
}

void outputParticleDist(size_t numberOfParticles, const double tau){
   // Plotting Final position (distribution) of many particles
   std::ofstream outStream("Output/distribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      exit(1);
   }

   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   for(int i=0; i<numberOfParticles; ++i){
      // Initialize position and gaussian distributed random kicks for particle i:
      x = 0;
      t = 0;
      xiArray = createGaussianRNGArray(timeSteps);
      // Find final position using Euler scheme:
      for (int m = 0; m < timeSteps; ++m){
         t = double(m)*dt;
         langevinEuler(x, t, tau, xiArray[m]);
      }
      outStream << x << std::endl;
   }

   outStream.close();
}

void outputEnergyDist(size_t numberOfParticles, const double tau){
   // Plotting final distribution of potential energy of the particles
   std::ofstream outStream("Output/energydistribution.txt");
   if (outStream.fail()){
      std::cout << "Outputfile 'distribution.txt' opening failed. \n ";
      exit(1);
   }

   // Declaring particle position and vector for gaussian distributed random numbers
   double x;
   double t;
   std::vector<double> xiArray;
   xiArray.resize(timeSteps);

   // Making Boltzmann distribution for comparison:
   std::vector<double> boltz_dist;
   boltz_dist.resize(numberOfParticles);
   std::vector<double> rangeU;
   rangeU.resize(numberOfParticles);
   for(int i=0; i<numberOfParticles; ++i){
      rangeU[i] = double(i)/( double(numberOfParticles)-1.0 );
      boltz_dist[i] = exp( -rangeU[i]/D )/( D*(1-exp(-1.0/D)) );
   }
   
   for(int i=0; i<numberOfParticles; ++i){
      // Initialize position and gaussian distributed random kicks for particle i:
      x = 0;
      xiArray = createGaussianRNGArray(timeSteps);
      // Find final position using Euler scheme:
      for (int m = 0; m < timeSteps; ++m){
         t = double(m)*dt;
         langevinEuler(x, t, tau, xiArray[m]);
      }
      outStream << potential(x, 0.75*tau, tau) << " " << rangeU[i] << " " << boltz_dist[i] << std::endl;
   }

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
         xiArray = createGaussianRNGArray(timeSteps);
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

//void outputTimeCriteria(){
   //std::ofstream outStream("Output/timecriteria.txt");
   //if (outStream.fail()){
      //std::cout << "Outputfile 'timecriteria.txt' opening failed. \n ";
      //exit(1);
   //}
   //size_t iterations=1000;
   //double idt;
   //double LHS;
   //for(int i = 0; i < iterations; ++i){
      ////idt = i*(alpha/(iterations-1));
      //idt = (double)i*0.0015/(double)(iterations-1);
      //LHS = idt/(alpha*alpha) + (4*sqrt(2*D*idt))/alpha;
      //outStream << idt << " " << LHS << std::endl;
   //}
   //outStream.close();
//}
