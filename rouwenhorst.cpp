// This is C++ code for discretizing AR(1) process by using Rouwenhorst method (1995)
// written on 2017/10/17
// z' = rho*z+e e~N(0,(sig_e)^2)
//
// INPUTS (through Terminal)
// N:     Number of grids
// rho:   Coeffecient of AR(1) process z'=mew+rho*z+e, e~N(0,sigmasq)
// sig_e: Std of innovation, e, in AR(1) process
//
// OUTPUT
// vZ:    grids
// mPI:   transition prob. matrix

#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <armadillo>  // YOU SHOULD INSTALL "armadillo" BEFORE USING THIS FILE


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//                          MAIN PART
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  ///////////////////////////////////////////////////////////
  // 1. Basic Parameters
  ///////////////////////////////////////////////////////////

  // Ex. Standard quarterly calibration for pruductivity
  // process of the US economy

  int     N = 5;          // Number of grids, N
  double  rho = 0.95;     // AR(1) coefficient
  double  sig_e = 0.007;  // std of innovation

  ///////////////////////////////////////////////////////////
  // 2. Generating grids && Parameterization
  ///////////////////////////////////////////////////////////

  double  sig_z = sqrt(pow(sig_e,2)/(1-pow(rho,2)));  // std of AR(1) process r.v.
  double  z_max = sig_z*sqrt(N-1);                    // Max grid
  double  z_min = -z_max;                             // Min grid
  double  p = (1+rho)/2, q = p;                       // prob. parameters

  arma::vec vZ = arma::linspace<arma::vec>(z_min,z_max,N); // Grid vector 1*N

  ////////////////////////////////////////////////////////////
  // 3. Transition matrix
  ////////////////////////////////////////////////////////////

  arma::mat mPI_old = {{p,1-p},{1-q,q}}; // Initial PI

  // for recursive computation
  arma::mat mPI_new(3,3);
  arma::mat mPI_new1(3,3);
  arma::mat mPI_new2(3,3);
  arma::mat mPI_new3(3,3);
  arma::mat mPI_new4(3,3);

  for (int i=3;i<N+1; i++)
  {

    mPI_new.zeros(i,i);
    mPI_new1.zeros(i,i);
    mPI_new2.zeros(i,i);
    mPI_new3.zeros(i,i);
    mPI_new4.zeros(i,i);

    mPI_new1.submat(0,0,i-2,i-2) = mPI_old;
    mPI_new2.submat(0,1,i-2,i-1) = mPI_old;
    mPI_new3.submat(1,0,i-1,i-2) = mPI_old;
    mPI_new4.submat(1,1,i-1,i-1) = mPI_old;

    mPI_new = p*mPI_new1 + (1-p)*mPI_new2 + (1-q)*mPI_new3 + q*mPI_new4;
    mPI_new.submat(1,0,i-2,i-1) = 0.5 * mPI_new.submat(1,0,i-2,i-1); // divide by 2

    mPI_old.resize(i,i);
    mPI_old = mPI_new;

  }

  ////////////////////////////////////////////////////////////
  // 4. Display result with 4 decimals
  ////////////////////////////////////////////////////////////

  std::cout << "Grids are\n";
  for (int i=0; i<N; i++)
  {
    std::cout << std::fixed << std::setprecision(4) << vZ(i) << " ";
  }

  std::cout << "\n \nTransition Prob Matrix is \n";
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      if (j == N-1)
      {
        std::cout << std::fixed << std::setprecision(4) << mPI_new(i,j) << '\n';
      }
      else
      {
        std::cout << std::fixed << std::setprecision(4) << mPI_new(i,j) << ' ';
      }
    }
  }
  return 0;
}
