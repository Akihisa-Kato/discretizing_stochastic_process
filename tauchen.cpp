// This is C++ code for discretizing AR(1) process by using Tauchen method (1986)
// 2017/10/17
//
// Model:
// z' = rho*z+e  where  e~N(0,(sig_e)^2)
//
// INPUTS (through Terminal)
// N:     Number of grids
// m:     Max number of std devs from mean
// rho:   Coeffecient of AR(1) process
// sig_e: Std of innovation, e, in AR(1) process
//
// OUTPUT
// vZ:    grids
// mPI:   transition prob. matrix

#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

double phi(double x); // Standard normal CDF (defined in section 5 in this file)

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//                          MAIN PART
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{

  ///////////////////////////////////////////////////////////
  // 1. Parameters
  ///////////////////////////////////////////////////////////

  // Ex. Standard quarterly calibration for pruductivity
  // process of the US economy

  int     N = 5;          // Number of grids, N
  double  m = 3;          // Max number of std devs from mean
  double  rho = 0.95;     // AR(1) coefficient
  double  sig_e = 0.007;  // std of innovation


  ///////////////////////////////////////////////////////////
  // 2. Generating grids
  ///////////////////////////////////////////////////////////

  double  sig_z = sqrt(pow(sig_e,2)/(1-pow(rho,2)));  // std of AR(1) process r.v.
  double  z_max = m*sig_z;                            // Max grid
  double  z_min = -z_max;                             // Min grid
  double  d = (z_max-z_min) / (N-1);                  // incriment of grid
  double vZ[1][N];                                    // Grid vector 1*N

  vZ[0][0]=z_min;
  vZ[0][N-1]=z_max;

  for (int i=1; i<N-1; i++)
  {
    vZ[0][i] = vZ[0][i-1]+d;
  }

  ////////////////////////////////////////////////////////////
  // 3. Transition matrix
  ////////////////////////////////////////////////////////////

  double mPI[N][N]; // Transition matrix N*N

  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      if (j==0)
      {
        mPI[i][j]=phi((vZ[0][0]+d/2-rho*vZ[0][i])/sig_e);       // If z' is min grid
      }
      else if (j==N-1)
      {
        mPI[i][j]=1-phi((vZ[0][N-1]-d/2-rho*vZ[0][i])/sig_e);   // z' is the max grid
      }
      else
      {
        mPI[i][j]=phi((vZ[0][j]+d/2-rho*vZ[0][i])/sig_e)-phi((vZ[0][j]-d/2-rho*vZ[0][i])/sig_e);
      }

    }
  }

  ////////////////////////////////////////////////////////////
  // 4. Display result with 4 decimals
  ////////////////////////////////////////////////////////////

  std::cout << "Grids are\n";
  for (int i=0; i<N; i++)
  {
    std::cout << std::fixed << std::setprecision(4) << vZ[0][i] << " ";
  }

  std::cout << "\n \nTransition Prob Matrix is \n";
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<N; j++)
    {
      if (j==N-1)
      {
        std::cout << std::fixed << std::setprecision(4) << mPI[i][j] << '\n';
      }
      else
      {
        std::cout << std::fixed << std::setprecision(4) << mPI[i][j] << ' ';
      }
    }
  }

  return 0;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//                          SUB PART
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
// 5. Define CDF function (from https://www.johndcook.com/blog/cpp_phi/)
///////////////////////////////////////////////////////////////

double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}
