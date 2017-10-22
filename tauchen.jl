# This is Julia code for discretizing AR(1) process by using Tauchen method (1986)
# written on 2017/10/17
# z' = rho*z+e e~N(0,(sig_e)^2)
#
# INPUTS
# N:     Number of grids
# m:     Max number of std devs from mean
# rho:   Coeffecient of AR(1) process z'=mew+rho*z+e, e~N(0,sigmasq)
# sig_e: Std of innovation, e, in AR(1) process
#
# OUTPUT
# vZ:    grids
# mPI:   transition prob. matrix


## Preliminary

###############################################################
## 0. Define CDF function (from https://www.johndcook.com/blog/cpp_phi/)
###############################################################

function phi(x::Real)
    # constants
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;

    # Save the sign of x
    sign = 1;
    if (x < 0)
        sign = -1;
    end
    x = abs(x)/sqrt(2.0);

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
end


##############################################################
##############################################################
##                          MAIN PART
##############################################################
##############################################################


###########################################################
## 1. Define Parameters
###########################################################

# Ex. Standard quarterly calibration for productivity process
# of US economy

const N = 5;            # Number of grids, N
const m = 3;            # Max number of std devs from mean
const rho = 0.95;       # AR(1) coefficient
const sig_e = 0.007;    # std of innovation

###########################################################
## 2. Generating grids
###########################################################

sig_z = sqrt((sig_e^2)/(1-rho^2));  # std of AR(1) r.v.
z_max = m*sig_z;                    # Max grid
z_min = -z_max;                     # Min grid
d = (z_max-z_min) / (N-1);          # incriment of grids
vZ = linspace(z_min, z_max, N);     # Grid vector 1*N

###########################################################
## 3. Transition matrix
###########################################################

mPI = zeros(N,N); # Transition matrix N*N

for i in 1:N
    for j in 1:N
        if j == 1
            mPI[i,j] = phi((vZ[1]+d/2-rho*vZ[i])/sig_e);       # If z' is min grid
        elseif j == N-1
            mPI[i,j] = 1 - phi((vZ[N]-d/2-rho*vZ[i])/sig_e);   # z' is the max grid
        else
            mPI[i,j] = phi((vZ[j]+d/2-rho*vZ[i])/sig_e) - phi((vZ[j]-d/2-rho*vZ[i])/sig_e);
        end
    end
end

###########################################################
## 4. Display result with 4 decimals
###########################################################
print("Grid points are \n")
for i in 1:N
    @printf "%6.4f" vZ[i]
    print(" ")
end
print("\n")
print("Transition prob. matrix is\n")
for i in 1:N
    for j in 1:N
        @printf "%6.4f" mPI[i,j]
        print(" ")
    end
    print("\n")
end
