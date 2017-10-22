# This is Julia code for discretizing AR(1) process by using Rouwenhorst method (1995)
# written on 2017/10/17
#
# Model:
# z' = rho*z+e  where e~N(0,(sig_e)^2)
#
# INPUTS
# N:     Number of grids
# rho:   Coeffecient of AR(1) process z'=mew+rho*z+e, e~N(0,sigmasq)
# sig_e: Std of innovation, e, in AR(1) process
#
# OUTPUT
# vZ:    grids
# mPI:   transition prob. matrix

#########################################################################
## 1. Parameters
#########################################################################

# Ex. Standard quarterly calibration for productivity process
# of US economy

const N = 5;            # Number of grids, N
const rho = 0.95;       # AR(1) coefficient
const sig_e = 0.007;    # std of innovation

#########################################################################
## 2. Generating grids && Parameterization
#########################################################################

sig_z = sqrt(sig_e^2/(1-rho^2));    # std of AR(1) r.v.
z_max = sig_z*sqrt(N-1);            # max grid
z_min = -z_max;                     # min grid
p = (1+rho)/2;                      # initial prob
q = p;

vZ = linspace(z_min, z_max, N)      # Grids

#########################################################################
## 3. Transition matrix
#########################################################################

mPI_old = [p 1-p; 1-q q]; # Initial PI
mPI_new = zeros(3,3);

for i in 3:N
    mPI_new  = zeros(i,i);
    mPI_new1 = zeros(i,i);  mPI_new1[1:i-1,1:i-1]    = mPI_old;
    mPI_new2 = zeros(i,i);  mPI_new2[1:i-1,2:i]      = mPI_old;
    mPI_new3 = zeros(i,i);  mPI_new3[2:i,1:i-1]      = mPI_old;
    mPI_new4 = zeros(i,i);  mPI_new4[2:i,2:i]        = mPI_old;

    mPI_new = p*mPI_new1 + (1-p)*mPI_new2 + (1-q)*mPI_new3 + q*mPI_new4;
    mPI_new[2:i-1,:] = 0.5*mPI_new[2:i-1,:]; # devide by 2

    mPI_old = mPI_new;
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
        @printf "%6.4f" mPI_new[i,j]
        print(" ")
    end
    print("\n")
end
