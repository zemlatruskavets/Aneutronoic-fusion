//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////                                                                  //////
//////                              FUSION                              //////
//////                                                                  //////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////

//                Include headers and assign global variables               //

//////////////////////////////////////////////////////////////////////////////



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



//////////////////////////////////////////////////////////////////////////////

//                      The units are as follows (SI):                      //
//                                                                          //
//                                                                          //
//                      Energy:             J                               //    
//                      Mass:               kg                              //
//                      Temperature:        K                               //
//                      Distance            m                               //
//                      Time                s                               //
//                      Angle               radians                         //
//                      Velocity            m / s                           //
//                      Force               N                               //
//                      Charge              C                               //
//                                                                          //

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//                          Initialize all variables                        //

//////////////////////////////////////////////////////////////////////////////


int      debug = 0;
double   eps0 = (8.99e9);  // (N * m^2) / C^2
char     *element;
double   *coords, *forces, *old_forces, *velocities, *mass, *charge;
int      numatoms = 200;
int      max_numatoms = 3 * numatoms;
int      *type;

double   radius_charge = 100; 
double   theta;
double   rf;
double   pi = 3.141592653;
double   KB = 1.381e-23; // J / K
double   ec = 1.602e-19;  // C
double   amu = 1.66e-27; // kg
double   AR = (pi / 180.0);  // Angles to radians
double   F1 = (((-1.0 / 360.0) * theta) + 0.5);
double   F2 = ((1.0 / 180.0) * theta);
double   d_threshold;


double   UpdateEnergyForces();
void     CompareAnalyticNumericalGradients();
double   IterateDynamics(double dt);
void     WriteCoords(FILE *fp, double *coords);
double   KineticEnergy();
double   Temperature();
double   Pressure(double rf);
void     InitializeCapps();
double   rand_uni();
double   rand_exp(double lambda);
double   rand_gauss();
void     FuseNearbyNuclei();
void     FuseNuclei(int i, int j);
void     RandomTriton(double *x, double *y, double *z, double R);
double   NormalizeQ(double *q0, double *q1, double *q2, double *q3);
double   norm(double q0, double q1, double q2, double q3);
void     RotatePoint(double x, double y, double z, double *new_x, double *new_y, double *new_z, double q0, double q1, double q2, double q3);
void     SetMass(int i, double m);
void     InitializeRandomVelocities(double temperature);
void     ReflectiveBoundaries();
double   DotProduct();

// PROGRAM OPTIONS
enum     {MINIMIZE, DYNAMICS};
int      calculation_type = DYNAMICS;
int      reflective_boundaries = 0;
int      check_fusion = 0;
double   **capp_distance;
int      **capp_time;
int      iteration;
double   dt;



int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Usage: ./Fusion <input.xyz> <radius of background charge> <temperature>\n");
        exit(1);
    }
    
    
//////////////////////////////////////////////////////////////////////////////
    
//      Create various arrays using heap memory in which to store data      //
    
//////////////////////////////////////////////////////////////////////////////
    
    
    // Allocate arrays
    element    = (char   *)  malloc(sizeof(char)       * max_numatoms);  
    coords     = (double *)  malloc(sizeof(double) * 3 * max_numatoms);
    forces     = (double *)  malloc(sizeof(double) * 3 * max_numatoms);
    old_forces = (double *)  malloc(sizeof(double) * 3 * max_numatoms);
    velocities = (double *)  malloc(sizeof(double) * 3 * max_numatoms);
    mass       = (double *)  malloc(sizeof(double) * 3 * max_numatoms);
    type       = (int    *)  malloc(sizeof(int   )     * max_numatoms);
    charge     = (double *)  malloc(sizeof(double)     * max_numatoms);
    
    
    int i, j, k;
    capp_distance = (double **) malloc(sizeof(double *) * max_numatoms);
    capp_time     = (int    **) malloc(sizeof(int    *) * max_numatoms);
    
    for (i = 0; i < max_numatoms; i++)
    {
        capp_distance[i] = (double *) malloc(sizeof(double) * max_numatoms);
        capp_time[i]     = (int    *) malloc(sizeof(int   ) * max_numatoms);
        type[i] = -1;
    }
    
    
    
//////////////////////////////////////////////////////////////////////////////
    
//       Read in data into the various arrays from the geometry file        //
    
//////////////////////////////////////////////////////////////////////////////
    
    
    
    // Read in uniform charge radius
    //sscanf(argv[2], "%lf", &radius_charge);
    
    
    // Read in temperature
    double temperature;
    sscanf(argv[2], "%lf", &temperature);
    
    
    // Open geometry file
    FILE *fp = fopen(argv[1], "r");
    
    
    // Fill out the distance and time arrays
    InitializeCapps();
    
    
    // Open the input file and get all sorts of data
    int EndOfFile;
    i = 0;
    EndOfFile = 1;
    
    while (EndOfFile != -1)
    {
        EndOfFile = fscanf(fp, "%s %lf %lf %lf", &element[i], &coords[3 * i], &coords[3 * i + 1], &coords[3 * i + 2]);
        i++;
    }
    
    
    // Fill out type, charge, and mass and velocity arrays
    for (i = 0; i < numatoms; i++)
    {
        if (element[i] == 'B')
        {
            type[i] = 2;
            charge[i] = 3 * ec;
            mass[3 * i    ] = 11 * amu;
            mass[3 * i + 1] = 11 * amu;
            mass[3 * i + 2] = 11 * amu;
        }
        
        
        else if (element[i] == 'H')
        {
            type[i] = 1;
            charge[i] = 1 * ec;
            mass[3 * i    ] = 1 * amu;
            mass[3 * i + 1] = 1 * amu;
            mass[3 * i + 2] = 1 * amu;
        }
    }  
    
    
    fclose(fp);
    
    
    
    
//////////////////////////////////////////////////////////////////////////////
    
//       Thus far, there is a coordinate array in which each atom is       //
//       characterized by three elements, a type array that describes      //
//       each particle based on whether it is a boron or hydrogen, a       //
//       charge array that assigns each particle a charge, and a mass     //
//       array that assigns each particle a mass, and velocity arrays      //
//       that store the velocities of each particle                        //
    
//////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
//////////////////////////////////////////////////////////////////////////////
    
//           Now, a file is opened in which the data is output              //      
    
//////////////////////////////////////////////////////////////////////////////
    
    
    
    // Initialize velocities
    InitializeRandomVelocities(temperature);
    
    
    fp = fopen("output.xyz", "w");
    
    // Do dynamics
    if (calculation_type == DYNAMICS)
        dt = 2e-18;           
    else
        dt = 1.0;           
    
    void ReflectiveBoundaries();
    // Calculate the energy of the system
    UpdateEnergyForces();  
    //for (theta = 0; theta < 2 * AR; theta += AR)
    for (iteration = 0; iteration < 1e4; iteration++)
    {
        // Calculate the potential and kinetic energies
        double PE = IterateDynamics(dt);
        double KE = KineticEnergy();
        
        // Calculate the temperature and pressure
        double T = Temperature();
        double rf = DotProduct();
        double P = Pressure(rf);
        
        // Make the array storing the boron-hydrogen 
        for (i = 0; i < numatoms; i++)
            for (j = 0; j < i; j++)
                capp_time[i][j]++;
        
        
        // Check for fusion every 50 iterations
        if (check_fusion && iteration % 50 == 0)
        {
            FuseNearbyNuclei();
            InitializeCapps();
        }
        
        // Print various information every 250 iterations
        if (iteration % 250 == 0) 
        {
            double time_step = 1e-18; // 1 attosecond
            printf("%c %i %g %g %g %g %g %e %g\n", element[99], type[99], charge[99], mass[99* 3], coords[99 * 3], velocities[99 * 3], forces[99 * 3], temperature, radius_charge);
            //printf("\n t (s) = %.3e  T = %.3e P = %.3g KE (J) = %.3e  PE (J) = %.3e  total (J) = %.3e ", (double) /*(theta / AR),*/ iteration * time_step, T, P, KE, PE, KE + PE);
            WriteCoords(fp, coords);
        }
        
        void ReflectiveBoundaries();
        void CompareAnalyticNumericalGradients();
        //if (reflective_boundaries) ReflectiveBoundaries();
        
    }
    
    
    fclose(fp);
    
}




//////////////////////////////////////////////////////////////////////////////

//      This sections contains the various functions called in              //      
//      the body of the program                                             //

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////

//            This creates the arrays for the distance and time             //

//////////////////////////////////////////////////////////////////////////////


void InitializeCapps()
{ 
    int i, j;
    for (i = 0; i < numatoms; i++)
    {
        for (j = 0; j < i; j++)
        {
            capp_distance[i][j] = 1e50;
            capp_time[i][j]     = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////

//      This ensures the particles are contained within a certain area      //

//////////////////////////////////////////////////////////////////////////////


void ReflectiveBoundaries()
{
    int i;
    for (i = 0; i < numatoms; i++)
    {
        double x, y, z;
        x = coords[3 * i + 0];
        y = coords[3 * i + 1];
        z = coords[3 * i + 2];
        
        if (x * x + y * y + z * z > pow(1.5 * radius_charge, 2))
        {
            velocities[3 * i + 0] *= -1;
            velocities[3 * i + 1] *= -1;
            velocities[3 * i + 2] *= -1;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////

//       This section assigns velocities and uses a Gamow probability       //
//       type tunneling formulation to assign criteria for fusion           //

//////////////////////////////////////////////////////////////////////////////


// Assign a velocity based on a predefined distribution
void InitializeRandomVelocities(double temperature)
{
    int i;
    double KB = 1.381e-23; // J / K
    for (i = 0; i < numatoms; i++)
    {
        double sigma = sqrt((KB * temperature) / mass[i]);
        velocities[i] = rand_gauss() * sigma;
    }
}


double rand_uni()
{
    return (double) rand() / (double) RAND_MAX;
}

double rand_exp(double lambda)
{
    return -log(rand_uni()) / lambda;
}

double rand_gauss()
{
    double r = 2.0, v1, v2;
    while (r > 1.0)
    {
        v1 = 2.0 * rand_uni() - 1.0;
        v2 = 2.0 * rand_uni() - 1.0;
        r = v1 * v1 + v2 * v2;
    }
    return v1 * sqrt(-2.0 * log(r) / r);
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////                                                                  //////
//////               This section is related to fusion                  //////
//////                                                                  //////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//          This creates a coordinate array for the α particles             //

//////////////////////////////////////////////////////////////////////////////


void FuseNearbyNuclei()
{
    // Check for fusion events based on distances of closest approach
    double alpha = 0.065;
    int i, j, k;
    for (i = 0; i < numatoms; i++)
    {
        for (j = 0; j < i; j++)
        { 
            //double d_threshold = 1;
            double d_threshold = rand_exp(alpha);
            if ((type[i] == 1 && type[j] == 2) || (type[i] == 2 && type[j] == 1))
            {
                if (capp_distance[i][j] < d_threshold)
                {
                    // Rewind simulation to point of closest approach
                    int iter;
                    int numsteps = capp_time[i][j];
                    for (k = 0; k < numatoms; k++)
                        velocities[k] = -velocities[k];
                    
                    for (iter = 0; iter < numsteps; iter++)
                        IterateDynamics(-dt);
                    
                    for (k = 0; k < numatoms; k++)
                        velocities[k] = -velocities[k];
                    
                    iteration -= numsteps;
                    
                    // Reset all points of closest approach
                    InitializeCapps();
                    
                    // Fuse nuclei
                    FuseNuclei(i, j);
                    printf("Fused %i and %i\n", i, j);
                    return;
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//                      This fuses particles together                       //

//////////////////////////////////////////////////////////////////////////////


void FuseNuclei(int i, int j)
{
    // Energy before any changes
    double E_before = UpdateEnergyForces();
    
    // KE of particles that will disappear
    double ke_before = 0;
    int idx;
    for (idx = 3 * i; idx < 3 * i + 3; idx++)
        ke_before += 0.5 * mass[idx] * velocities[idx] * velocities[idx];
    for (idx = 3 * j; idx < 3 * j + 3; idx++)
        ke_before += 0.5 * mass[idx] * velocities[idx] * velocities[idx];
    
    // Centre of mass velocity and position
    double x0, y0, z0, vx0, vy0, vz0;
    
    x0  = (mass[i] * coords[3 * i    ] + mass[j] * coords[3 * j    ]) / (mass[i] + mass[j]);
    y0  = (mass[i] * coords[3 * i + 1] + mass[j] * coords[3 * j + 1]) / (mass[i] + mass[j]);
    z0  = (mass[i] * coords[3 * i + 2] + mass[j] * coords[3 * j + 2]) / (mass[i] + mass[j]);
    
    vx0 = (mass[i] * velocities[3 * i    ] + mass[j] * velocities[3 * j    ]) / (mass[i] + mass[j]);
    vy0 = (mass[i] * velocities[3 * i + 1] + mass[j] * velocities[3 * j + 1]) / (mass[i] + mass[j]);
    vz0 = (mass[i] * velocities[3 * i + 2] + mass[j] * velocities[3 * j + 2]) / (mass[i] + mass[j]);
    
    
    
    // Set position of CM and then split into three α particles
    double R = eps0 * (12 * amu) / sqrt(3) / ((8.7e6) * ec);
    double x[3], y[3], z[3];
    RandomTriton(x, y, z, R);
    
    // Delete two particles and set their array elements to 0
    charge[i] = charge[j] = 0;
    type[i] = type[j] = -1;
    velocities[3 * i + 0] = velocities[3 * j + 0] = 0;
    velocities[3 * i + 1] = velocities[3 * j + 1] = 0;
    velocities[3 * i + 2] = velocities[3 * j + 2] = 0;
    
    // Make three new particles
    i     = numatoms;
    j     = numatoms + 1;
    int k = numatoms + 2;
    numatoms += 3;
    
    coords[3 * i + 0] = x[0] + x0; 
    coords[3 * j + 0] = x[1] + x0; 
    coords[3 * k + 0] = x[2] + x0;  
    
    coords[3 * i + 1] = y[0] + y0; 
    coords[3 * j + 1] = y[1] + y0; 
    coords[3 * k + 1] = y[2] + y0;
    
    coords[3 * i + 2] = z[0] + z0; 
    coords[3 * j + 2] = z[1] + z0; 
    coords[3 * k + 2] = z[2] + z0;
    
    
    // Set to helium type and mass
    charge[i] = charge[j] = charge[k] = 2;
    type[i] = type[j] = type[k] = 0;
    SetMass(i, (4 * amu)); SetMass(j, (4 * amu)); SetMass(k, (4 * amu));
    
    // Energy after change; also updates forces
    double E_after = UpdateEnergyForces();  
    double ke_after = 0.5 * (4 * amu) * (vx0 * vx0 + vy0 * vy0 + vz0 * vz0) * 3;
    
    double delta_E = ((8.7e6) * ec) - (E_after + ke_after - E_before - ke_before);
    double v_needed = sqrt((2.0 * delta_E / 3.0) / (4 * amu));
    
    // Velocities to make up energy difference
    double vx[3], vy[3], vz[3];    
    for (idx = 0; idx < 3; idx++)
    {
        double norm = 1.0 / sqrt(x[idx] * x[idx] + y[idx] * y[idx] + z[idx] * z[idx]);
        vx[idx] = x[idx] * norm * v_needed;
        vy[idx] = y[idx] * norm * v_needed;
        vz[idx] = z[idx] * norm * v_needed;
    }
    
    // Set velocities
    velocities[3 * i + 0] = F2 * vx[0] + vx0; 
    velocities[3 * j + 0] = F1 * vx[1] + vx0; 
    velocities[3 * k + 0] = F1 * vx[2] + vx0;
    
    velocities[3 * i + 1] = F2 * vy[0] + vy0; 
    velocities[3 * j + 1] = F1 * vy[1] + vy0; 
    velocities[3 * k + 1] = F1 * vy[2] + vy0;
    
    velocities[3 * i + 2] = vz[0] + vz0; 
    velocities[3 * j + 2] = vz[1] + vz0; 
    velocities[3 * k + 2] = vz[2] + vz0;
    
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//          This creates a coordinate array for the α particles             //

//////////////////////////////////////////////////////////////////////////////


void RandomTriton(double *x, double *y, double *z, double R)
{
    // Make a triangle of α particles
    double xc[3], yc[3], zc[3];
    int i;
    
    
    for (i = 0; i < 3; i++)
    {
        xc[i] = R * cos(theta) * (i - 1);
        yc[i] = R * sin(theta) * (i - 1);
        zc[i] = 0;
    }
    
    
    // Make a random quaternion
    double q0, q1, q2, q3;
    q0 = rand_uni();
    q1 = rand_uni();
    q2 = rand_uni();
    q3 = rand_uni();
    NormalizeQ(&q0, &q1, &q2, &q3);
    
    // Rotate triangle by random quaternion
    for (i = 0; i < 3; i++)
        RotatePoint(xc[i], yc[i], zc[i], &(x[i]), &(y[i]), &(z[i]), q0, q1, q2, q3);
    
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//                        I do not know what this is                        //

//////////////////////////////////////////////////////////////////////////////


double NormalizeQ(double *q0, double *q1, double *q2, double *q3)
{
    double N = 1 / norm(*q0, *q1, *q2, *q3);
    *q0 *= N;
    *q1 *= N;
    *q2 *= N;
    *q3 *= N;
    return N;
}

double norm(double q0, double q1, double q2, double q3)
{
    return sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
}

void RotatePoint(double x, double y, double z, double *new_x, double *new_y, double *new_z, double q0, double q1, double q2, double q3)
{
    double t2, t3, t4, t5, t6, t7, t8, t9, t10;
    t2 =   q0 * q1;
    t3 =   q0 * q2;
    t4 =   q0 * q3;
    t5 =  -q1 * q1;
    t6 =   q1 * q2;
    t7 =   q1 * q3;
    t8 =  -q2 * q2;
    t9 =   q2 * q3;
    t10 = -q3 * q3;
    
    double nx, ny, nz;
    nx = 2 * ((t8 + t10) * x + (t6 -  t4) * y + (t3 + t7) * z) + x;
    ny = 2 * ((t4 + t6) *  x + (t5 + t10) * y + (t9 - t2) * z) + y;
    nz = 2 * ((t7 - t3) *  x + (t2 +  t9) * y + (t5 + t8) * z) + z;
    *new_x = nx; *new_y = ny; *new_z = nz;
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//        This puts the masses of each α particle in the mass array         //

//////////////////////////////////////////////////////////////////////////////


void SetMass(int i, double m)
{
    mass[3 * i    ] = m;
    mass[3 * i + 1] = m;
    mass[3 * i + 2] = m;
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//                  This calculates the kinetic energy                      //

//////////////////////////////////////////////////////////////////////////////


double KineticEnergy()
{
    double ke = 0;
    int i;
    for (i = 0; i < 3 * numatoms; i++)
        ke += 0.5 * mass[i] * velocities[i] * velocities[i];
    return ke; 
}


double AverageKineticEnergy()
{
    double ke = KineticEnergy();
    double avg = ke / numatoms;
    return avg;
    
}

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

//      The potential energy of the system is calculated using              //      
//      Velocity Verlet algorithm                                           //

//////////////////////////////////////////////////////////////////////////////


double IterateDynamics(double dt)
{
    
    // Velocity Verlet algorithm
    int i;
    for (i = 0; i < 3 * numatoms; i++)
    {
        // s = vt + 0.5at^2
        coords[i] += velocities[i] * dt + 0.5 * (forces[i] / mass[i]) * dt * dt;
    }
    
    // Copy the force array  
    for (i = 0; i < 3 * numatoms; i++)
    {
        old_forces[i] = forces[i];      
    }
    // Update energy with these new coordinates  
    double energy = UpdateEnergyForces();
    
    // Update the velocity array  
    for (i = 0; i < 3 * numatoms; i++)
    {
        velocities[i] += 0.5 * ((old_forces[i] + forces[i]) * dt) / mass[i];
    }
    
    
    if (calculation_type == MINIMIZE)  
        for (i = 0; i < 3 * numatoms; i++)
            velocities[i] *= 0.999;
    
    return energy;
}

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////

//     This section contains terms to calculate the Coulomb repulsion       //
//     and the Coulomb attraction to the uniform electron gas               //

//////////////////////////////////////////////////////////////////////////////


double UpdateEnergyForces()
{
    double energy = 0;
    int i, j;
    
    // Create a force array for the particles with initial values of 0
    for (i = 0; i < 3 * numatoms; i++)
        forces[i] = 0;
    
    // Coulomb repulsion between ions
    double dx, dy, dz, r2, r;
    for (i = 0; i < numatoms; i++)
    {
        for (j = 0; j < i; j++)
        {
            // Calculate the distance between the particles
            dx = coords[3 * i + 0] - coords[3 * j + 0];
            dy = coords[3 * i + 1] - coords[3 * j + 1];
            dz = coords[3 * i + 2] - coords[3 * j + 2];
            r2 = dx * dx + dy * dy + dz * dz;
            r  = sqrt(r2);
            
            // Store this distance of closest approach 
            capp_distance[i][j] = r;
            capp_time[i][j] = 0;
            
            // The energy due to the Coulombic interactions
            double en = eps0 * charge[i] * charge[j] / r;
            
            // Define the forces
            double fx = en * dx / r2;
            double fy = en * dy / r2;
            double fz = en * dz / r2;
            
            energy += en;
            
            forces[3 * i    ] += fx; 
            forces[3 * j    ] -= fx;
            forces[3 * i + 1] += fy; 
            forces[3 * j + 1] -= fy;
            forces[3 * i + 2] += fz; 
            forces[3 * j + 2] -= fz;
            
            /*//Print the distances of closest approach less than a certain value
             if (r < (5e-8)) 
             {
             printf("%i %i %.3e\n", i, j, r);
             } */
        }
    }
    
    // Coulomb attraction to uniform electron gas
    double total_q = 0;
    for (i = 0; i < numatoms; i++)
    {    
        // The charge of the electron gas is the negative of the total charge of the borons and hydrogens
        total_q -= charge[i];
    }
    
    for (i = 0; i < numatoms; i++)
    {
        if ((type[i] == 1) || (type[i] == 2)) 
        {
            // Calculate the position vector of each particle
            double dx = coords[3 * i    ];
            double dy = coords[3 * i + 1];
            double dz = coords[3 * i + 2];
            double r2 = dx * dx + dy * dy + dz * dz;
            double r  = sqrt(r2);
        }
        
        // Calculate the attraction to the cloud of electrons if it is within the cloud
        double en, fx, fy, fz;
        
        if (r < radius_charge)
        {
            double factor = charge[i] * eps0 * 0.5 * total_q / radius_charge;
            en = factor * (3 - r * r  / (radius_charge * radius_charge));
            fx = factor * (    2 * dx / (radius_charge * radius_charge));
            fy = factor * (    2 * dy / (radius_charge * radius_charge));
            fz = factor * (    2 * dz / (radius_charge * radius_charge));
        }
        
        else
        {
            en = charge[i] * eps0 * total_q / r;
            fx = en * dx / r2;
            fy = en * dy / r2;
            fz = en * dz / r2;
        }
        
        energy += en;
        forces[3 * i    ] += fx;
        forces[3 * i + 1] += fy;
        forces[3 * i + 2] += fz;
    }
    
    return energy;
}

//////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////

//      The section contains the functions calculating the temperature      //      
//                      and pressure of the system                          //

//////////////////////////////////////////////////////////////////////////////



double Temperature()
{
    double ke = KineticEnergy();
    double avg = AverageKineticEnergy();
    double temp = 0;
    temp += ((2.0 / 3.0) * avg) / (numatoms * KB);
    return temp;
}


double DotProduct()
{
    int i;   
    double rf = 0;
    for (i = 0; i < 3 * numatoms; i++)
        rf += (coords[i] * forces[i]);
    return rf;
}  


double Pressure(double rf)
{
    double pres = 0;
    double temp = Temperature();
    double ke = KineticEnergy();
    pres += (((numatoms * KB * temp) / ((4.0 / 3.0) * pi * pow(radius_charge, 3))) + (rf / (4.0 * pi * pow(radius_charge, 3))));
    return pres;
}

//////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////                                                                  //////
//////            This section is related to the α particles            //////
//////                                                                  //////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////







void CompareAnalyticNumericalGradients()
{
    int i;
    for (i = 0; i < 3 * numatoms; i++)
    {
        double delta = 0.0001;
        coords[i] += delta;
        double E1 = UpdateEnergyForces();
        coords[i] -= 2 * delta;
        double E2 = UpdateEnergyForces();
        coords[i] += delta;
        
        printf("%e %e\n", -(E1 - E2) / (2 * delta), forces[i]);
    }
}



//////////////////////////////////////////////////////////////////////////////

//                     This outputs the final coordinates                   //

//////////////////////////////////////////////////////////////////////////////


void WriteCoords(FILE *fp, double *coords)
{
    fprintf(fp, "%i\n\n", max_numatoms);
    int i;
    for (i = 0; i < max_numatoms; i++)
    {
        if (type[i] == 1)
            fprintf(fp, "H %f %f %f\n", coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
        else if (type[i] == 2)
            fprintf(fp, "B %f %f %f\n", coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
        else if (type[i] == 0)
            fprintf(fp, "He %f %f %f\n", coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]);
        else
            fprintf(fp, "He 10000 0 0\n"); 
    }
}

//////////////////////////////////////////////////////////////////////////////








