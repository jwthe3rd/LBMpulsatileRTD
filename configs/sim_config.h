#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using olb::util::Randomizer;
using olb::util::Timer;
using namespace particles;
using namespace particles::subgrid;
using namespace particles::communication;
using namespace particles::dynamics;
using namespace particles::creators;
using namespace particles::io;

using T = FLOATING_POINT_TYPE;
typedef D3Q19<> DESCRIPTOR;
typedef SubgridParticle3D PARTICLETYPE;
typedef SmagorinskyBGKdynamics<T,DESCRIPTOR> BulkDynamics;

#define BOUZIDI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Physical Settings
const T tau = 0.51;
const T rhoP = 1045.0; // Physical Density (kg/m3)
const T muP = 0.00287;
const T nuP = muP / rhoP; // Physical Nu
const T nuL = (tau - nuP) / 3; // Lattice Nu (LU)
const T conversionFactor = 0.001;
const T smagCoeff = 0.12;


// Sim Resolution Determination
const int N = 25;
const T charMinL =0.001; // m
const T deltaX = charMinL / N; // LU
const T deltaT = (nuL * deltaX * deltaX) / nuP; // s
const T maxAllowableLatticeVel = 5.;
char stlFileName[] = "geo/p89.stl";

// Sim Time Settings
const T fluidMaxPhysT = T( 2 );     // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 2 ); // max. particle simulation time in s, SI unit
const T fluidPeriod = 0.25;

// Average Velocity Determination
const T flowRate = 0.5 * 1.0e-6; // m3/s
const T radInlet = 0.0014;//0.00161; // m
const T avgRad = (2.*radInlet+charMinL) / 4.; // for calculating the inlet velocity
const T avgVel = flowRate / (M_PI * avgRad* avgRad); // m/s
const T avgLVel = (avgVel*deltaT) / deltaX; // LU
const T radOutlet = 0.001;//0.00075; // m
const T distForCalcInFlowOutFlow = 2.; //LU
bool sinInletVelocityBoundaryCondition = false;
bool javadInletVelocityBoundaryCondition = true;

// center of inflow and outflow regions [*units based on conversionFactor (i.e. 1.=m, 0.001=mm, 1000.=km, etc.)]
Vector<T, 3> outletCenter( 0.00454033663, 0.0086697265, 0.00470741802 );
Vector<T, 3> inletCenter( -0.01931674079, 0.02270769342, -0.02595900419 );


// normals of inflow and outflow regions
Vector<T, 3> outletNormal(0.1834349009, -0.4027303427, 0.8967496352 );
Vector<T, 3> inletNormal( -0.9670365813, 0.04464441815, -0.2506932916);

// Particle Settings
std::size_t noOfParticles = 500;   // total number of inserted particles
const T radius = 7.5e-5;            // particles radius
const T partRho = 1045.0;
const T charLIn = 0.0002;
Vector<T, 3> particleInjectionP(inletCenter[0]-(charLIn*inletNormal[0]), inletCenter[1]-(charLIn*inletNormal[1]), inletCenter[2]-(charLIn*inletNormal[2]));
const T injectionRadius = 5.0; // LU away from the wall
const T injectionLength = 0.5; // LU, length of the injection cylinder (must be cylinder because 3D, small to make it like a surface)
const T injectionPeriod = 0.25;

//Set capture method:
// materialCapture: based on material number
// wallCapture:     based on more accurate stl description
typedef enum {materialCapture, wallCapture} ParticleDynamicsSetup;
const ParticleDynamicsSetup particleDynamicsSetup = materialCapture;


//writing data constants
char vtkFileName[] = "rtdVal";
char vtkParticleFileName[] = "particle";
const T physVTKiter = 0.01;
const T physStatiter = 0.005;
