#include "../configs/sim_config.h"

T getVelocityPulse(UnitConverter<T,DESCRIPTOR> const& converter,
                    int iT)
{

    T vel_char;

    T t = (converter.getPhysTime(iT)*1.111-floor(converter.getPhysTime(iT)*1.111))/1.111;


    if (t <= 0.15)
    {
      vel_char =1.322515215*0.722378598771807*((0.349999999999992+9.39592352092547*t-191.515151515187*(pow(t,2))+871.933621933773*(pow(t,3))));
    }
    else
    { 
      vel_char =1.322515215*0.722378598771807*((-8.87872100015116+158.853244453899*t-1009.813663592*(pow(t,2))+3083.59699791101*(pow(t,3))-4840.042878547*(pow(t,4))+3765.53701950272*(pow(t,5))-1142.95591854388*(pow(t,6))-8.95498219229757*(pow(t,7))));
    }

    return vel_char;
}

void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader,
                      SuperGeometry<T,3>& superGeometry)
{
    OstreamManager clout( std::cout,"prepareGeometry" );
    clout << "Prepare Geometry ..." << std::endl;

    //take all values in the 3d volume of the cylinder.stl and assign them
    //material number 2
    superGeometry.rename(0, 2, indicator);
    //Now rename all values on the surface of the stl to material 1
    superGeometry.rename(2, 1, stlReader);
    //cleans auxillary voxels
    superGeometry.clean();

    //setup the inlet and outlet faces

    IndicatorCircle3D<T> inflow(inletCenter, inletNormal, radInlet);
    IndicatorCylinder3D<T> layerInflow(inflow, 2. * converter.getConversionFactorLength() );
    superGeometry.rename(2,3,1,layerInflow);

    IndicatorCircle3D<T> outflow(outletCenter, outletNormal, radOutlet);
    IndicatorCylinder3D<T> layerOutflow(outflow, 2. * converter.getConversionFactorLength() );
    superGeometry.rename(2, 4, 1, layerOutflow);

    superGeometry.clean();
    superGeometry.innerClean( 3 );
    superGeometry.checkForErrors();

    superGeometry.print();
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )


{
    OstreamManager clout( std::cout,"prepareLattice" );
    clout << "Prepare Lattice ..." << std::endl;

    const T omega = converter.getLatticeRelaxationFrequency();

    //set material 1 to bulk dynamics
    lattice.defineDynamics<BulkDynamics>(superGeometry, 1);

    setBounceBackBoundary(lattice, superGeometry, 2);
    setInterpolatedVelocityBoundary<T, DESCRIPTOR>(lattice, omega, superGeometry, 3);
    lattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator(4));
    setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(4));

    // Initial conditions
    AnalyticalConst3D<T, T> rhoF( 1 );
    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF( velocity );

    lattice.defineRhoU( superGeometry.getMaterialIndicator({1, 3, 4}),rhoF,uF );
    lattice.iniEquilibrium( superGeometry.getMaterialIndicator({1, 3, 4}),rhoF,uF );
    lattice.setParameter<descriptors::OMEGA>(omega);
    lattice.setParameter<collision::LES::Smagorinsky>(smagCoeff);

    lattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;


}

void prepareParticles( UnitConverter<T,DESCRIPTOR> const& converter,
  SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
  SolidBoundary<T,3>& wall,
  SuperIndicatorMaterial<T,3>& materialIndicator,
  ParticleDynamicsSetup particleDynamicsSetup,
  Randomizer<T>& randomizer)
{
  OstreamManager clout( std::cout, "prepareParticles" );
  clout << "Prepare Particles ..." << std::endl;

  //Add selected particle dynamics
  if (particleDynamicsSetup==wallCapture){
    //Create verlet dynamics with material aware wall capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialAwareWallCapture<T,PARTICLETYPE>>(
        wall, materialIndicator);
  } else {
    //Create verlet dynamics with material capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialCapture<T,PARTICLETYPE>>(materialIndicator);
  }

  // particles generation at inlet3
  //Vector<T, 3> c( inletCenter );
 //c[0] = particleInjectionX;
  IndicatorCircle3D<T> inflowCircle( particleInjectionP, inletNormal, radInlet -
                                     converter.getConversionFactorLength() * injectionRadius);
  IndicatorCylinder3D<T> inletCylinder( inflowCircle, injectionLength*
                                        converter.getConversionFactorLength() );
                  

  //Add particles
  addParticles( superParticleSystem, inletCylinder, partRho, radius, noOfParticles, randomizer );

  //Print super particle system summary
  superParticleSystem.print();
  clout << "Prepare Particles ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{
    OstreamManager clout( std::cout,"setBoundaryValues" );

    // No of time steps for smooth start-up
    const int iTmaxStart = converter.getLatticeTime( fluidMaxPhysT*0.5 );
    const int iTupdate = converter.getLatticeTime( 0.001 );
    const int iTperiod = converter.getLatticeTime(fluidPeriod);
    //Vector <T, 3> maxVelocity(T(), T(), T());

    if ( iT%iTupdate==0 && iT<=iTmaxStart) {//if ( iT%iTupdate==0) {
        sLattice.setProcessingContext(ProcessingContext::Evaluation);

        // Smooth start curve, sinus
        // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

        // Smooth start curve, polynomial

        // if (sinInletVelocityBoundaryCondition)
        // {
        //   SinusStartScale<T,int> nSinusStartScale( iTperiod,2.*converter.getLatticeVelocity(avgVel) );

        //   // Creates and sets the Poiseuille inflow profile using functors
        //   int iTvec[1]= {iT};
        //   //T maxVelocity[1]= {T()};
        //   nSinusStartScale( maxVelocity[0],maxVelocity[1],iTvec );
        //   CirclePoiseuille3D<T> computedVelBC( superGeometry,3,maxVelocity[0], T() );
        // }
        // else if (javadInletVelocityBoundaryCondition)
        // {
        //  PolynomialStartScale<T, int> startScale( iTmaxStart, T(1));
        //  int iTvec[1] = {iT};
        //  T frac[1]={};
        //  startScale(frac, iTvec);
        Vector<T, 3> maxVelocity(2.*converter.getLatticeVelocity(getVelocityPulse(converter, iT)), T(), T());
         

        CirclePoiseuille3D<T> computedVelBC(superGeometry, 3, maxVelocity[0], T());
        //}
        // else
        // {
		// clout << "poopy" << std::endl;
		// std::exit(0);
    //         PolynomialStartScale<T,int> startScale( iTmaxStart, T( 1 ) );
    //       // Creates and sets the Poiseuille inflow profile using functors
    //         int iTvec[1]= {iT};
    //         T frac[1]= {};
    //         startScale( frac,iTvec );
    //         std::vector<T> maxVelocity( 3,0 );
    //         maxVelocity[0] = 2.*frac[0]*converter.getLatticeVelocity(getVelocityPulse(converter, iT));

    //         T distance2Wall = converter.getConversionFactorLength()/2.;
    //         CirclePoiseuille3D<T> computedVelBC( superGeometry, 3, maxVelocity[0], distance2Wall );
    //       }
        // }

        sLattice.defineU( superGeometry, 3, computedVelBC );
        sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);

        if ( iT % (10*iTupdate) == 0 && iT <= iTmaxStart ) {
        clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
        }

        sLattice.setProcessingContext(ProcessingContext::Simulation);
    }

}

bool getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR>& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, STLreader<T>& stlReader,
                 SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem)
{
  OstreamManager clout( std::cout, "getResults" );
  const int vtkIter  = converter.getLatticeTime( physVTKiter );
  const int statIter = converter.getLatticeTime( physStatiter );

  SuperVTMwriter3D<T> vtmWriter(vtkFileName);
  VTUwriter<T,PARTICLETYPE,true> superParticleWriter(vtkParticleFileName,false,false);

  if ( iT==0 ) {


    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
    superParticleWriter.createMasterFile();


  }

  // Writes output to vtk files for viewing in Paraview
  if ( iT%vtkIter==0 && iT <= converter.getLatticeTime(fluidMaxPhysT) ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter(vtkFileName);
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      SuperLatticeGeometry3D<T, DESCRIPTOR> materials(sLattice, superGeometry);

      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      vtmWriter.addFunctor(materials);

      task(vtmWriter, iT);

    });
  }
  // Writes output on the console
  if ( iT%statIter==0 && iT <= converter.getLatticeTime(fluidMaxPhysT) ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    clout << sLattice.getStatistics().getMaxU() << std::endl;
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Flux at the inflow and outflow region
    std::vector<int> materials = { 1, 3, 4};

    //Surface Areas for computing inflow and outflow
    Vector<T, 3> outletCircle(outletCenter[0]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*outletNormal[0]), outletCenter[1]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*outletNormal[1]), outletCenter[2]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*outletNormal[2]));
    Vector<T, 3> inletCircle(inletCenter[0]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*inletNormal[0]), inletCenter[1]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*inletNormal[1]), inletCenter[2]-(distForCalcInFlowOutFlow*converter.getConversionFactorLength()*inletNormal[2]));

    IndicatorCircle3D<T> inflow(inletCircle, inletNormal, radInlet+2.*converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Analytical);
    vFluxInflow.print( "inflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Analytical );
    pFluxInflow.print( "inflow","N","mmHg" );

    IndicatorCircle3D<T> outflow(outletCircle,  outletNormal, radOutlet+2.*converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow( sLattice, converter, superGeometry, outflow, materials, BlockDataReductionMode::Analytical );
    vFluxOutflow.print( "outflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow( sLattice, converter, superGeometry, outflow, materials, BlockDataReductionMode::Analytical);
    pFluxOutflow.print( "outflow","N","mmHg" );


    if ( sLattice.getStatistics().getMaxU() > maxAllowableLatticeVel ) {
      clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
      std::exit(0);
    }
  }

  if (iT%vtkIter==0 && iT >= converter.getLatticeTime(injectionPeriod)){
    SuperParticleGroupedFieldF<T,PARTICLETYPE,GENERAL,POSITION> particlePosition( superParticleSystem );
    SuperParticleGroupedFieldF<T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE> particleActivity( superParticleSystem );

    superParticleWriter.addFunctor( particlePosition );
    superParticleWriter.addFunctor( particleActivity );

    superParticleWriter.write(iT - converter.getLatticeTime(fluidMaxPhysT));
  }
  if (iT%statIter==0 && iT >= converter.getLatticeTime(injectionPeriod)){

    //Purge invalid particles (delete invalidated particles)
    purgeInvalidParticles<T,PARTICLETYPE>( superParticleSystem );

    //Define materials for capture rate
    std::vector<int> materialsOutput {4};
    SuperIndicatorMaterial<T,3> materialIndicatorOutput (superGeometry, materialsOutput);

    //Perform capture statistics
    std::size_t noActive;
    captureStatistics( superParticleSystem, materialIndicatorOutput, noActive );

    // true as long as certain amount of active particles
    if ( noActive < 0.0001 * noOfParticles
         && iT > 0.9*converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT )) {
      return false;
    }
    //Additional criterion added 02.02.23
    if ( noActive==0) {
      return false;
    }
  }
    return true;
}
