/** \file proto.h

Declaration of all the functions of the FARGO code

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

void masterprint (const char *template, ...);
void mastererr (const char *template, ...);
real GetGlobalIFrac ();
void prs_exit ();
void prs_error();
void message ();
PolarGrid    *CreatePolarGrid();
void MultiplyPolarGridbyConstant ();
void DumpSources ();
void UpdateLog ();
void ReadfromFile ();
void InitLabel ();
void Initialization ();
void var();
void ReadVariables();
void PrintUsage ();
real TellNbOrbits ();
real TellNbOutputs ();
void TellEverything ();
void GiveTimeInfo ();
void InitSpecificTime ();
void GiveSpecificTime ();
void EmptyPlanetSystemFile ();
void WritePlanetFile ();
void WritePlanetSystemFile ();
void WriteBigPlanetFile ();
void WriteBigPlanetSystemFile ();
real GetfromPlanetFile ();
void RestartPlanetarySystem ();
void WriteDiskPolar();
void WriteDim ();
void SendOutput ();
void FillForcesArrays ();
void AdvanceSystemFromDisk ();
real ConstructSequence ();
void InitGas ();
void AccreteOntoPlanets ();
void FindOrbitalElements ();
int FindNumberOfPlanets ();
PlanetarySystem *AllocPlanetSystem ();
void FreePlanetary ();
PlanetarySystem *InitPlanetarySystem ();
void ListPlanets ();
real GetPsysInfo ();
void RotatePsys ();
real GasTotalMass ();
real GasMomentum ();
void DivisePolarGrid ();
void InitComputeAccel ();
void OpenBoundary ();
void NonReflectingBoundary ();
void ApplyOuterSourceMass ();
void ApplyBoundaryCondition ();
void CorrectVtheta ();
boolean DetectCrash ();
void FillPolar1DArrays ();
void InitEuler ();
real min2 ();
real max2 ();
void ActualiseGas ();
void AlgoGas ();
void SubStep1 ();
void SubStep2 ();
int ConditionCFL ();
real Sigma();
void FillSigma();
void RefillSigma ();
void Transport ();
void OneWindRad ();
void ComputeThetaElongations ();
void ComputeAverageThetaVelocities ();
void ComputeResiduals ();
void ComputeConstantResidual ();
void AdvectSHIFT ();
void OneWindTheta ();
void QuantitiesAdvection ();
void ComputeExtQty ();
void ComputeSpeQty ();
void InitTransport () ;
void ComputeStarRad ();
void ComputeStarTheta ();
void ComputeLRMomenta ();
void ComputeVelocities ();
real VanLeerRadial ();
void VanLeerTheta ();
void InitViscosity ();
void ViscousTerms ();
void AllocateComm ();
void CommunicateBoundaries ();
void handfpe();
void setfpe ();
void merge ();
void ReadPrevDim ();
void CheckRebin ();
void SplitDomain ();
void InitVariables();
real FViscosity ();
real AspectRatio ();
void MakeDir ();
FILE *fopenp ();
/* #THORIN: energy equation implementation */
void SubStep3 ();
void ComputeSoundSpeed ();
void ComputeTemperatureField ();
void ComputePressureField ();
real ThicknessSmoothing ();
void mpi_make1Dprofile ();
void InitGasDensityEnergy ();
real GasTotalEnergy ();
real Energy ();
void FillEnergy ();
void RefillEnergy ();
void FillVtheta ();
void InitGasVelocity ();
real InitCoolingTime ();
void FillCoolingTime ();
real InitQplus ();
void FillQplus ();
void UpdateDivVelocAndStressTensor ();
void UpdateVelocityWithViscousTerms ();
void ImposeKeplerianEdges ();
void ReadfromAsciiFile ();
/* #THORIN: radiative diffusion */
void InitRadiatDiffusionFields ();
void CalculateQminus ();
void CalculateFlaring ();
void CalculateQirr ();
void ImplicitRadiativeDiffusion ();
void TemperatureGradient ();
void MidplaneVolumeDensity ();
void OpacityProfile ();
real FluxLimiterValue ();
real EffectiveOpticalDepth ();
void IterateRelaxationParameter ();
int SuccessiveOverrelaxation ();
void DiffusionCoefs ();
void SynchronizeOverlapFields ();
void ChessBoardIndexing ();
void SetWaveKillingZones ();
void DampingBoundary ();
void ActualizeQbalance ();
/* #THORIN: FARGO-REBOUND interface */
struct reb_simulation *SetupReboundSimulation ();
void SetupIntegratorParams ();
void AdvanceSystemRebound ();
void AdditionalForces ();
void OutputElements ();
void OutputNbodySimulation ();
boolean ChkCloseEncWithPl ();
void DiscardParticlesDist ();
void DiscardParticlesUnbound ();
int ResolveCollisions ();
struct reb_simulation *RestartReboundSimulation ();
void SynchronizeFargoRebound ();
void MinStepForRebound ();
real DampingTW04 ();
real GetPsysInfoFromRsim ();
void DumpOmegaFrame ();
real GetOmegaFrame ();
/* #THORIN: pebble accretion */
void InitPebbleArrays ();
void EquilPebbleDisk ();
void InitPebblesViaFlux ();
void RestartPebbleDisk ();
void PebbleStokesNumbers ();
void AccretePebblesOntoPlanets ();
void CorrectPebblesVtheta ();
void EvolvePebbleDisk ();
void WritePebbles ();
real Trapzd ();
real IntegrateColumnMass ();
void EtaPressureSupport ();
void DampPebbles ();
void TransportPebbles ();
void OneWindRadPebbles ();
void OneWindThetaPebbles ();
void QuantitiesAdvectionPebbles ();
void SourceTermsPebbles ();
void SubStep1Pebbles ();
boolean DetectCrashPebbles ();
void SynchronizePebbleDisc ();
void CriticalCharTime ();
void ParticleDiffusion ();
void BckpFieldsForBC ();
void ParametricAccretion ();
/* #THORIN: tools */
void CreateTorqueMapInfile ();
void DiffusionTimestep ();

