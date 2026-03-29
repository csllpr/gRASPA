
#include "axpy.h"

#include "mc_single_particle.h"
#include "mc_swap_moves.h"
#include "mc_box.h"



#include "write_data.h"

#include "print_statistics.cuh"

//#include "lambda.h"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <optional>

#include <fstream>

//#include <format>

inline void Copy_AtomData_from_Device(Atoms* System, Atoms* d_a, Components& SystemComponents, Boxsize& HostBox, Simulations& Sims)
{
  cudaMemcpy(System, d_a, SystemComponents.NComponents.x * sizeof(Atoms), cudaMemcpyDeviceToHost);
  for(size_t ijk=0; ijk < SystemComponents.NComponents.x; ijk++)
  {
    if(SystemComponents.HostSystem[ijk].Allocate_size != System[ijk].Allocate_size)
    {
      // if the host allocate_size is different from the device, allocate more space on the host
      SystemComponents.HostSystem[ijk].pos       = (double3*) malloc(System[ijk].Allocate_size*sizeof(double3));
      SystemComponents.HostSystem[ijk].scale     = (double*)  malloc(System[ijk].Allocate_size*sizeof(double));
      SystemComponents.HostSystem[ijk].charge    = (double*)  malloc(System[ijk].Allocate_size*sizeof(double));
      SystemComponents.HostSystem[ijk].scaleCoul = (double*)  malloc(System[ijk].Allocate_size*sizeof(double));
      SystemComponents.HostSystem[ijk].Type      = (size_t*)  malloc(System[ijk].Allocate_size*sizeof(size_t));
      SystemComponents.HostSystem[ijk].MolID     = (size_t*)  malloc(System[ijk].Allocate_size*sizeof(size_t));
      SystemComponents.HostSystem[ijk].Allocate_size = System[ijk].Allocate_size;
    }
  
    cudaMemcpy(SystemComponents.HostSystem[ijk].pos, System[ijk].pos, sizeof(double3)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(SystemComponents.HostSystem[ijk].scale, System[ijk].scale, sizeof(double)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(SystemComponents.HostSystem[ijk].charge, System[ijk].charge, sizeof(double)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(SystemComponents.HostSystem[ijk].scaleCoul, System[ijk].scaleCoul, sizeof(double)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(SystemComponents.HostSystem[ijk].Type, System[ijk].Type, sizeof(size_t)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(SystemComponents.HostSystem[ijk].MolID, System[ijk].MolID, sizeof(size_t)*System[ijk].Allocate_size, cudaMemcpyDeviceToHost);
    SystemComponents.HostSystem[ijk].size = System[ijk].size;
  }
  HostBox.Cell = (double*) malloc(9 * sizeof(double));
  HostBox.InverseCell = (double*) malloc(9 * sizeof(double));
  cudaMemcpy(HostBox.Cell,        Sims.Box.Cell,        sizeof(double)*9, cudaMemcpyDeviceToHost);
  cudaMemcpy(HostBox.InverseCell, Sims.Box.InverseCell, sizeof(double)*9, cudaMemcpyDeviceToHost);
  HostBox.Cubic = Sims.Box.Cubic;
}

inline void GenerateRestartMovies(Variables& Vars, size_t systemId, PseudoAtomDefinitions& PseudoAtom, int SimulationMode)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations& Sims = Vars.Sims[systemId];
  Boxsize& HostBox = Vars.Box[systemId];
  //Generate Restart file during the simulation, regardless of the phase
  Atoms device_System[SystemComponents.NComponents.x];
  Copy_AtomData_from_Device(device_System, Sims.d_a, SystemComponents, HostBox, Sims);
  create_Restart_file(0, SystemComponents.HostSystem, SystemComponents, SystemComponents.FF, HostBox, PseudoAtom.Name, systemId);
  Write_All_Adsorbate_data(0, SystemComponents.HostSystem, SystemComponents, SystemComponents.FF, HostBox, PseudoAtom.Name, systemId);
  //Only generate LAMMPS data movie for production phase
  if(SimulationMode == PRODUCTION)  create_movie_file(SystemComponents.HostSystem, SystemComponents, HostBox, PseudoAtom.Name, systemId);
}

static inline size_t CeilDivide(size_t numerator, size_t denominator)
{
  if(denominator == 0) throw std::runtime_error("Attempted division by zero while preparing adaptive production");
  return numerator / denominator + ((numerator % denominator) > 0 ? 1 : 0);
}

static inline double FrameworkDensityFromVolume(const Components& SystemComponents, double volume, const Units& Constants)
{
  if(volume <= 0.0) return std::numeric_limits<double>::quiet_NaN();

  double CellMass = 0.0;
  size_t NCell = static_cast<size_t>(SystemComponents.NumberofUnitCells.x) *
                 static_cast<size_t>(SystemComponents.NumberofUnitCells.y) *
                 static_cast<size_t>(SystemComponents.NumberofUnitCells.z);
  for(size_t j = 0; j < static_cast<size_t>(SystemComponents.NComponents.y); j++)
    CellMass += SystemComponents.MolecularWeight[j] * NCell;

  return CellMass * 1.0e-3 / (Constants.Avogadro * volume * 1.0e-30);
}

static inline void ResetAdaptiveProductionState(Components& SystemComponents)
{
  SystemComponents.AdaptiveState = AdaptiveProductionState{};
  SystemComponents.AdaptiveState.HenryStatus.assign(SystemComponents.NComponents.x, AdaptiveObservableStatus{});
  SystemComponents.AdaptiveState.HeatStatus.assign(SystemComponents.NComponents.x, AdaptiveObservableStatus{});
  for(size_t comp = 0; comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
  {
    bool monitorHenry = comp < SystemComponents.AdaptiveTargets.size() && SystemComponents.AdaptiveTargets[comp].MonitorHenryCoefficient;
    bool monitorHeat  = comp < SystemComponents.AdaptiveTargets.size() && SystemComponents.AdaptiveTargets[comp].MonitorHeatOfAdsorption;
    SystemComponents.AdaptiveState.HenryStatus[comp].Enabled = monitorHenry;
    SystemComponents.AdaptiveState.HeatStatus[comp].Enabled = monitorHeat;
  }
}

static inline void InitializeProductionStatisticsStorage(Components& SystemComponents, size_t nblocks)
{
  SystemComponents.Nblock = nblocks;
  SystemComponents.BookKeepEnergy.assign(nblocks, MoveEnergy{});
  SystemComponents.BookKeepEnergy_SQ.assign(nblocks, MoveEnergy{});
  SystemComponents.AverageEnergy = MoveEnergy{};
  SystemComponents.AverageEnergy_Errorbar = MoveEnergy{};

  std::vector<double> zeroDoubleBlocks(nblocks, 0.0);
  SystemComponents.EnergyTimesNumberOfMolecule.assign(SystemComponents.NComponents.x, zeroDoubleBlocks);
  SystemComponents.VolumeAverage.assign(nblocks, {0.0, 0.0});
  SystemComponents.DensityPerComponent.assign(SystemComponents.NComponents.x, std::vector<double2>(nblocks, {0.0, 0.0}));
  if(SystemComponents.AmountOfExcessMolecules.size() > 0)
    SystemComponents.ExcessLoading.assign(SystemComponents.NComponents.x, std::vector<double2>(nblocks, {0.0, 0.0}));
  else
    SystemComponents.ExcessLoading.clear();

  for(size_t comp = 0; comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
  {
    Move_Statistics& MoveStats = SystemComponents.Moves[comp];
    MoveStats.BlockID = 0;
    MoveStats.MolAverage.assign(nblocks, {0.0, 0.0});
    MoveStats.Rosen.assign(nblocks, RosenbluthWeight{});
    MoveStats.MolSQPerComponent.assign(SystemComponents.NComponents.x, std::vector<double>(nblocks, 0.0));
    MoveStats.WidomEnergy = MoveEnergy{};
    MoveStats.WidomEnergy_ERR = MoveEnergy{};
  }
}

static inline bool AdaptiveProductionEnabled(const Variables& Vars)
{
  return Vars.SimulationMode == PRODUCTION && Vars.AdaptiveProduction.Enabled;
}

static inline size_t AdaptiveBlockSize(size_t blockId, size_t totalCycles, size_t batchCycles)
{
  size_t blockStart = blockId * batchCycles;
  if(blockStart >= totalCycles) return 0;
  size_t remainingCycles = totalCycles - blockStart;
  return std::min(batchCycles, remainingCycles);
}

static inline double StudentT95TwoSided(size_t degreesOfFreedom)
{
  static constexpr double t_table[] = {
    0.0,
    12.706, 4.303, 3.182, 2.776, 2.571,
    2.447, 2.365, 2.306, 2.262, 2.228,
    2.201, 2.179, 2.160, 2.145, 2.131,
    2.120, 2.110, 2.101, 2.093, 2.086,
    2.080, 2.074, 2.069, 2.064, 2.060,
    2.056, 2.052, 2.048, 2.045, 2.042
  };
  if(degreesOfFreedom < sizeof(t_table) / sizeof(t_table[0]))
    return t_table[degreesOfFreedom];
  return 1.96;
}

static inline AdaptiveObservableStatus BuildAdaptiveObservableStatus(const std::vector<double>& estimates, const AdaptiveProductionSettings& Settings, double absoluteTolerance)
{
  AdaptiveObservableStatus Status;
  Status.ValidBatches = estimates.size();
  if(estimates.empty()) return Status;

  double sum = 0.0;
  double sqsum = 0.0;
  for(double value : estimates)
  {
    sum += value;
    sqsum += value * value;
  }
  Status.Mean = sum / static_cast<double>(estimates.size());
  if(estimates.size() < 2) return Status;

  double numerator = sqsum - sum * sum / static_cast<double>(estimates.size());
  numerator = std::max(0.0, numerator);
  double sampleVariance = numerator / static_cast<double>(estimates.size() - 1);
  double standardError = std::sqrt(sampleVariance / static_cast<double>(estimates.size()));
  Status.HalfWidth = StudentT95TwoSided(estimates.size() - 1) * standardError;
  Status.RelativeHalfWidth = Status.HalfWidth / std::max(std::abs(Status.Mean), Settings.RelativeFloor);
  Status.Available = std::isfinite(Status.HalfWidth) && std::isfinite(Status.RelativeHalfWidth);
  if(!Status.Available) return Status;

  bool relativeEnabled = Settings.RelativeTolerance > 0.0;
  bool absoluteEnabled = absoluteTolerance >= 0.0;
  bool hasCriterion = relativeEnabled || absoluteEnabled;
  bool relativePass = relativeEnabled && Status.RelativeHalfWidth <= Settings.RelativeTolerance;
  bool absolutePass = absoluteEnabled && Status.HalfWidth <= absoluteTolerance;

  bool criteriaPassed = false;
  if(hasCriterion && Settings.CriteriaMode == ADAPTIVE_CRITERIA_ALL)
  {
    criteriaPassed = true;
    if(relativeEnabled) criteriaPassed = criteriaPassed && relativePass;
    if(absoluteEnabled) criteriaPassed = criteriaPassed && absolutePass;
  }
  else if(hasCriterion)
  {
    criteriaPassed = relativePass || absolutePass;
  }
  Status.Passed = Status.ValidBatches >= Settings.MinimumBatches && criteriaPassed;
  return Status;
}

static inline bool ComputeHenryEstimateForBlock(const Components& SystemComponents, const Units& Constants, const AdaptiveProductionSettings& Settings, size_t component, size_t blockId, size_t totalCycles, size_t batchCycles, double& henry)
{
  size_t cyclesInBlock = AdaptiveBlockSize(blockId, totalCycles, batchCycles);
  if(cyclesInBlock == 0) return false;

  const RosenbluthWeight& Rosen = SystemComponents.Moves[component].Rosen[blockId];
  if(Rosen.Total.z < static_cast<double>(Settings.MinimumWidomSamples)) return false;

  double averageWr = Rosen.Total.x / Rosen.Total.z;
  if(!std::isfinite(averageWr) || averageWr <= 0.0) return false;

  double averageVolume = SystemComponents.VolumeAverage[blockId].x / static_cast<double>(cyclesInBlock);
  double rhoFramework = FrameworkDensityFromVolume(SystemComponents, averageVolume, Constants);
  if(!std::isfinite(rhoFramework) || rhoFramework <= 0.0) return false;

  henry = averageWr / (Constants.gas_constant * SystemComponents.Temperature * rhoFramework);
  return std::isfinite(henry);
}

static inline bool ComputeHeatEstimateForBlock(Components& SystemComponents, const Units& Constants, size_t component, size_t blockId, size_t totalCycles, size_t batchCycles, double& heat)
{
  if(component < static_cast<size_t>(SystemComponents.NComponents.y)) return false;

  size_t cyclesInBlock = AdaptiveBlockSize(blockId, totalCycles, batchCycles);
  if(cyclesInBlock == 0) return false;

  size_t NumberOfAdsorbateComponents = static_cast<size_t>(SystemComponents.NComponents.x - SystemComponents.NComponents.y);
  if(NumberOfAdsorbateComponents == 0) return false;

  double Average_E = SystemComponents.BookKeepEnergy[blockId].total();
  Average_E -= SystemComponents.BookKeepEnergy[blockId].HHVDW;
  Average_E -= SystemComponents.BookKeepEnergy[blockId].HHReal;
  Average_E -= SystemComponents.BookKeepEnergy[blockId].HHEwaldE;
  Average_E /= static_cast<double>(cyclesInBlock);

  std::vector<std::vector<double>> matrix(NumberOfAdsorbateComponents, std::vector<double>(NumberOfAdsorbateComponents, 0.0));
  std::vector<std::vector<double>> tempMatrix(NumberOfAdsorbateComponents, std::vector<double>(NumberOfAdsorbateComponents, 0.0));

  for(size_t compi = static_cast<size_t>(SystemComponents.NComponents.y); compi < static_cast<size_t>(SystemComponents.NComponents.x); compi++)
  {
    for(size_t compj = static_cast<size_t>(SystemComponents.NComponents.y); compj < static_cast<size_t>(SystemComponents.NComponents.x); compj++)
    {
      size_t adjustCompi = compi - static_cast<size_t>(SystemComponents.NComponents.y);
      size_t adjustCompj = compj - static_cast<size_t>(SystemComponents.NComponents.y);
      double Average_N = SystemComponents.Moves[compi].MolAverage[blockId].x / static_cast<double>(cyclesInBlock);
      double Average_Nj = SystemComponents.Moves[compj].MolAverage[blockId].x / static_cast<double>(cyclesInBlock);
      double Average_NxNj = SystemComponents.Moves[compi].MolSQPerComponent[compj][blockId] / static_cast<double>(cyclesInBlock);
      matrix[adjustCompi][adjustCompj] = Average_NxNj - Average_N * Average_Nj;
    }
  }
  GaussJordan(matrix, tempMatrix);

  size_t adjustedComponent = component - static_cast<size_t>(SystemComponents.NComponents.y);
  heat = 0.0;
  for(size_t compj = static_cast<size_t>(SystemComponents.NComponents.y); compj < static_cast<size_t>(SystemComponents.NComponents.x); compj++)
  {
    double Average_N = SystemComponents.Moves[compj].MolAverage[blockId].x / static_cast<double>(cyclesInBlock);
    double Average_ExN = SystemComponents.EnergyTimesNumberOfMolecule[compj][blockId] / static_cast<double>(cyclesInBlock);
    size_t adjustedCompj = compj - static_cast<size_t>(SystemComponents.NComponents.y);
    double inverseVariance = matrix[adjustedCompj][adjustedComponent];
    if(!std::isfinite(inverseVariance)) return false;
    heat += Constants.energy_to_kelvin * (Average_ExN - Average_E * Average_N) * inverseVariance;
  }

  double kelvin_to_kjmol = 0.01 / Constants.energy_to_kelvin;
  heat -= SystemComponents.Temperature;
  heat *= kelvin_to_kjmol;
  return std::isfinite(heat);
}

static inline void PrintAdaptiveConvergenceStatus(Components& SystemComponents, const AdaptiveProductionSettings& Settings)
{
  fprintf(SystemComponents.OUTPUT, "Adaptive convergence check: mode=%s, cycles=%zu, batches=%zu, consecutive_passes=%zu/%zu\n",
          Settings.CriteriaMode == ADAPTIVE_CRITERIA_ALL ? "All" : "Any",
          SystemComponents.AdaptiveState.CyclesCompleted,
          SystemComponents.AdaptiveState.CompletedBlocks,
          SystemComponents.AdaptiveState.ConsecutivePasses,
          Settings.ConsecutivePasses);

  for(size_t comp = static_cast<size_t>(SystemComponents.NComponents.y); comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
  {
    if(comp < SystemComponents.AdaptiveState.HenryStatus.size() && SystemComponents.AdaptiveState.HenryStatus[comp].Enabled)
    {
      const AdaptiveObservableStatus& Status = SystemComponents.AdaptiveState.HenryStatus[comp];
      fprintf(SystemComponents.OUTPUT,
              "  Henry[%s]: valid_batches=%zu, mean=%.10g, halfwidth=%.10g, rel_halfwidth=%.10g, status=%s\n",
              SystemComponents.MoleculeName[comp].c_str(),
              Status.ValidBatches,
              Status.Mean,
              Status.HalfWidth,
              Status.RelativeHalfWidth,
              Status.Passed ? "pass" : "keep-running");
    }
    if(comp < SystemComponents.AdaptiveState.HeatStatus.size() && SystemComponents.AdaptiveState.HeatStatus[comp].Enabled)
    {
      const AdaptiveObservableStatus& Status = SystemComponents.AdaptiveState.HeatStatus[comp];
      fprintf(SystemComponents.OUTPUT,
              "  HeatOfAdsorption[%s]: valid_batches=%zu, mean=%.10g, halfwidth=%.10g, rel_halfwidth=%.10g, status=%s\n",
              SystemComponents.MoleculeName[comp].c_str(),
              Status.ValidBatches,
              Status.Mean,
              Status.HalfWidth,
              Status.RelativeHalfWidth,
              Status.Passed ? "pass" : "keep-running");
    }
  }
}

static inline bool EvaluateAdaptiveProduction(Variables& Vars, size_t systemId, size_t totalCycles)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  const AdaptiveProductionSettings& Settings = Vars.AdaptiveProduction;
  size_t completedBlocks = CeilDivide(totalCycles, Settings.BatchCycles);

  SystemComponents.AdaptiveState.CyclesCompleted = totalCycles;
  SystemComponents.AdaptiveState.CompletedBlocks = completedBlocks;

  bool hasTargets = false;
  bool allTargetsPassed = true;
  for(size_t comp = static_cast<size_t>(SystemComponents.NComponents.y); comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
  {
    if(comp < SystemComponents.AdaptiveState.HenryStatus.size() && SystemComponents.AdaptiveState.HenryStatus[comp].Enabled)
    {
      hasTargets = true;
      std::vector<double> estimates;
      estimates.reserve(completedBlocks);
      for(size_t blockId = 0; blockId < completedBlocks; blockId++)
      {
        double henry = 0.0;
        if(ComputeHenryEstimateForBlock(SystemComponents, Vars.Constants, Settings, comp, blockId, totalCycles, Settings.BatchCycles, henry))
          estimates.push_back(henry);
      }
      AdaptiveObservableStatus Status = BuildAdaptiveObservableStatus(estimates, Settings, Settings.HenryAbsoluteTolerance);
      Status.Enabled = true;
      SystemComponents.AdaptiveState.HenryStatus[comp] = Status;
      allTargetsPassed = allTargetsPassed && Status.Passed;
    }

    if(comp < SystemComponents.AdaptiveState.HeatStatus.size() && SystemComponents.AdaptiveState.HeatStatus[comp].Enabled)
    {
      hasTargets = true;
      std::vector<double> estimates;
      estimates.reserve(completedBlocks);
      for(size_t blockId = 0; blockId < completedBlocks; blockId++)
      {
        double heat = 0.0;
        if(ComputeHeatEstimateForBlock(SystemComponents, Vars.Constants, comp, blockId, totalCycles, Settings.BatchCycles, heat))
          estimates.push_back(heat);
      }
      AdaptiveObservableStatus Status = BuildAdaptiveObservableStatus(estimates, Settings, Settings.HeatAbsoluteTolerance);
      Status.Enabled = true;
      SystemComponents.AdaptiveState.HeatStatus[comp] = Status;
      allTargetsPassed = allTargetsPassed && Status.Passed;
    }
  }

  if(!hasTargets)
    throw std::runtime_error("Adaptive production is enabled, but no components requested MonitorHenryCoefficient or MonitorHeatOfAdsorption");

  if(totalCycles >= static_cast<size_t>(Vars.NumberOfProductionCycles) && allTargetsPassed)
    SystemComponents.AdaptiveState.ConsecutivePasses++;
  else
    SystemComponents.AdaptiveState.ConsecutivePasses = 0;

  SystemComponents.AdaptiveState.Converged = SystemComponents.AdaptiveState.ConsecutivePasses >= Settings.ConsecutivePasses;
  PrintAdaptiveConvergenceStatus(SystemComponents, Settings);
  return SystemComponents.AdaptiveState.Converged;
}

static inline void PrintAdaptiveProductionSummary(Components& SystemComponents, const AdaptiveProductionSettings& Settings)
{
  fprintf(SystemComponents.OUTPUT, "=====================ADAPTIVE PRODUCTION SUMMARY=====================\n");
  fprintf(SystemComponents.OUTPUT, "Adaptive Production: %s\n", Settings.Enabled ? "enabled" : "disabled");
  if(Settings.Enabled)
  {
    fprintf(SystemComponents.OUTPUT, "Criteria Mode: %s\n", Settings.CriteriaMode == ADAPTIVE_CRITERIA_ALL ? "All" : "Any");
    fprintf(SystemComponents.OUTPUT, "Cycles Completed: %zu\n", SystemComponents.AdaptiveState.CyclesCompleted);
    fprintf(SystemComponents.OUTPUT, "Batches Completed: %zu (batch size: %zu cycles)\n", SystemComponents.AdaptiveState.CompletedBlocks, Settings.BatchCycles);
    fprintf(SystemComponents.OUTPUT, "Converged: %s\n", SystemComponents.AdaptiveState.Converged ? "yes" : "no");
    fprintf(SystemComponents.OUTPUT, "Stop Reason: %s\n", SystemComponents.AdaptiveState.StopReason.c_str());
    for(size_t comp = static_cast<size_t>(SystemComponents.NComponents.y); comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
    {
      if(comp < SystemComponents.AdaptiveState.HenryStatus.size() && SystemComponents.AdaptiveState.HenryStatus[comp].Enabled)
      {
        const AdaptiveObservableStatus& Status = SystemComponents.AdaptiveState.HenryStatus[comp];
        fprintf(SystemComponents.OUTPUT,
                "Henry[%s]: valid_batches=%zu, mean=%.10g, halfwidth=%.10g, rel_halfwidth=%.10g, passed=%s\n",
                SystemComponents.MoleculeName[comp].c_str(),
                Status.ValidBatches,
                Status.Mean,
                Status.HalfWidth,
                Status.RelativeHalfWidth,
                Status.Passed ? "yes" : "no");
      }
      if(comp < SystemComponents.AdaptiveState.HeatStatus.size() && SystemComponents.AdaptiveState.HeatStatus[comp].Enabled)
      {
        const AdaptiveObservableStatus& Status = SystemComponents.AdaptiveState.HeatStatus[comp];
        fprintf(SystemComponents.OUTPUT,
                "HeatOfAdsorption[%s]: valid_batches=%zu, mean=%.10g, halfwidth=%.10g, rel_halfwidth=%.10g, passed=%s\n",
                SystemComponents.MoleculeName[comp].c_str(),
                Status.ValidBatches,
                Status.Mean,
                Status.HalfWidth,
                Status.RelativeHalfWidth,
                Status.Passed ? "yes" : "no");
      }
    }
  }
  fprintf(SystemComponents.OUTPUT, "========================================================================\n");
}

///////////////////////////////////////////////////////////
// Wrapper for Performing a move for the selected system //
///////////////////////////////////////////////////////////
void Select_Box_Component_Molecule(Variables& Vars, size_t box_index)
{
  Components& SystemComponents = Vars.SystemComponents[box_index];
  WidomStruct& Widom = Vars.Widom[box_index];
  SystemComponents.TempVal.Initialize();
  size_t& comp                   = SystemComponents.TempVal.component;
  size_t& SelectedMolInComponent = SystemComponents.TempVal.molecule;
  
  //Randomly Select an Adsorbate Molecule and determine its Component: MoleculeID --> Component
  //Zhao's note: The number of atoms can be vulnerable, adding throw error here//
  if(SystemComponents.TotalNumberOfMolecules < SystemComponents.NumberOfFrameworks)
    throw std::runtime_error("There is negative number of adsorbates. Break program!");

  size_t NumberOfImmobileFrameworkMolecules = 0; size_t ImmobileFrameworkSpecies = 0;
  for(size_t i = 0; i < SystemComponents.NComponents.y; i++)
    if(SystemComponents.Moves[i].TotalProb < 1e-10)
    {
      ImmobileFrameworkSpecies++;
      NumberOfImmobileFrameworkMolecules += SystemComponents.NumberOfMolecule_for_Component[i];
    }
  while(SystemComponents.Moves[comp].TotalProb < 1e-10)
  {
    comp = (size_t) (Get_Uniform_Random() * SystemComponents.NComponents.x);
  }
  SelectedMolInComponent = (size_t) (Get_Uniform_Random() * SystemComponents.NumberOfMolecule_for_Component[comp]);

  Vars.RandomNumber = Get_Uniform_Random();
}
void RunMoves(Variables& Vars, size_t box_index, int Cycle)
{
  MC_MOVES MOVES;

  Components& SystemComponents = Vars.SystemComponents[box_index];
  Simulations& Sims = Vars.Sims[box_index];
  ForceField& FF = Vars.device_FF;
  //RandomNumber& Random = Vars.Random;
  WidomStruct& Widom = Vars.Widom[box_index];

  //variables that affects the selection of a move, written into TempVal//
  Select_Box_Component_Molecule(Vars, box_index);
  double& RANDOMNUMBER = Vars.RandomNumber;
  size_t& comp         = SystemComponents.TempVal.component;
  size_t& SelectedMolInComponent = SystemComponents.TempVal.molecule;
  //printf("Step %zu, selected Comp %zu, Mol %zu, RANDOM: %.5f", Cycle, comp, SelectedMolInComponent, RANDOMNUMBER);

  MoveEnergy DeltaE;
  int& MoveType = SystemComponents.TempVal.MoveType;
  if(RANDOMNUMBER < SystemComponents.Moves[comp].TranslationProb)
  {
    MoveType = TRANSLATION;
    //////////////////////////////
    // PERFORM TRANSLATION MOVE //
    //////////////////////////////
    //printf(" Translation\n");
    if(SystemComponents.NumberOfMolecule_for_Component[comp] > 0)
    {
      DeltaE = SingleBodyMove(Vars, box_index);
    }
    else
    {
      SystemComponents.Tmmc[comp].Update(1.0, SystemComponents.NumberOfMolecule_for_Component[comp], TRANSLATION);
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].RotationProb) //Rotation
  {
    MoveType = ROTATION;
    ///////////////////////////
    // PERFORM ROTATION MOVE //
    ///////////////////////////
    //printf(" Rotation\n");
    if(SystemComponents.NumberOfMolecule_for_Component[comp] > 0)
    {
      DeltaE = SingleBodyMove(Vars, box_index);
    }
    else
    {
      SystemComponents.Tmmc[comp].Update(1.0, SystemComponents.NumberOfMolecule_for_Component[comp], ROTATION);
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].SpecialRotationProb) //Special Rotation for Framework Components
  {
    MoveType = SPECIAL_ROTATION;
    ///////////////////////////////////
    // PERFORM SPECIAL ROTATION MOVE //
    ///////////////////////////////////
    //printf(" Special Rotation\n");
    if(SystemComponents.NumberOfMolecule_for_Component[comp] > 0)
      DeltaE = SingleBodyMove(Vars, box_index);
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].WidomProb)
  {
    MoveType = WIDOM;
    //////////////////////////////////
    // PERFORM WIDOM INSERTION MOVE //
    //////////////////////////////////
    //printf(" Widom Insertion\n");
    double2& newScale = SystemComponents.TempVal.Scale; 
    newScale = SystemComponents.Lambda[comp].SET_SCALE(1.0); //Set scale for full molecule (lambda = 1.0)//
    double Rosenbluth = MOVES.INSERTION.WidomMove(Vars, box_index);
    //Also record move energy (delta energy)
    //MoveEnergy widom_e = MOVES.INSERTION.energy;
    //printf("Widom E: "); widom_e.print();

    if(Vars.SimulationMode == PRODUCTION)
    {
      size_t blockID = Cycle/Vars.BlockAverageSize;
      if(blockID >= SystemComponents.Nblock) blockID --;
      SystemComponents.Moves[comp].BlockID = blockID;
      SystemComponents.Moves[comp].RecordRosen(Rosenbluth, WIDOM);
      //weight with Rosenbluth (heavy weight on low (attractive) energy, lower on high energy)
      SystemComponents.Moves[comp].Rosen[blockID].widom_energy += MOVES.INSERTION.energy * Rosenbluth;
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].ReinsertionProb)
  {
    //////////////////////////////
    // PERFORM REINSERTION MOVE //
    //////////////////////////////
    //printf(" Reinsertion\n");
    MoveType = REINSERTION;
    if(SystemComponents.NumberOfMolecule_for_Component[comp] > 0)
    {
      //DeltaE = Reinsertion(Vars, box_index);
      DeltaE = MOVES.REINSERTION.Run(Vars, box_index);
    }
    else
    {
      SystemComponents.Tmmc[comp].Update(1.0, SystemComponents.NumberOfMolecule_for_Component[comp], REINSERTION);
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].IdentitySwapProb)
  {
    MoveType = IDENTITY_SWAP;
    //printf(" Identity Swap\n");
    DeltaE = IdentitySwapMove(Vars, box_index);
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].CBCFProb && SystemComponents.hasfractionalMolecule[comp])
  {
    ///////////////////////
    // PERFORM CBCF MOVE //
    ///////////////////////
    //printf(" CBCF\n");
    SelectedMolInComponent = SystemComponents.Lambda[comp].FractionalMoleculeID;
    DeltaE = CBCFMove(Vars, box_index);
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].SwapProb)
  {
    ////////////////////////////
    // PERFORM GCMC INSERTION //
    ////////////////////////////
    if(Get_Uniform_Random() < 0.5)
    {
      //printf(" Swap Insertion\n");
      if(!SystemComponents.SingleSwap)
      {
        MoveType = INSERTION;
        DeltaE = MOVES.INSERTION.Run(Vars, box_index);
        //DeltaE = Insertion(Vars, box_index);
      }
      else
      {
        MoveType = SINGLE_INSERTION;
        DeltaE = SingleBodyMove(Vars, box_index);
        //DeltaE = SingleSwapMove(SystemComponents, Sims, Widom, FF, Random, SelectedMolInComponent, comp, SINGLE_INSERTION);
      }
    }
    else
    {
      ///////////////////////////
      // PERFORM GCMC DELETION //
      ///////////////////////////
      //printf(" Swap Deletion\n");
      //Zhao's note: Do not do a deletion if the chosen molecule is a fractional molecule, fractional molecules should go to CBCFSwap moves//
      if(!((SystemComponents.hasfractionalMolecule[comp]) && SelectedMolInComponent == SystemComponents.Lambda[comp].FractionalMoleculeID))
      {
        if(SystemComponents.NumberOfMolecule_for_Component[comp] > 0)
        {
          if(!SystemComponents.SingleSwap)
          {
            MoveType = DELETION;
            DeltaE = MOVES.DELETION.Run(Vars, box_index);
            //DeltaE = Deletion(Vars, box_index);
          }
          else
          {
            MoveType = SINGLE_DELETION;
            DeltaE = SingleBodyMove(Vars, box_index);
          }
        }
        else
        {
          MoveType = DELETION;
          SystemComponents.Tmmc[comp].Update(0.0, SystemComponents.NumberOfMolecule_for_Component[comp], DELETION);
        }
      }
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].VolumeMoveProb)
  {
    //printf(" VOLUME MOVE\n");
    double start = omp_get_wtime();
    ForceField& FF = Vars.device_FF;
    VolumeMove(SystemComponents, Sims, FF);
    double end = omp_get_wtime();
    SystemComponents.VolumeMoveTime += end - start;
  }
  //Gibbs Xfer//
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].GibbsSwapProb)
  {
    //if(Vars.GibbsStatistics.DoGibbs)
    //printf(" Gibbs SWAP\n");
    if(Vars.SystemComponents.size() == 2)
    {
      //GibbsParticleTransfer(Vars, comp, Vars.GibbsStatistics);
      MOVES.GIBBS_PARTICLE_XFER.Run(Vars, box_index, Vars.GibbsStatistics);
    }
  }
  else if(RANDOMNUMBER < SystemComponents.Moves[comp].GibbsVolumeMoveProb)
  {
    //printf(" Gibbs VOLUME\n");
    if(Vars.SystemComponents.size() == 2)
      NVTGibbsMove(Vars.SystemComponents, Vars.Sims, FF, Vars.GibbsStatistics);
  }
  SystemComponents.deltaE += DeltaE;
}

double CreateMolecule_InOneBox(Variables& Vars, size_t systemId, bool AlreadyHasFractionalMolecule)
{
  MC_MOVES MOVES; 
  Components& SystemComponents = Vars.SystemComponents[systemId];
  //Simulations& Sims = Vars.Sims[systemId];
  //ForceField& FF = Vars.device_FF;
  //RandomNumber& Random = Vars.Random;
  //WidomStruct& Widom = Vars.Widom[systemId];
  double running_energy = 0.0;
  // Create Molecules in the Box Before the Simulation //
  for(size_t comp = SystemComponents.NComponents.y; comp < SystemComponents.NComponents.x; comp++)
  {
    size_t CreateFailCount = 0; size_t Created = 0; size_t SelectedMol = 0;
    CreateFailCount = 0;
    fprintf(SystemComponents.OUTPUT, "Component %zu, Need to create %zu full molecule\n", comp, SystemComponents.NumberOfCreateMolecules[comp]);
    //Create Fractional Molecule first//
    if(SystemComponents.hasfractionalMolecule[comp])
    {
      //Zhao's note: If we need to create fractional molecule, then we initialize WangLandau Histogram//
      size_t FractionalMolToCreate = 1;
      if(AlreadyHasFractionalMolecule) FractionalMolToCreate = 0;
      if(FractionalMolToCreate > 0) Initialize_WangLandauIteration(SystemComponents.Lambda[comp]);
      while(FractionalMolToCreate > 0)
      {
        fprintf(SystemComponents.OUTPUT, "Creating Fractional Molecule for Component %zu; There are %zu Molecules of that component in the System\n", comp, SystemComponents.NumberOfMolecule_for_Component[comp]);
        SelectedMol = Created; if(Created > 0) SelectedMol = Created - 1; 
        //Zhao's note: this is a little confusing, but when number of molecule for that species = 0 or 1, the chosen molecule is zero. This is creating from zero loading, need to change in the future, when we read from restart file//
        size_t OldVal = SystemComponents.NumberOfMolecule_for_Component[comp];

        size_t NewBin = 5;
        MoveEnergy DeltaE;
        if(SystemComponents.Tmmc[comp].DoTMMC) NewBin = 0;
        double newLambda = static_cast<double>(NewBin) * SystemComponents.Lambda[comp].delta;
        SystemComponents.TempVal.Initialize();
	SystemComponents.TempVal.Scale = SystemComponents.Lambda[comp].SET_SCALE(newLambda);
        SystemComponents.TempVal.MoveType = INSERTION;
        SystemComponents.TempVal.component = comp;
        SystemComponents.TempVal.molecule  = SelectedMol;
        DeltaE = MOVES.INSERTION.CreateMolecule(Vars, systemId);
        running_energy += DeltaE.total();
        SystemComponents.CreateMoldeltaE += DeltaE;
        if(SystemComponents.NumberOfMolecule_for_Component[comp] == OldVal)
        {
          CreateFailCount ++;
        }
        else
        {
          FractionalMolToCreate --; Created ++; SystemComponents.Lambda[comp].FractionalMoleculeID = SelectedMol;
          SystemComponents.Lambda[comp].currentBin = NewBin;
        }
        if(CreateFailCount > 1e20) throw std::runtime_error("Bad Insertions When Creating Fractional Molecules!");
      }
    }
    while(SystemComponents.NumberOfCreateMolecules[comp] > 0)
    {
      fprintf(SystemComponents.OUTPUT, "Creating %zu Molecule for Component %zu; There are %zu Molecules of that component in the System\n", Created, comp, SystemComponents.NumberOfMolecule_for_Component[comp]);
      SelectedMol = Created; if(Created > 0) SelectedMol = Created - 1; //Zhao's note: this is a little confusing, but when number of molecule for that species = 0 or 1, the chosen molecule is zero. This is creating from zero loading, need to change in the future, when we read from restart file//
      size_t OldVal    = SystemComponents.NumberOfMolecule_for_Component[comp];
      MoveEnergy DeltaE;
      SystemComponents.TempVal.Initialize();
      SystemComponents.TempVal.Scale = SystemComponents.Lambda[comp].SET_SCALE(1.0); //Set scale for full molecule (lambda = 1.0)//
      SystemComponents.TempVal.MoveType = INSERTION;
      SystemComponents.TempVal.component = comp;
      SystemComponents.TempVal.molecule  = SelectedMol;
      DeltaE = MOVES.INSERTION.CreateMolecule(Vars, systemId);
      //printf("Creating %zu molecule\n", SelectedMol);
      //DeltaE.print();
      running_energy += DeltaE.total();
      SystemComponents.CreateMoldeltaE += DeltaE;
      fprintf(SystemComponents.OUTPUT, "Delta E in creating molecules:\n"); DeltaE.print();
      if(SystemComponents.NumberOfMolecule_for_Component[comp] == OldVal)
      {CreateFailCount ++;} else {SystemComponents.NumberOfCreateMolecules[comp] --; Created ++;}
      if(CreateFailCount > 1e10) throw std::runtime_error("Bad Insertions When Creating Molecules!");
    }
  }
  return running_energy;
}

void GatherStatisticsDuringSimulation(Variables& Vars, size_t systemId, size_t cycle)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  Simulations&  Sims           = Vars.Sims[systemId];
  size_t& i = cycle;
  int& BlockAverageSize        = Vars.BlockAverageSize;
  int& SimulationMode          = Vars.SimulationMode;
  std::string& Mode = Vars.Mode;
  //////////////////////////////////////////////
  // SAMPLE (EQUILIBRATION) CBCF BIASING TERM //
  //////////////////////////////////////////////
  if(SimulationMode == EQUILIBRATION && i%50==0)
  {
    for(size_t icomp = 0; icomp < SystemComponents.NComponents.x; icomp++)
    { //Try to sample it if there are more CBCF moves performed//
      if(SystemComponents.hasfractionalMolecule[icomp] && !SystemComponents.Tmmc[icomp].DoTMMC)
      {
        Sample_WangLandauIteration(SystemComponents.Lambda[icomp]);
        SystemComponents.CBCFPerformed[icomp] = SystemComponents.Moves[icomp].CBCFTotal; 
        SystemComponents.WLSampled++;
      }
    }
  }

  if(i%500==0)
  {
    for(size_t comp = 0; comp < SystemComponents.NComponents.x; comp++)
    {  
      if(SystemComponents.Moves[comp].TranslationTotal > 0)
        Update_Max_Translation(SystemComponents, comp);
      if(SystemComponents.Moves[comp].RotationTotal > 0)
        Update_Max_Rotation(SystemComponents, comp);
      if(SystemComponents.Moves[comp].SpecialRotationTotal > 0)
        Update_Max_SpecialRotation(SystemComponents, comp);
      if(SystemComponents.VolumeMoveAttempts > 0) Update_Max_VolumeChange(SystemComponents);
    }
  }
  if(i%SystemComponents.PrintStatsEvery==0) Print_Cycle_Statistics(i, SystemComponents, Mode);
  ////////////////////////////////////////////////
  // ADJUST CBCF BIASING FACTOR (EQUILIBRATION) //
  ////////////////////////////////////////////////
  if(i%5000==0 && SimulationMode == EQUILIBRATION)
  {
    for(size_t icomp = 0; icomp < SystemComponents.NComponents.x; icomp++)
      if(SystemComponents.hasfractionalMolecule[icomp] && !SystemComponents.Tmmc[icomp].DoTMMC)//Try not to use CBCFC + TMMC//
      {  Adjust_WangLandauIteration(SystemComponents.Lambda[icomp]); SystemComponents.WLAdjusted++;}
  }
  if(SimulationMode == PRODUCTION)
  {
    //Record values for Number of atoms//
    for(size_t comp = 0; comp < SystemComponents.NComponents.x; comp++)
    {
      Gather_Averages_Types(SystemComponents.Moves[comp].MolAverage, SystemComponents.NumberOfMolecule_for_Component[comp], 0.0, i, BlockAverageSize, SystemComponents.Nblock);
      //Gather total energy * number of molecules for each adsorbate component//
      if(comp >= SystemComponents.NComponents.y)
      {
        double deltaE_Adsorbate = SystemComponents.deltaE.total() - SystemComponents.deltaE.HHVDW - SystemComponents.deltaE.HHEwaldE - SystemComponents.deltaE.HHReal;
        double ExN = SystemComponents.createmol_energy + deltaE_Adsorbate * SystemComponents.NumberOfMolecule_for_Component[comp];
        Gather_Averages_double(SystemComponents.EnergyTimesNumberOfMolecule[comp], ExN, i, BlockAverageSize, SystemComponents.Nblock);
        //Calculate Average Excess Loading//
        //AmountOfExcessMolecules only be resized during EOS calculation, don't have that? then no excess loading because excess loading needs compressibility from EOS//
        if(SystemComponents.AmountOfExcessMolecules.size() > 0)
          Gather_Averages_Types(SystemComponents.ExcessLoading[comp], SystemComponents.NumberOfMolecule_for_Component[comp] - SystemComponents.AmountOfExcessMolecules[comp], 0.0, i, BlockAverageSize, SystemComponents.Nblock);
      }
      for(size_t compj = 0; compj < SystemComponents.NComponents.x; compj++)
      {
        if(comp >= SystemComponents.NComponents.y && compj >= SystemComponents.NComponents.y)
        {
          double NxNj = SystemComponents.NumberOfMolecule_for_Component[comp] * SystemComponents.NumberOfMolecule_for_Component[compj];
          Gather_Averages_double(SystemComponents.Moves[comp].MolSQPerComponent[compj], NxNj, i, BlockAverageSize, SystemComponents.Nblock);
        }
      }
      Gather_Averages_Types(SystemComponents.DensityPerComponent[comp], SystemComponents.NumberOfMolecule_for_Component[comp] / Sims.Box.Volume, 0.0, i, BlockAverageSize, SystemComponents.Nblock);
    }
    Gather_Averages_Types(SystemComponents.VolumeAverage, Sims.Box.Volume, 0.0, i, BlockAverageSize, SystemComponents.Nblock);
    Gather_Averages_MoveEnergy(SystemComponents, i, BlockAverageSize, SystemComponents.deltaE);
  }
  if(SimulationMode != INITIALIZATION && i > 0)
  {
    for(size_t comp = 0; comp < SystemComponents.NComponents.x; comp++)
      if(i % SystemComponents.Tmmc[comp].UpdateTMEvery == 0)
        SystemComponents.Tmmc[comp].AdjustTMBias();
  }
  if(i % SystemComponents.MoviesEvery == 0)//Generate restart file and movies 
    GenerateRestartMovies(Vars, systemId, SystemComponents.PseudoAtoms, SimulationMode);
}

void InitialMCBeforeMoves(Variables& Vars, size_t systemId)
{
  Components& SystemComponents = Vars.SystemComponents[systemId];
  size_t NumberOfSimulations   = Vars.SystemComponents.size();
  int&   BlockAverageSize      = Vars.BlockAverageSize;
  int&   Cycles                = Vars.Cycles;

  SystemComponents.CBCFPerformed.resize(SystemComponents.NComponents.x);
  SystemComponents.WLSampled = 0; SystemComponents.WLAdjusted = 0;

  fprintf(SystemComponents.OUTPUT, "==================================\n");
  std::string& Mode   = Vars.Mode;
  int& SimulationMode = Vars.SimulationMode;
  switch(SimulationMode)
  {
    case INITIALIZATION: {Mode = "INITIALIZATION"; fprintf(SystemComponents.OUTPUT, "== RUNNING INITIALIZATION PHASE ==\n"); Cycles = Vars.NumberOfInitializationCycles; break;}
    case EQUILIBRATION:  {Mode = "EQUILIBRATION";  fprintf(SystemComponents.OUTPUT, "== RUNNING EQUILIBRATION PHASE ==\n");  Cycles = Vars.NumberOfEquilibrationCycles; break;}
    case PRODUCTION:
    {
      Mode = "PRODUCTION";
      fprintf(SystemComponents.OUTPUT, "==  RUNNING PRODUCTION PHASE   ==\n");
      if(AdaptiveProductionEnabled(Vars))
      {
        if(Vars.RunTogether)
          throw std::runtime_error("Adaptive production is currently supported only for one-box-at-a-time production runs");
        if(Vars.AdaptiveProduction.MaximumCycles == 0)
          Vars.AdaptiveProduction.MaximumCycles = static_cast<size_t>(Vars.NumberOfProductionCycles);
        if(Vars.AdaptiveProduction.BatchCycles == 0)
          Vars.AdaptiveProduction.BatchCycles = CeilDivide(std::max(1, Vars.NumberOfProductionCycles), std::max<size_t>(1, SystemComponents.ConfiguredNblock));
        if(Vars.AdaptiveProduction.MaximumCycles < static_cast<size_t>(Vars.NumberOfProductionCycles))
          throw std::runtime_error("MaximumProductionCycles must be greater than or equal to NumberOfProductionCycles when adaptive production is enabled");
        if(Vars.AdaptiveProduction.BatchCycles == 0)
          throw std::runtime_error("AdaptiveBatchCycles must be greater than ZERO");
        if(Vars.AdaptiveProduction.MinimumBatches == 0)
          throw std::runtime_error("AdaptiveMinimumBatches must be greater than ZERO");
        if(Vars.AdaptiveProduction.ConsecutivePasses == 0)
          throw std::runtime_error("AdaptiveConsecutivePasses must be greater than ZERO");
        if(Vars.AdaptiveProduction.RelativeTolerance <= 0.0 &&
           Vars.AdaptiveProduction.HenryAbsoluteTolerance < 0.0 &&
           Vars.AdaptiveProduction.HeatAbsoluteTolerance < 0.0)
          throw std::runtime_error("Adaptive production needs either AdaptiveRelativeTolerance or an absolute tolerance for Henry/Heat");
        Cycles = CheckedSizeToInt(Vars.AdaptiveProduction.MaximumCycles, "MaximumProductionCycles");
      }
      else
      {
        Cycles = Vars.NumberOfProductionCycles;
      }
      break;
    }
  }
  fprintf(SystemComponents.OUTPUT, "==================================\n");

  fprintf(SystemComponents.OUTPUT, "CBMC Uses %zu trial positions and %zu trial orientations\n", Vars.Widom[systemId].NumberWidomTrials, Vars.Widom[systemId].NumberWidomTrialsOrientations);

  if(SimulationMode == INITIALIZATION)
  {
    fprintf(SystemComponents.OUTPUT, "Box %zu, Volume: %.5f\n", systemId, Vars.Sims[systemId].Box.Volume);
    Vars.GibbsStatistics.TotalVolume += Vars.Sims[systemId].Box.Volume;

    fprintf(SystemComponents.OUTPUT, "Total Volume: %.5f\n", Vars.GibbsStatistics.TotalVolume);
  }
  // Kaihang Shi: Record initial energy but exclude the host-host Ewald
  SystemComponents.createmol_energy = SystemComponents.CreateMol_Energy.total() - SystemComponents.CreateMol_Energy.HHVDW - SystemComponents.CreateMol_Energy.HHEwaldE - SystemComponents.CreateMol_Energy.HHReal;

  if(SimulationMode == PRODUCTION)
  {
    if(AdaptiveProductionEnabled(Vars))
    {
      bool hasAdaptiveTarget = false;
      for(size_t comp = static_cast<size_t>(SystemComponents.NComponents.y); comp < static_cast<size_t>(SystemComponents.NComponents.x); comp++)
      {
        if(comp < SystemComponents.AdaptiveTargets.size() && SystemComponents.AdaptiveTargets[comp].MonitorHenryCoefficient)
        {
          hasAdaptiveTarget = true;
          if(Vars.AdaptiveProduction.RelativeTolerance <= 0.0 && Vars.AdaptiveProduction.HenryAbsoluteTolerance < 0.0)
            throw std::runtime_error("MonitorHenryCoefficient requires AdaptiveRelativeTolerance or AdaptiveAbsoluteToleranceHenry");
          if(SystemComponents.Moves[comp].WidomProb < 1e-10)
            throw std::runtime_error("MonitorHenryCoefficient requires a non-zero WidomProbability for the same component");
        }
        if(comp < SystemComponents.AdaptiveTargets.size() && SystemComponents.AdaptiveTargets[comp].MonitorHeatOfAdsorption)
        {
          hasAdaptiveTarget = true;
          if(Vars.AdaptiveProduction.RelativeTolerance <= 0.0 && Vars.AdaptiveProduction.HeatAbsoluteTolerance < 0.0)
            throw std::runtime_error("MonitorHeatOfAdsorption requires AdaptiveRelativeTolerance or AdaptiveAbsoluteToleranceHeat");
        }
      }
      if(!hasAdaptiveTarget)
        throw std::runtime_error("Adaptive production is enabled, but no adaptive observables were requested in component blocks");

      size_t adaptiveBlocks = CeilDivide(static_cast<size_t>(Cycles), Vars.AdaptiveProduction.BatchCycles);
      BlockAverageSize = CheckedSizeToInt(Vars.AdaptiveProduction.BatchCycles, "AdaptiveBatchCycles");
      InitializeProductionStatisticsStorage(SystemComponents, adaptiveBlocks);
      ResetAdaptiveProductionState(SystemComponents);
      fprintf(SystemComponents.OUTPUT,
              "Adaptive production enabled: minimum cycles=%d, maximum cycles=%d, batch cycles=%d, allocated batches=%zu\n",
              Vars.NumberOfProductionCycles,
              Cycles,
              BlockAverageSize,
              adaptiveBlocks);
    }
    else
    {
      size_t requestedBlocks = std::max<size_t>(1, SystemComponents.ConfiguredNblock);
      size_t usableBlocks = std::min(requestedBlocks, static_cast<size_t>(std::max(1, Cycles)));
      if(usableBlocks != requestedBlocks)
        fprintf(SystemComponents.OUTPUT, "Warning! Number of production cycles is smaller than NumberOfBlocks. Reducing block count from %zu to %zu\n", requestedBlocks, usableBlocks);
      BlockAverageSize = std::max(1, Cycles / CheckedSizeToInt(usableBlocks, "NumberOfBlocks"));
      if(Cycles % CheckedSizeToInt(usableBlocks, "NumberOfBlocks") != 0)
        fprintf(SystemComponents.OUTPUT, "Warning! Number of Cycles cannot be divided by Number of blocks. Residue values go to the last block\n");
      InitializeProductionStatisticsStorage(SystemComponents, usableBlocks);
      SystemComponents.AdaptiveState = AdaptiveProductionState{};
    }
  }

  /////////////////////////////////////////////
  // FINALIZE (PRODUCTION) CBCF BIASING TERM //
  /////////////////////////////////////////////
  if(SimulationMode == PRODUCTION)
  {
    for(size_t icomp = 0; icomp < SystemComponents.NComponents.x; icomp++)
      if(SystemComponents.hasfractionalMolecule[icomp] && !SystemComponents.Tmmc[icomp].DoTMMC)
        Finalize_WangLandauIteration(SystemComponents.Lambda[icomp]);
  }

  ///////////////////////////////////////////////////////////////////////
  // FORCE INITIALIZING CBCF BIASING TERM BEFORE INITIALIZATION CYCLES //
  ///////////////////////////////////////////////////////////////////////
  if(SimulationMode == INITIALIZATION && Cycles > 0)
  {
    for(size_t icomp = 0; icomp < SystemComponents.NComponents.x; icomp++)
      if(SystemComponents.hasfractionalMolecule[icomp])
        Initialize_WangLandauIteration(SystemComponents.Lambda[icomp]);
  }
  
  if(SimulationMode == EQUILIBRATION) //Rezero the TMMC stats at the beginning of the Equilibration cycles//
  {
    for(size_t comp = 0; comp < SystemComponents.NComponents.x; comp++)
    {
      //Clear TMMC data in the collection matrix//
      SystemComponents.Tmmc[comp].ClearCMatrix();
      //Clear Rosenbluth weight statistics after Initialization//
      for(size_t i = 0; i < SystemComponents.Nblock; i++)
        SystemComponents.Moves[comp].ClearRosen(i);
    }
  }
}

inline void MCEndOfPhaseSummary(Variables& Vars)
{
  std::vector<Components>& SystemComponents = Vars.SystemComponents;
  Simulations*&  Sims   = Vars.Sims;
  Units& Constants = Vars.Constants;
  std::string& Mode = Vars.Mode;

  size_t NumberOfSimulations = SystemComponents.size();
  int& Cycles = Vars.Cycles;

  //print statistics
  if(Cycles > 0)
  {
    for(size_t sim = 0; sim < NumberOfSimulations; sim++)
    {
      if(Vars.SimulationMode == EQUILIBRATION) fprintf(SystemComponents[sim].OUTPUT, "Sampled %zu WangLandau, Adjusted WL %zu times\n", SystemComponents[sim].WLSampled, SystemComponents[sim].WLAdjusted);
      PrintAllStatistics(SystemComponents[sim], Sims[sim], Cycles, Vars.SimulationMode, Vars.BlockAverageSize, Constants);
      if(Vars.SimulationMode == PRODUCTION)
      {
        Calculate_Overall_Averages_MoveEnergy(SystemComponents[sim], Vars.BlockAverageSize, Cycles);
        if(Vars.AdaptiveProduction.Enabled)
          PrintAdaptiveProductionSummary(SystemComponents[sim], Vars.AdaptiveProduction);
      }
    }
    PrintSystemMoves(Vars);
  }
  for(size_t i = 0; i < Vars.SystemComponents.size(); i++)
  {
    fprintf(SystemComponents[i].OUTPUT, "===============================\n");
    fprintf(SystemComponents[i].OUTPUT, "== %s PHASE ENDS ==\n", Mode.c_str());
    fprintf(SystemComponents[i].OUTPUT, "===============================\n");
  }
}
//Default is 20 steps per cycle//
//If # of molecules > 20, use # of molecules//
//If a max limit is imposed, use the max limit if it exceeds//
size_t Determine_Number_Of_Steps(Variables& Vars, size_t systemId, size_t current_cycle)
{ 
  //Record current step//
  Vars.SystemComponents[systemId].CURRENTCYCLE = current_cycle;
  size_t Steps = 20;
  if(Steps < Vars.SystemComponents[systemId].TotalNumberOfMolecules)
  {
    Steps = Vars.SystemComponents[systemId].TotalNumberOfMolecules;
  }
  if(Vars.SetMaxStep && Steps > Vars.MaxStepPerCycle) Steps = Vars.MaxStepPerCycle;
  return Steps;
}

void Run_Simulation_MultipleBoxes(Variables& Vars)
{
  if(AdaptiveProductionEnabled(Vars))
    throw std::runtime_error("Adaptive production is currently unsupported when simulations are run together");

  std::vector<Components>&   SystemComponents = Vars.SystemComponents;
  size_t NumberOfSimulations = SystemComponents.size();

  for(size_t sim = 0; sim < NumberOfSimulations; sim++)
    InitialMCBeforeMoves(Vars, sim);

  ///////////////////////////////////////////////////////
  // Run the simulations for different boxes IN SERIAL //
  ///////////////////////////////////////////////////////
  for(size_t i = 0; i < Vars.Cycles; i++)
  {
    for(size_t sim = 0; sim < NumberOfSimulations; sim++)
    {
      double RNM = Get_Uniform_Random();
      size_t selectedSim = static_cast<size_t>(RNM * static_cast<double>(NumberOfSimulations));
      size_t Steps = Determine_Number_Of_Steps(Vars, selectedSim, i);
      //printf("STEPS: %zu, RNM: %.5f, selectedSim: %zu\n", Steps, RNM, selectedSim);
      for(size_t j = 0; j < Steps; j++)
      {
        RunMoves(Vars, selectedSim, i);
      }
    }
    for(size_t sim = 0; sim < NumberOfSimulations; sim++)
    {
      GatherStatisticsDuringSimulation(Vars, sim, i);
    }
    if(i > 0 && i % 500 == 0)
      Update_Max_GibbsVolume(Vars.GibbsStatistics);
  }
  MCEndOfPhaseSummary(Vars);
}

void Run_Simulation_ForOneBox(Variables& Vars, size_t box_index)
{
  InitialMCBeforeMoves(Vars, box_index);
  Components& SystemComponents = Vars.SystemComponents[box_index];
  bool adaptiveProduction = AdaptiveProductionEnabled(Vars);
  size_t cyclesCompleted = 0;

  for(size_t i = 0; i < Vars.Cycles; i++)
  {
    size_t Steps = Determine_Number_Of_Steps(Vars, box_index, i);
    for(size_t j = 0; j < Steps; j++)
    {
      RunMoves(Vars, box_index, i);
    }
    GatherStatisticsDuringSimulation(Vars, box_index, i);
    cyclesCompleted = i + 1;

    if(adaptiveProduction)
    {
      bool completedBatch = (cyclesCompleted % Vars.AdaptiveProduction.BatchCycles) == 0;
      bool reachedMaximumCycles = cyclesCompleted == static_cast<size_t>(Vars.Cycles);
      if(completedBatch || reachedMaximumCycles)
      {
        bool converged = EvaluateAdaptiveProduction(Vars, box_index, cyclesCompleted);
        if(converged)
        {
          SystemComponents.AdaptiveState.StopReason = "Convergence criteria satisfied";
          break;
        }
        if(reachedMaximumCycles)
          SystemComponents.AdaptiveState.StopReason = "Reached MaximumProductionCycles before satisfying convergence criteria";
      }
    }
  }

  if(adaptiveProduction)
  {
    size_t usedBlocks = cyclesCompleted > 0 ? CeilDivide(cyclesCompleted, Vars.AdaptiveProduction.BatchCycles) : 0;
    SystemComponents.Nblock = usedBlocks;
    Vars.Cycles = CheckedSizeToInt(cyclesCompleted, "completed production cycles");
    if(SystemComponents.AdaptiveState.StopReason.empty())
      SystemComponents.AdaptiveState.StopReason = "Production phase finished";
  }
  MCEndOfPhaseSummary(Vars);
}
