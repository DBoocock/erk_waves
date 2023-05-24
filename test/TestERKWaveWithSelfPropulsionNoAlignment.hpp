/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#ifndef TESTERKWAVEWITHSELFPROPULSIONNOALIGNMENT_HPP_
#define TESTERKWAVEWITHSELFPROPULSIONNOALIGNMENT_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "ToroidalHoneycombVertexMeshGenerator2.hpp"    // Modified to give unit size cells
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

#include "TargetAreaAndPerimeterForce.hpp"
#include "SelfPropulsionForce.hpp"
#include "ErkPropulsionSrnModelNoAlignment.hpp"
#include "ErkPropulsionModifierNoAlignment.hpp"
#include "ErkPropulsionWriterNoAlignment.hpp"

#include "CommandLineArguments.hpp"
#include <iostream>
#include <boost/filesystem.hpp>


std::string get_timestamp()
{
  auto now = std::time(nullptr);
  char buf[sizeof("YYYY-MM-DD  HH:MM:SS")];
  return std::string(buf, buf + std::strftime(buf, sizeof(buf), "%F_%T", std::gmtime(&now)));
}

class TestERKWaveWithSelfPropulsion : public AbstractCellBasedTestSuite
{
public:
  void TestRunSimulation()
    {
      // Vertex based simulations cannot be run in parallel
      EXIT_IF_PARALLEL;

      // Read in parameters from command line

      // A random seed for the initial noise and dynamics of
      // persistent random walk
      int seed = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-seed");
      RandomNumberGenerator::Instance()->Reseed(seed);

      // Name for output directory (with root directory set by
      // environment variable CHASTE_TEST_OUTPUT)
      std::string outdir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-outdir");

      // System size in number of cells
      int nx = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-nx");
      int ny = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-ny");

      // Timesteps for vertex model and ODE solver (I set these the same)
      double dt = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");
      double dt_ode = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt_ode");

      // Simulation end time
      double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");
      // How often write data to file
      int sampling_timestep_multiple = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-sampling_timestep_multiple");

      // We avoid recording the initial transitiant we run one
      // simulation until end_time (saving only at start and end) then
      // load and simulate again for bonus_time, sampling at the rate
      // set by the sampling_timestep_multiple.
      double bonus_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-bonus_time");

      // Initial mean values of variables
      double init_erk = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-init_erk");
      double init_A0 = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-init_A0");
      // Non dimensionalize length using the average cell area
      double init_A = 1.0;

      // Timescales: in units of tau_E (Note that taur enters through
      // the area elasticity KA and is not predicted to affect the
      // onset of instability set by ab_crit)
      double taul = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-taul");
      // Self-propulsion timescales: in units of tau_E
      double tau_p = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-tau_p");

      // Standard deviation of random noise in self-propulsion angle
      // calculated from tau_p below. Corrections corresponging to a
      // change in time-step are accounted for in the ODE system.
      double eta_std = 1/sqrt(tau_p);

      // Values of area and preimneter elasticities KA and KP
      // pre-scaled by the substrate friction coefficient.
      double KA = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-KA");
      double KP = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-KP");
      // Preferred perimeter (dimensionles shape parameter with this
      // scaling) determing the solid-fluid / rigid-floppy rheology of
      // the tissue.
      double P0 = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-P0");
      // Self-propulsion force pre-scaled by the subrate friction coefficient
      double F0 = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-F0");

      // Value of mechanochemical coupling strength alpha*beta in
      // terms of the critical value for onset of pattern formaion
      // predicted by linear stability
      double n_abcrit = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-n_abcrit");
      // Ratio of mechanochemical coupling strengths
      double ab_ratio = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-ab_ratio");

      // Calculate the critical value of alpha*beta for onset of
      // pattern formaion predicted by linear stability.  All times
      // should be non-dimensionalized by taue
      double abcrit = (taul + 1.0)*pow(sqrt(taul) + sqrt(1.0), 2) / (1.0 * taul);
      double ab = n_abcrit*abcrit;    // Alpha*beta

      // Determine alpha and beta from alpha*beta and alpha/beta
      double beta;
      double alpha;
      if (ab != 0.0) {
	beta = sqrt(ab/ab_ratio);
	alpha = ab/beta;
      } else {
	// Broken mechanochemcial coupling is specified by input value
	// ab=0. Set beta=1 so that we can still output ERK's slave
	// response to area fluctuations induced by the persistent
	// random walk.
	alpha = 0.0;
	beta = 1.0;
      }

      // Size of initial pertubation in vertex positions and ERK. An
      // initial perturbation is required in order to destabilize the
      // system away from the steady state and allow onset of pattern
      // formation (especially when there is no self-propulsion).
      double init_noise = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-init_noise");
      double noiseSD_pos;
      double noiseSD_erk;
      if (beta >= alpha) {
	// Large beta leaves ERK very sensitive to changes in cell
	// area so scale initial niose in vertex positions by the
	// value of beta. Alternatively we could have set noise on
	// vertex positions (areas) to zero and just pertubed ERK.
	noiseSD_pos = init_noise/beta;
	noiseSD_erk = 1E-10;
      }
      else if (alpha != 0.0) {
	// Large alpha leaves A0 very sensitive to changes in ERK so
	// scale initial niose in ERK by the value of
	// alpha. Alternatively we could have set noise on ERK to zero
	// and just pertubed vertex positions (areas).
	noiseSD_pos = 1E-10;    // Small noise to escape regular hexagonal tiling
	noiseSD_erk = init_noise/alpha;
      }
      else {    // alpha=beta=0 (i.e. ERK slave to random persistent walk)
	noiseSD_pos = init_noise;
	noiseSD_erk = 0.0;
      }
      // Easier to account for changes in mechanochemical couplings
      // alpha/beta through noise in vertex positions and ERK so set
      // initial perturbation on A0 to zero.
      double noiseSD_A0 = 0.0;

      // Checks for geometric errors in simulation, e.g. missing T1s
      // may lead to overlappling (intersecting) cells.
      bool check_for_internal_intersections = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-check_for_internal_intersections");

      // Save all parameters to file.
      std::string outdirpath = std::string(getenv("CHASTE_TEST_OUTPUT")) + "/" + outdir;
      boost::filesystem::create_directories(outdirpath);
      std::string outpath = outdirpath + "/params.txt";
      std::ofstream myfile;
      myfile.open(outpath);
      myfile << "seed " << std::to_string(seed) << std::endl
	     << "nx " << std::to_string(nx) << std::endl
	     << "ny " << std::to_string(ny) << std::endl
	     << "dt " << std::to_string(dt) << std::endl
	     << "dt_ode " << std::to_string(dt_ode) << std::endl
	     << "end_time " << std::to_string(end_time) << std::endl
	     << "sampling_timestep_multiple " << std::to_string(sampling_timestep_multiple) << std::endl
	     << "bonus_time " << std::to_string(bonus_time) << std::endl
	     << "init_erk " << std::to_string(init_erk) << std::endl
	     << "init_A0 " << std::to_string(init_A0) << std::endl
	     << "init_A " << std::to_string(init_A) << std::endl
	     << "init_noise " << std::to_string(init_noise) << std::endl
	     << "taue " << std::to_string(1.0) << std::endl
	     << "taul " << std::to_string(taul) << std::endl
	     << "tau_p " << std::to_string(tau_p) << std::endl
	     << "eta_std " << std::to_string(eta_std) << std::endl
	     << "F0 " << std::to_string(F0) << std::endl
	     << "KA " << std::to_string(KA) << std::endl
	     << "KP " << std::to_string(KP) << std::endl
	     << "P0 " << std::to_string(P0) << std::endl
	     << "n_abcrit " << std::to_string(n_abcrit) << std::endl
	     << "abcrit " << std::to_string(abcrit) << std::endl
	     << "ab_ratio " << std::to_string(ab_ratio) << std::endl
	     << "ab " << std::to_string(ab) << std::endl
	     << "alpha " << std::to_string(alpha) << std::endl
	     << "beta " << std::to_string(beta) << std::endl
	     << "check_for_internal_intersections " << std::to_string(check_for_internal_intersections) << std::endl;
      myfile.close();

      // Set up the vertex model

      // Create a vertex mesh of regular hexagonal cells with unit
      // area and perturb by adding noise to the initial vertex
      // positions.
      ToroidalHoneycombVertexMeshGenerator2 generator(nx, ny, init_A, noiseSD_pos);
      Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

      // Add a check for overlapping cells (can happen when T1s aren't
      // properly detected.
      p_mesh->SetCheckForInternalIntersections(check_for_internal_intersections);

      // Create and initialize cells
      std::vector<CellPtr> cells;
      // We are required to specify a mutation and cell type so choose
      // the most basic
      MAKE_PTR(WildTypeCellMutationState, p_state);
      // Differentiated cells without cell cycle (no proliferation / cell division)
      MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

      for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
	  // Set initial values of variables for this cell in the SRN
	  // model. Must push back in correct order (theta, ERK, A0)
	  // corresponding to the order of variables in the
	  // ErkPropulsionOdeSystemNoAlignment.
	  std::vector<double> initial_conditions;
	  // Initial random value for self propulsion angle between
	  // -Pi and Pi
	  double theta;
	  theta = (RandomNumberGenerator::Instance()->ranf()*2 - 1)*M_PI;
	  initial_conditions.push_back(theta);    // Theta
	  // Initial values of ERK and A0 plus optional noise
	  initial_conditions.push_back(RandomNumberGenerator::Instance()->NormalRandomDeviate(init_erk, noiseSD_erk));    // Erk
	  initial_conditions.push_back(RandomNumberGenerator::Instance()->NormalRandomDeviate(init_A0, noiseSD_A0));    // Target area
	  ErkPropulsionSrnModelNoAlignment* p_srn_model = new ErkPropulsionSrnModelNoAlignment();
	  p_srn_model->SetDt(dt_ode);
	  p_srn_model->SetInitialConditions(initial_conditions);

	  // Create the cell. Simulations require specification of a cell cycle model
	  // but we set so create one but set non-proliferative cells
	  UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
	  p_cc_model->SetDimension(2);
	  CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
	  p_cell->SetCellProliferativeType(p_diff_type);
	  // No cell division or death so birth time should not matter
	  p_cell->SetBirthTime(0.0);

	  // Set initial values of variables for the cell.
	  p_cell->GetCellData()->SetItem("Theta", initial_conditions[0]);    // Variable
	  p_cell->GetCellData()->SetItem("Erk", initial_conditions[1]);    // Variable
	  p_cell->GetCellData()->SetItem("Target Area", initial_conditions[2]);    // Variable
	  // Set parameters for the cell. These are passed to the ODE
	  // solver by the ErkPropulsionSrnModel. TODO: Find a way to
	  // specify global parameters only once rather than for every
	  // cell.
	  p_cell->GetCellData()->SetItem("taul", taul);
	  p_cell->GetCellData()->SetItem("alpha", alpha);
	  p_cell->GetCellData()->SetItem("beta", beta);
	  // Self-propulsion
	  p_cell->GetCellData()->SetItem("Eta Std", eta_std);
	  // The size of the timestep is used to scale the noise
	  // variance in the SDE so that the timestep does not affect
	  // the persistence time or the persistent random walk.
	  p_cell->GetCellData()->SetItem("dt_ode", dt_ode);
	  cells.push_back(p_cell);
        }

      // Create a cell-based population object, and specify which
      // results to output to file.
      VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
      cell_population.AddCellWriter<ErkPropulsionWriterNoAlignment>();
      for (typename VertexBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
	   cell_iter != cell_population.End();
	   ++cell_iter)
	{
	  // Get the area (2D volume) of this cell and initialize the
	  // value in the CellData.
	  double cell_volume = cell_population.GetVolumeOfCell(*cell_iter);
	  cell_iter->GetCellData()->SetItem("volume", cell_volume);
	}

      // Create and configure the cell-based simulation object
      OffLatticeSimulation<2> simulator(cell_population);
      simulator.SetOutputDirectory(outdir);
      simulator.SetDt(dt);
      // Only save the start and end points of the burn-in period
      simulator.SetSamplingTimestepMultiple(end_time/dt);
      simulator.SetEndTime(end_time);

      // Add a modifier that keeps track of variables and updates them
      // between the solver and the CellData.
      MAKE_PTR(ErkPropulsionModifierNoAlignment<2>, p_modifier);
      simulator.AddSimulationModifier(p_modifier);

      // Add self propulsion force to use with the persistent random
      // walk in the ODE system / SRN model
      MAKE_PTR(SelfPropulsionForce<2>, p_propulsion_force);
      p_propulsion_force->SetF0(F0);
      simulator.AddForce(p_propulsion_force);

      // Add force with target area and perimeter
      MAKE_PTR(TargetAreaAndPerimeterForce<2>, p_force);
      p_force->SetKA(KA);
      p_force->SetKP(KP);
      p_force->SetP0(P0);
      simulator.AddForce(p_force);

      // Run the simulation over the burn-in period
      simulator.Solve();
      CellBasedSimulationArchiver<2, OffLatticeSimulation<2>>::Save(&simulator);

      // Now load from the checkpoint created after the burn-in period
      // and record data at more frequenct intervals
      OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2>>::Load(outdir, end_time);
      p_simulator->SetSamplingTimestepMultiple(sampling_timestep_multiple);
      p_simulator->SetEndTime(end_time+bonus_time);
      p_simulator->Solve();
      CellBasedSimulationArchiver<2, OffLatticeSimulation<2>>::Save(p_simulator);

    }
};

#endif /*TESTERKWAVEWITHSELFPROPULSIONNOALIGNMENT_HPP_*/
