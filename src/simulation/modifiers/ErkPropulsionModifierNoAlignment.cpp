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

#include "ErkPropulsionModifierNoAlignment.hpp"
#include "ErkPropulsionSrnModelNoAlignment.hpp"

template<unsigned DIM>
ErkPropulsionModifierNoAlignment<DIM>::ErkPropulsionModifierNoAlignment()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ErkPropulsionModifierNoAlignment<DIM>::~ErkPropulsionModifierNoAlignment()
{
}

template<unsigned DIM>
void ErkPropulsionModifierNoAlignment<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
    double current_time = SimulationTime::Instance()->GetTime();
    std::cout << "time " << current_time << std::endl;
}

template<unsigned DIM>
void ErkPropulsionModifierNoAlignment<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * Update CellData in SetupSolve(), to assure that it has been
     * fully initialised before we enter the main time loop.
     */
    // UpdateCellData(rCellPopulation);    // Already doing this explicitly during setup in TestERKWaveWithSelfPropulsionNoAlignment.hpp
}

template<unsigned DIM>
void ErkPropulsionModifierNoAlignment<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    rCellPopulation.Update();

    // Iterate over the population to compute and update each cell's
    // data from the ODE solver to the population CellData.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      // Get the current values of ERK and theta for this cell
      ErkPropulsionSrnModelNoAlignment* p_model = static_cast<ErkPropulsionSrnModelNoAlignment*>(cell_iter->GetSrnModel());

      // Get the ERK value for this cell from the solver and update in
      // CellData
      double this_erk = p_model->GetErk();
      cell_iter->GetCellData()->SetItem("Erk", this_erk);

      // Get the target area for this cell from the solver and update
      // in CellData
      double this_target_area = p_model->GetTargetArea();
      cell_iter->GetCellData()->SetItem("Target Area", this_target_area);

      // Get the area for this cell from the solver and update in
      // CellData
      double cell_volume = rCellPopulation.GetVolumeOfCell(*cell_iter);
      cell_iter->GetCellData()->SetItem("volume", cell_volume);

      // Get the self-propulsion angle for this cell from the solver and update
      // in CellData
      double this_theta = p_model->GetTheta();
      // Map theta to the range -Pi and Pi for better visualization in
      // VTK output. Makes no difference for the ODE system which
      // takes the sine and cosine.
      this_theta = atan2(sin(this_theta), cos(this_theta));
      cell_iter->GetCellData()->SetItem("Theta", this_theta);
    }
}

template<unsigned DIM>
void ErkPropulsionModifierNoAlignment<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ErkPropulsionModifierNoAlignment<1>;
template class ErkPropulsionModifierNoAlignment<2>;
template class ErkPropulsionModifierNoAlignment<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ErkPropulsionModifierNoAlignment)
