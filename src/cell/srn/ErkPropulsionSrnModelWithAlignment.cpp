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
#include "ErkPropulsionSrnModelWithAlignment.hpp"

ErkPropulsionSrnModelWithAlignment::ErkPropulsionSrnModelWithAlignment(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
  : AbstractOdeSrnModel(3, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
      // Evolution of self propulsion angle follows a stochastic
      // differential equation so use a basic Euler solver rather than
      // one solver with an addaptive timestep.
      mpOdeSolver = CellCycleModelOdeSolver<ErkPropulsionSrnModelWithAlignment, EulerIvpOdeSolver>::Instance();
      mpOdeSolver->Initialise();
 }
    assert(mpOdeSolver->IsSetUp());
}

ErkPropulsionSrnModelWithAlignment::ErkPropulsionSrnModelWithAlignment(const ErkPropulsionSrnModelWithAlignment& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
  assert(rModel.GetOdeSystem());
  SetOdeSystem(new ErkPropulsionOdeSystemWithAlignment(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* ErkPropulsionSrnModelWithAlignment::CreateSrnModel()
{
    return new ErkPropulsionSrnModelWithAlignment(*this);
}

void ErkPropulsionSrnModelWithAlignment::SimulateToCurrentTime()
{
    // Update areas from CellData to the ode system
    UpdateSrnAreas();
    // Update the average theta (angle of self propulsion) of the neighbours
    // from CellData to the ode system
    UpdateSrnMeanThetas();
    // Run the ODE simulation
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void ErkPropulsionSrnModelWithAlignment::Initialise()
{
    AbstractOdeSrnModel::Initialise(new ErkPropulsionOdeSystemWithAlignment);
    SetSrnParams();    // Pass parameters from CellData to the OdeSystem
}

void ErkPropulsionSrnModelWithAlignment::UpdateSrnAreas()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    // Send the cell area (2D volume) to the ode solver
    double cell_area = mpCell->GetCellData()->GetItem("volume");
    mpOdeSystem->SetParameter("Cell Area", cell_area);
}

void ErkPropulsionSrnModelWithAlignment::UpdateSrnMeanThetas()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    // Send the cell volumes to the ode solver
    double mean_theta = mpCell->GetCellData()->GetItem("Mean Theta");
    mpOdeSystem->SetParameter("Mean Theta", mean_theta);
}

void ErkPropulsionSrnModelWithAlignment::SetSrnParams()
{
  assert(mpOdeSystem != nullptr);
  assert(mpCell != nullptr);

  double taul = mpCell->GetCellData()->GetItem("taul");
  mpOdeSystem->SetParameter("taul", taul);

  // Send the cell volumes to the ode solver
  double alpha = mpCell->GetCellData()->GetItem("alpha");
  mpOdeSystem->SetParameter("alpha", alpha);

  double beta = mpCell->GetCellData()->GetItem("beta");
  mpOdeSystem->SetParameter("beta", beta);

  double eta_std = mpCell->GetCellData()->GetItem("Eta Std");
  mpOdeSystem->SetParameter("Eta Std", eta_std);

  double taua = mpCell->GetCellData()->GetItem("taua");
  mpOdeSystem->SetParameter("taua", taua);

  double mean_theta = mpCell->GetCellData()->GetItem("Mean Theta");
  mpOdeSystem->SetParameter("Mean Theta", mean_theta);

  double dt_ode = mpCell->GetCellData()->GetItem("dt_ode");
  mpOdeSystem->SetParameter("dt_ode", dt_ode);

}

double ErkPropulsionSrnModelWithAlignment::GetTheta()
{
    assert(mpOdeSystem != nullptr);
    double theta = mpOdeSystem->rGetStateVariables()[0];
    return theta;
}

double ErkPropulsionSrnModelWithAlignment::GetErk()
{
    assert(mpOdeSystem != nullptr);
    double erk = mpOdeSystem->rGetStateVariables()[1];
    return erk;
}

double ErkPropulsionSrnModelWithAlignment::GetTargetArea()
{
    assert(mpOdeSystem != nullptr);
    double target_area = mpOdeSystem->rGetStateVariables()[2];
    return target_area;
}

// Is this an unused and broken function that was only relevant for
// the original ERK simulation?
double ErkPropulsionSrnModelWithAlignment::GetCellArea()
{
    assert(mpOdeSystem != nullptr);
    double cell_area = mpOdeSystem->GetParameter("Cell Area");
    return cell_area;
}

void ErkPropulsionSrnModelWithAlignment::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ErkPropulsionSrnModelWithAlignment)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(ErkPropulsionSrnModelWithAlignment)
