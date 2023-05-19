#include "SelfPropulsionForce.hpp"

template<unsigned DIM>
SelfPropulsionForce<DIM>::SelfPropulsionForce()
  : AbstractForce<DIM>(),
    // Self propulsion force magnitude. Each vertex has contributions
    // with this magnitude (but different directions) from the three
    // cells which it connects.
    mF0(0.0)
{
}

template<unsigned DIM>
SelfPropulsionForce<DIM>::~SelfPropulsionForce()
{
}

template<unsigned DIM>
void SelfPropulsionForce<DIM>::SetF0(double newValue)
{
  mF0 = newValue;
}

template<unsigned DIM>
double SelfPropulsionForce<DIM>::GetF0()
{
  return mF0;
}

template<unsigned DIM>
void SelfPropulsionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
  // Throw an exception message if not using a VertexBasedCellPopulation
  if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
      EXCEPTION("SelfPropulsionForce is to be used with a VertexBasedCellPopulation only");
    }

  // Define some helper variables
  VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
  unsigned num_nodes = p_cell_population->GetNumNodes();
  unsigned num_elements = p_cell_population->GetNumElements();

  // Get the current propulsion angle for each cell from CellData
  std::vector<double> propulsion_theta(num_elements);
  for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
       elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
       ++elem_iter)
    {
      unsigned elem_index = elem_iter->GetIndex();
      try
        {
	  // Throw an exception if self propulsion angle Theta isn't
	  // specified in the CellData. See
	  // TestERKWaveWithSelfPropulsionNoAlignment for an example
	  // of how to initialize. There is a slight chance that the
	  // exception is thrown although the error is not about the
	  // propulsion thetas.
	  propulsion_theta[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("Theta");
        }
      catch (Exception&)
        {
	  EXCEPTION("CellData needs to contain an entry for 'Theta' for SelfPropulsionForce to work");
        }
    }

  // Iterate over vertices in the cell population
  for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {

      c_vector<double, DIM> propulsion_contribution = zero_vector<double>(DIM);

      // Find the indices of the elements (cells) owned by this node
      std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

      // Iterate over these elements
      for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
	   iter != containing_elem_indices.end();
	   ++iter)
        {
	  // Get this element and its index
	  VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
	  unsigned elem_index = p_element->GetIndex();

	  // Add the propulsion contributions from each cell (element)
	  // associated with the node
	  propulsion_contribution[0] += cos(propulsion_theta[elem_index]);
	  propulsion_contribution[1] += sin(propulsion_theta[elem_index]);
	}
      // Multiply the unit vector by the force magnitude
      propulsion_contribution = GetF0()*propulsion_contribution;
      p_cell_population->GetNode(node_index)->AddAppliedForceContribution(propulsion_contribution);
    }
}

template<unsigned DIM>
void SelfPropulsionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<F0>" << mF0 << "</F0>\n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SelfPropulsionForce<2>;
template class SelfPropulsionForce<3>;    // Might work in three dimensions but only adds force in x-y plane

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SelfPropulsionForce)
