#ifndef SELFPROPULSIONFORCE_HPP_
#define SELFPROPULSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * A 2d 'self propulsion force' class. Also works for 3d but will only
 * add propulsion in the x-y plane.
 *
 * The angle theta of self propulsion is read from CellData. This
 * force class was developed as part of simulations of a persistent
 * random walk which used an SDE system
 * (ErkPropulsionOdeSystemNoAlignment) to update the self propulsion
 * angle theta and a modifier class (ErkPropulsionModifierNoAlignment)
 * to update theta in CellData.
 *
 * Works with vertex based cell populations by averaging the
 * propulsion forces from the cells connected by a cell junction onto
 * that cell junction.
 */
template<unsigned DIM>
class SelfPropulsionForce : public AbstractForce<DIM>
{
private :
  /**
   * propulsion force magnitude
   */
  double mF0;

  /**
   * Archiving.
   */
  friend class boost::serialization::access;
  /**
   * Boost Serialization method for archiving/checkpointing.
   * Archives the object and its member variables.
   *
   * @param archive  The boost archive.
   * @param version  The current version of this class.
   */
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
    archive & mF0;
  }

public :

  /**
   * Constructor.
   */
  SelfPropulsionForce();

  /**
   * Destructor.
   */
  ~SelfPropulsionForce();

  /**
   * Set the propulsion force magnitude F0.
   *
   * @param force
   */
  void SetF0(double force);

  /**
   * Get the propulsion force magnitude F0.
   *
   * @return mF0
   */
  double GetF0();

  /**
   * Overridden AddForceContribution() method.
   *
   * @param rCellPopulation reference to the tissue
   */
  void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

  /**
   * Overridden OutputForceParameters() method.
   *
   * @param rParamsFile the file stream to which the parameters are output
   */
  void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SelfPropulsionForce)

#endif /*SELFPROPULSIONFORCE_HPP_*/
