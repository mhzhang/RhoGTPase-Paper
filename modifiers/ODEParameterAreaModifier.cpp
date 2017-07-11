/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 */

#include "ODEParameterAreaModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "AbstractOdeSrnModel.hpp"
#include "AbstractSrnModel.hpp"

template<unsigned DIM>
ODEParameterAreaModifier<DIM>::ODEParameterAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ODEParameterAreaModifier<DIM>::~ODEParameterAreaModifier()
{
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the volume of this cell
		double cell_volume = cell_iter->GetCellData()->GetItem("volume");
		
        AbstractSrnModel* model_ptr_base = cell_iter->GetSrnModel();
		AbstractOdeSrnModel* model_ptr = dynamic_cast<AbstractOdeSrnModel*>(model_ptr_base);
		model_ptr->GetStateVariables()[2]=cell_volume;
    }
}

template<unsigned DIM>
void ODEParameterAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ODEParameterAreaModifier<1>;
template class ODEParameterAreaModifier<2>;
template class ODEParameterAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ODEParameterAreaModifier)