/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 * Only supports vertex based simulations
 */

#include "OneCellGTPaseWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VertexElement.hpp"
#include <VertexBasedCellPopulation.hpp>
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::OneCellGTPaseWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("OneCellData.csv")
{
}

// Writing CSV Header
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	
	*this->mpOutStream << "# TimeStamp,Cell_ID,Cell_Area,Target_Area,Cell_Perimeter,num_neighbours,num_edges,G,A_t,A\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{

    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << ",";
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OneCellGTPaseWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
	
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		
		double cell_id = cell_iter->GetCellId();
		
		if (cell_id == 1){
			
			double volume = cell_iter->GetCellData()->GetItem("volume");
			unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
			double perimeter = pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(elem_index);
			std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringLocationIndices(*cell_iter);
			int num_neighbours = neighbour_indices.size();
			
			VertexElement < SPACE_DIM, SPACE_DIM > *VertexElement = pCellPopulation->GetElementCorrespondingToCell(*cell_iter);
			int num_edges = VertexElement->GetNumNodes();
			double G = cell_iter->GetCellData()->GetItem("G");
			double a = cell_iter->GetCellData()->GetItem("AREA");
			double target_area = cell_iter->GetCellData()->GetItem("target area");
			
		    *this->mpOutStream  << ","<< cell_id;
			*this->mpOutStream  << ","<< volume;
			*this->mpOutStream  << ","<< target_area;
			*this->mpOutStream  << ","<< perimeter;
			*this->mpOutStream  << ","<< num_neighbours;
			*this->mpOutStream  << ","<< num_edges;
			*this->mpOutStream  << ","<< G;
			*this->mpOutStream  << ","<< target_area;
			*this->mpOutStream  << ","<< a;
			
		}

	}
}

// Explicit instantiation
template class OneCellGTPaseWriter<1,1>;
template class OneCellGTPaseWriter<1,2>;
template class OneCellGTPaseWriter<2,2>;
template class OneCellGTPaseWriter<1,3>;
template class OneCellGTPaseWriter<2,3>;
template class OneCellGTPaseWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OneCellGTPaseWriter)