/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 * Only supports vertex based simulations
 */

#include "NumNeighboursWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCellPopulationWriter.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VertexElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::NumNeighboursWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("neighbour_data.csv")
{
}

/* Write CSV Header */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	*this->mpOutStream << "# TimeStamp,1Neighbour,2Neighbours,3Neighbours,4Neighbours,5Neighbours,6Neighbours,7Neighbours\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
	
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << " ";
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NumNeighboursWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
	
	pCellPopulation->Update();

	int neighbours[7] = {};

	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		std::set<unsigned> neighbour_indices = pCellPopulation->GetNeighbouringLocationIndices(*cell_iter);
		int num_neighbours = neighbour_indices.size();
		neighbours[num_neighbours-1] += 1;
	}
	
	for (int i = 0; i < sizeof(neighbours)/sizeof(neighbours[0]); i ++){
    	*this->mpOutStream  << ","<< neighbours[i];
	}	
	
}

// Explicit instantiation
template class NumNeighboursWriter<1,1>;
template class NumNeighboursWriter<1,2>;
template class NumNeighboursWriter<2,2>;
template class NumNeighboursWriter<1,3>;
template class NumNeighboursWriter<2,3>;
template class NumNeighboursWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NumNeighboursWriter)