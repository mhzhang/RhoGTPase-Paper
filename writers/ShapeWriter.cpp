/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 * Only supports vertex based simulations
 */

#include "ShapeWriter.hpp"
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
ShapeWriter<ELEMENT_DIM, SPACE_DIM>::ShapeWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("shape_data.csv")
{
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	*this->mpOutStream << "# TimeStamp,Triangle,Quadrilateral,Pentagon,Hexagon,Heptagon\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{
	
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << " ";
    
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ShapeWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
	
    int edges[5] = {};
	
	pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();

	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		VertexElement < SPACE_DIM, SPACE_DIM > *VertexElement = pCellPopulation->GetElementCorrespondingToCell(*cell_iter);
		int num_edges = VertexElement->GetNumNodes();
		edges[num_edges-3]+=1;
	}
	for (int i = 0; i < sizeof(edges)/sizeof(edges[0]); i ++){
    	*this->mpOutStream << "," << edges[i];
	}	

}

// Explicit instantiation
template class ShapeWriter<1,1>;
template class ShapeWriter<1,2>;
template class ShapeWriter<2,2>;
template class ShapeWriter<1,3>;
template class ShapeWriter<2,3>;
template class ShapeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ShapeWriter)