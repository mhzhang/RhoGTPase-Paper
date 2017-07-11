/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 * Only supports vertex based simulations
 */

#include "CsvWriter.hpp"
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
CsvWriter<ELEMENT_DIM, SPACE_DIM>::CsvWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("data.csv")
{
}

// Writing CSV Header
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::WriteHeader(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
	
	*this->mpOutStream << "# TimeStamp,Total_Number_Of_Cells,Average_Area,Average_Perimeter,Num_Labelled\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::VisitAnyPopulation(AbstractCellPopulation<SPACE_DIM, SPACE_DIM>* pCellPopulation)
{

    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	

	*this->mpOutStream << num_cells << ",";
    
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = static_cast<MutableMesh<ELEMENT_DIM,SPACE_DIM>&>((pCellPopulation->rGetMesh())).GetVolume();

    
    *this->mpOutStream << num_cells << " ";

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    VisitAnyPopulation(pCellPopulation);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CsvWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    pCellPopulation->Update();

    unsigned num_cells = pCellPopulation->GetNumRealCells();
	double total_area = 0;
	double total_perimeter = 0;
	double labelled_count = 0;
	
	for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
	cell_iter != pCellPopulation->End();
	++cell_iter){
		
		double volume = cell_iter->GetCellData()->GetItem("volume");
		total_area+=volume;
		
		unsigned elem_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
		double perimeter = pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(elem_index);
		total_perimeter += perimeter;
		
		if (cell_iter->template HasCellProperty<CellLabel>()){
			labelled_count +=1;
		}
	}
	
	double avg_area = total_area/num_cells;
	double avg_perimeter = total_perimeter/num_cells;
    *this->mpOutStream  << ","<< num_cells;
	*this->mpOutStream  << ","<< avg_area;
	*this->mpOutStream  << ","<< avg_perimeter;
	*this->mpOutStream  << ","<< labelled_count;
}

// Explicit instantiation
template class CsvWriter<1,1>;
template class CsvWriter<1,2>;
template class CsvWriter<2,2>;
template class CsvWriter<1,3>;
template class CsvWriter<2,3>;
template class CsvWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CsvWriter)