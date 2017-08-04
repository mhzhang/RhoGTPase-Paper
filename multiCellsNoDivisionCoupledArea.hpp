#ifndef MULTICELLSNODIVISIONCOUPLED_HPP_
#define MULTICELLSNODIVISIONCOUPLED_HPP_

/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: Aug 3, 2017
 * Do not reproduce this code without permission.
 */

#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "RepulsionForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "VoronoiDataWriter.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SimulationTime.hpp"

#include "CellLabel.hpp"

#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "ContactInhibitionCellCycleModel.hpp"

#include "VolumeTrackingModifier.hpp"
#include "ODEParameterAreaModifier.hpp"
#include "ShapeWriter.hpp"
#include "CsvWriter.hpp"
#include "NumNeighboursWriter.hpp"
#include "XMLCellWriter.hpp"
#include "OneCellGTPaseWriter.hpp"

#include "ODESRNCoupledArea.hpp"
#include "AbstractOdeSrnModel.hpp"

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include <cmath>

class multiCellsNoDivisionCoupled : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedMonolayer() throw (Exception)
    {
		

        /* 2500 cells */

        HoneycombVertexMeshGenerator generator(50, 50);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
		
	/* Seed RNG */
	RandomNumberGenerator* randGenerator = RandomNumberGenerator::Instance();
	randGenerator->Reseed(1);

        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
	    ODESrnModel* p_srn_model = new ODESrnModel;
	    
	    /* Random Initial GTPase Concentration */
	    double G_conc = (randGenerator->ranf());
	    
	    std::vector<double> initial_conditions;

	    /* Set initial conditions for ODE solver */
	    initial_conditions.push_back(G_conc);
	    // initial target area
	    initial_conditions.push_back(0.8);
	    // initial cell area
	    initial_conditions.push_back(0.866025);
	    p_srn_model->SetInitialConditions(initial_conditions);
			
            p_cycle_model->SetDimension(2);
	    /* Prevent cell division - all cells out of M phase */
            p_cycle_model->SetBirthTime(-(double)i - 2.0);
            p_cycle_model->SetQuiescentVolumeFraction(1.0);
            p_cycle_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
			
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
		
	// Record mean area, mean perimeter in a CSV file
	cell_population.AddPopulationWriter<CsvWriter>();
	// Record polygon geometry information in a CSV file
	cell_population.AddPopulationWriter<ShapeWriter>();
	// Record neighbour information in a CSV file
	cell_population.AddPopulationWriter<NumNeighboursWriter>();
	// Record Rho GTPase level in cell ID 1
	cell_population.AddPopulationWriter<OneCellGTPaseWriter>();
	// Generate XML file
	cell_population.AddCellWriter<XMLCellWriter>();
		
        OffLatticeSimulation<2> simulator(cell_population);
		

        simulator.SetOutputDirectory("50x50GTPAse_2500_0.2beta_medAdhesion_Random_G_scale_1point15_deformation100_surface_0");

        
	simulator.SetSamplingTimestepMultiple(200);
	simulator.SetDt(0.01);
        simulator.SetEndTime(2500.0);
        
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
		
	/* Couples cell target area to GTPase concentration */
	MAKE_PTR(ODEParameterAreaModifier<2>, p_ODE_modifier);
	simulator.AddSimulationModifier(p_ODE_modifier);

	/* Default values
	* Cell-Boundary Adhesion = 1
	* Cell-Cell Adhesion = 0.5
	* Deformation Energy = 100
	* Membrane Surface Energy = 10
	*/
	
	/* GTPase simulation values
	* Deformation Energy = 100
	* Membrane Surface Energy = 0 (no perimeter constraint)
	* High Adhesion Scenario: cell-cell: 0.75, cell-boundary: 1.0
	* Medium Adhesion Scenario: cell-cell: 1.0, cell-boundary: 1.0
	* Low Adhesion Scenario: cell-cell: 1.0, cell-boundary: 0.75
	*/
	MAKE_PTR(NagaiHondaForce<2>, p_force);
	p_force->SetNagaiHondaDeformationEnergyParameter(100.0);
	p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
	p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0);
	p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);

        simulator.Solve();

    }
};

#endif /*MULTICELLSNODIVISIONCOUPLED_HPP_*/
