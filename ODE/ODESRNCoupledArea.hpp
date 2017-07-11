#ifndef ODESRN_HPP_
#define ODESRN_HPP_

/*
 * Rho GTPase Simulation
 * Author: MoHan Zhang <mohan_z@hotmail.com>
 * Last Modified: July 11, 2017
 * Do not reproduce this code without permission.
 */

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"

/* The next header includes the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/* The next header defines a base class for ode-based SRN models.
 * Our new SRN model will inherit from this abstract class. */
#include "AbstractOdeSrnModel.hpp"

/* These headers specify the methods to solve the ODE system.*/
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/* This header specifies the ODE solvers.*/
#include "CellCycleModelOdeSolver.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based Chaste
 * tutorials. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

#include <cmath>

class ODESRN : public AbstractOdeSystem
{
	
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:
    ODESRN() : AbstractOdeSystem(3)
    {
        mpSystemInfo = OdeSystemInformation<ODESRN>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
		// Birfurcation parameter
		double beta = 0.2;
		// GTPase eqn
        rDY[0] = ((0.1 + beta * (pow(rY[2],10)/(pow(rY[1],10)+pow(rY[2],10))) + 1.5*(pow(rY[0],4)/(1+pow(rY[0],4))))*(2-rY[0])-rY[0])*0.25;
		// Target Area eqn
        rDY[1] = (-0.1*(rY[1]-1.15*(1-0.75*(pow(rY[0],4)/(pow(0.3,4)+pow(rY[0],4))))))*0.25;
		// Dummy eqn to get access to cell area, initially equal to 0.866025
		rDY[2] = 0;
    }
};

template<>
void OdeSystemInformation<ODESRN>::Initialise()
{
    this->mVariableNames.push_back("G");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1);

    this->mVariableNames.push_back("TARGET AREA");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.8);
	
    this->mVariableNames.push_back("AREA");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.8);

    this->mInitialised = true;
}


class ODESrnModel : public AbstractOdeSrnModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

public:

    ODESrnModel()
        : AbstractOdeSrnModel(3, boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
		// ODE solver
        mpOdeSolver = CellCycleModelOdeSolver<ODESrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.01);

        assert(mpOdeSolver->IsSetUp());
    }

    AbstractSrnModel* CreateSrnModel()
    {
        ODESrnModel* p_model = new ODESrnModel();

        p_model->SetOdeSystem(new ODESRN);

        return AbstractOdeSrnModel::CreateSrnModel(p_model);
    }

    void Initialise()
    {
        AbstractOdeSrnModel::Initialise(new ODESRN);
    }

    void SimulateToCurrentTime()
    {
        // run the ODE simulation as needed
        AbstractOdeSrnModel::SimulateToCurrentTime();

        /* Output the ODE system variable to {{{CellData}}}. */
        mpCell->GetCellData()->SetItem("G",mpOdeSystem->rGetStateVariables()[0]);
		mpCell->GetCellData()->SetItem("target area",mpOdeSystem->rGetStateVariables()[1]);
		mpCell->GetCellData()->SetItem("AREA",mpOdeSystem->rGetStateVariables()[2]);
    }
	
	void ResetForDivision(){
		AbstractOdeSrnModel::ResetForDivision();
	    std::vector<double> init_conds = mpOdeSystem->GetInitialConditions();
	    for (unsigned i=0; i<2; i++)
	    {
	        mpOdeSystem->rGetStateVariables()[i] = init_conds[i];
	    }
	}

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ODESRN)
CHASTE_CLASS_EXPORT(ODESrnModel)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ODESRN)
CHASTE_CLASS_EXPORT(ODESrnModel)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(ODESrnModel)

#endif /*ODESRN_HPP_*/