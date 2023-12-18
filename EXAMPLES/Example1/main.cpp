/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.6 $
// $Date: 2009-11-02 22:23:58 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/Example1/main.cpp,v $


// File: ~/model/main.C
//
// Written: fmk 08/99
//
// Purpose: this file contains a C++ main procedure to perform the analysis
// of example1 (found in most documents). In the main() procedure:
// 	1) each object of the domain, i.e. Nodes, Elements, Constraints,
//	   and LoadPattern objects are created and then added to the Domain.
//	2) the components of the analysis object are constructed and then
//	   the Analysis object is created.
//	3) the analysis is performed.
//	4) the results are printed - here the contents of Domain and end of
//	   the analysis operation.

// standard C++ includes

//#include <stdlib.h>
//
//#include <OPS_Globals.h>
//#include <StandardStream.h>
//
//#include <ArrayOfTaggedObjects.h>
//
//// includes for the domain classes
//#include <Domain.h>
//#include <Node.h>
//#include <Truss.h>
//#include <ElasticMaterial.h>
//#include <SP_Constraint.h>
//#include <LoadPattern.h>
//#include <LinearSeries.h>
//#include <NodalLoad.h>
//
//// includes for the analysis classes
//#include <StaticAnalysis.h>
//#include <AnalysisModel.h>
//#include <Linear.h>
//#include <PenaltyConstraintHandler.h>
//#include <DOF_Numberer.h>
//#include <RCM.h>
//#include <LoadControl.h>
//#include <BandSPDLinSOE.h>
//#include <BandSPDLinLapackSolver.h>
//
//
//// init the global variabled defined in OPS_Globals.h
//StandardStream sserr;
//OPS_Stream *opserrPtr = &sserr;
//
//
//
//
//// main routine
//int main(int argc, char **argv)
//{
//    //
//    //	now create a domain and a modelbuilder
//    //  and build the model
//    //
//
//    Domain *theDomain = new Domain();
//    
//    // create the nodes using constructor: 
//    //		Node(tag, ndof, crd1, crd2)
//    // and then add them to the domain
//    
//    Node *node1 = new Node(1, 2,   0.0,  0.0);
//    Node *node2 = new Node(2, 2, 144.0,  0.0);
//    Node *node3 = new Node(3, 2, 168.0,  0.0);    
//    Node *node4 = new Node(4, 2,  72.0, 96.0);        
//    theDomain->addNode(node1);
//    theDomain->addNode(node2);
//    theDomain->addNode(node3);
//    theDomain->addNode(node4);
//    
//    // create an elastic material using constriuctor:  
//    //		ElasticMaterialModel(tag, E)
//
//    UniaxialMaterial *theMaterial = new ElasticMaterial(1, 3000);
//    
//    // create the truss elements using constructor:
//    //		Truss(tag, dim, nd1, nd2, Material &,A)
//    // and then add them to the domain
//    
//    Truss *truss1 = new Truss(1, 2, 1, 4, *theMaterial, 10.0);
//    Truss *truss2 = new Truss(2, 2, 2, 4, *theMaterial,  5.0);    
//    Truss *truss3 = new Truss(3, 2, 3, 4, *theMaterial,  5.0);        
//    theDomain->addElement(truss1);
//    theDomain->addElement(truss2);
//    theDomain->addElement(truss3);    
//    
//    // create the single-point constraint objects using constructor:
//    //		SP_Constraint(tag, nodeTag, dofID, value)
//    // and then add them to the domain
//    
//    SP_Constraint *sp1 = new SP_Constraint(1, 0, 0.0);
//    SP_Constraint *sp2 = new SP_Constraint(1, 1, 0.0);    
//    SP_Constraint *sp3 = new SP_Constraint(2, 0, 0.0);
//    SP_Constraint *sp4 = new SP_Constraint(2, 1, 0.0);    
//    SP_Constraint *sp5 = new SP_Constraint(3, 0, 0.0);
//    SP_Constraint *sp6 = new SP_Constraint(3, 1, 0.0);        
//    theDomain->addSP_Constraint(sp1);
//    theDomain->addSP_Constraint(sp2);
//    theDomain->addSP_Constraint(sp3);
//    theDomain->addSP_Constraint(sp4);    
//    theDomain->addSP_Constraint(sp5);    
//    theDomain->addSP_Constraint(sp6);    
//
//    // construct a linear time series object using constructor:
//    //		LinearSeries()
//    
//    TimeSeries *theSeries = new LinearSeries();
//    
//    // construct a load pattern using constructor:
//    //		LoadPattern(tag)
//    // and then set it's TimeSeries and add it to the domain
//    
//    LoadPattern *theLoadPattern = new LoadPattern(1);
//    theLoadPattern->setTimeSeries(theSeries);
//    theDomain->addLoadPattern(theLoadPattern);
//    
//    // construct a nodal load using constructor:
//    //		NodalLoad(tag, nodeID, Vector &)
//    // first construct a Vector of size 2 and set the values NOTE C INDEXING
//    // then construct the load and add it to the domain
//    
//    Vector theLoadValues(2);
//    theLoadValues(0) = 100.0;
//    theLoadValues(1) = -50.0;
//    NodalLoad *theLoad = new NodalLoad(1, 4, theLoadValues);
//    theDomain->addNodalLoad(theLoad, 1);
//
//    // create an Analysis object to perform a static analysis of the model
//    //  - constructs:
//    //    AnalysisModel of type AnalysisModel,
//    //	  EquiSolnAlgo of type Linear
//    //	  StaticIntegrator of type LoadControl
//    //	  ConstraintHandler of type Penalty
//    //    DOF_Numberer which uses RCM
//    //    LinearSOE of type Band SPD
//    // and then the StaticAnalysis object
//    
//    AnalysisModel     *theModel = new AnalysisModel();
//    EquiSolnAlgo      *theSolnAlgo = new Linear();
//    StaticIntegrator  *theIntegrator = new LoadControl(1.0, 1, 1.0, 1.0);
//    ConstraintHandler *theHandler = new PenaltyConstraintHandler(1.0e8,1.0e8);
//    RCM               *theRCM = new RCM();
//    DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);    
//    BandSPDLinSolver  *theSolver = new BandSPDLinLapackSolver();       
//    LinearSOE         *theSOE = new BandSPDLinSOE(*theSolver);        
//
//    StaticAnalysis    theAnalysis(*theDomain,
//				  *theHandler,
//				  *theNumberer,
//				  *theModel,
//				  *theSolnAlgo,
//				  *theSOE,
//				  *theIntegrator);
//
//    // perform the analysis & print out the results for the domain
//    int numSteps = 1;
//    theAnalysis.analyze(numSteps);
//    opserr << *theDomain;
//
//    // Clean up memory before exit
//    theAnalysis.clearAll();
//    theDomain->clearAll();
//    delete theDomain;
//    delete theMaterial;
//    
//    exit(0);
//}	
	
/* *********************************************************************************************** **
**				Test para revisar el funcionamiento los metodos agregados para MEFI3D		       **
** **********************************************************************************************+ */
#include "Vector.h"
#include "ID.h"
#include "Matrix.h"

#include <Concrete02.h>
#include <Steel02.h>
#include <OrthotropicRotatingAngleConcreteT2DMaterial01/OrthotropicRotatingAngleConcreteT2DMaterial01.h>
#include <SmearedSteelDoubleLayerT2DMaterial01/SmearedSteelDoubleLayerT2DMaterial01.h>
#include <ReinforcedConcreteLayerMembraneSection/ReinforcedConcreteLayerMembraneSection01.h>
#include <FSAM.h>
#include <ReinforcedConcreteLayerMembraneSection/ReinforcedConcreteLayerMembraneSection02.h>

#include <OPS_Globals.h>
#include <StandardStream.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

//init the global variabled defined in OPS_Globals.h
StandardStream sserr;
OPS_Stream* opserrPtr = &sserr;


int main()
{
	// Concrete material propierties --------------------------------------------------------------------------------------------------
	// Unconfined propierties
	double fpc = -47.09;
	double ec0 = -0.00232;
	double ft = 2.13;
	double et = 0.00008;
	double Ec = 34766.59;
	double fcu = 0.0 * fpc;
	double ecu = -0.037;
	double Et = 0.05 * Ec;
	// Confined propierties
	double fpcc = -53.78;
	double ec0c = -0.00397;
	double Ecc = 36542.37;
	double fcuc = 0.2 * fpc;
	double ecuc = -0.047;
	double Etc = 0.05 * Ecc;

	double magnitudGSelfWeightLoad = 9800.0;
	double rhoConcreteMaterial = 2500.0 * pow(10, -9) / magnitudGSelfWeightLoad;

	// Steel properties ---------------------------------------------------------------------------------------------------------------
	double fyX = 469.93;
	double bx = 0.02;
	double fyYw = 409.71;
	double byw = 0.02;
	double fyYb = 429.78;
	double byb = 0.01;
	double Esy = 200000.0;
	double Esx = Esy;
	double R0 = 20.0;
	double A1 = 0.925;
	double A2 = 0.15;
	double rouXb = 0.0027;
	double rouYb = 0.0323;

	// ========================== BUILD CONCRETE MATERIALS ============================================================================
	UniaxialMaterial* concreteUniaxialMat1 = new Concrete02(1, fpc, ec0, fcu, ecu, 0.1, ft, Et);
	UniaxialMaterial* concreteUniaxialMat2 = new Concrete02(2, fpcc, ec0c, fcuc, ecuc, 0.1, ft, Etc);
	NDMaterial** ORAC = new NDMaterial * [2];
	ORAC[0] = new OrthotropicRotatingAngleConcreteT2DMaterial01(3, concreteUniaxialMat1, concreteUniaxialMat1, et, ec0, rhoConcreteMaterial);
	ORAC[1] = new OrthotropicRotatingAngleConcreteT2DMaterial01(4, concreteUniaxialMat2, concreteUniaxialMat2, et, ec0c, rhoConcreteMaterial);

	// ========================== BUILD STEEL MATERIALS ================================================================================
	UniaxialMaterial* steelUniaxialMat1 = new Steel02(5, fyX, Esx, bx, R0, A1, A2);
	UniaxialMaterial* steelUniaxialMat2 = new Steel02(6, fyYb, Esy, byb, R0, A1, A2);
	NDMaterial** SmearedSteel = new NDMaterial * [1];
	SmearedSteel[0] = new SmearedSteelDoubleLayerT2DMaterial01(7, steelUniaxialMat1, steelUniaxialMat2, rouXb, rouYb, 0.0);

	// FSAM
	//double nu = 0.35;
	//double alfadow = 0.005;

	//NDMaterial* FSAMmaterial = new FSAM(8, rhoConcreteMaterial, steelUniaxialMat1, steelUniaxialMat2, concreteUniaxialMat2, concreteUniaxialMat2, concreteUniaxialMat2, concreteUniaxialMat2, concreteUniaxialMat2, concreteUniaxialMat2, rouXb, rouYb, nu, alfadow);

	// ========================== BUILD SECTION ========================================================================================
	double* t = new double[2];
	t[0] = 81.0;
	t[1] = 71.4;

	double h = 152.4;

	ReinforcedConcreteLayerMembraneSection01 RCMembraneSection(9, 1, 2, SmearedSteel, ORAC, t);
	//ReinforcedConcreteLayerMembraneSection02 RCMembraneSectionFSAM(10, FSAMmaterial, h);

	// Perfil de deformaciones --------------------------------------------------------------------------------------------------------
	// Vamos a leer un archivo de texto que contiene un vector de vectores con los perfiles de deformaciones 
	string filename = "strainConcreteMatrix.txt";

	ifstream archivo_entrada(filename);			// Se abre el archivo de entrada en modo lectura
	if (archivo_entrada.is_open()) {
		vector<vector<double>>matriz;
		string linea;

		// Se lee cada linea del archivo
		while (getline(archivo_entrada, linea)) {
			vector<double>fila;
			stringstream ss(linea);
			double valor;
			// Se lee los valores de la línea y se almacena en la fila
			while (ss >> valor) {
				fila.push_back(valor);
			}
			// Se agrega la fila a la matriz
			matriz.push_back(fila);
		}

		// Se trabaja ahora con la matriz de valores
		for (int i = 0; i < matriz.size(); i++) {
			// Se obtiene la fila actual
			vector<double>& fila = matriz[i];

			// Se guarda cada fila en un objeto de tipo "Vector"
			Vector strain(3);
			for (int j = 0; j < fila.size(); j++) {
				strain(j) = fila[j];
			}
			// Se aplica el perfil de deformaciones y se obtienen tensiones y matrices tangentes
			opserr << "\n ********** SET TRIAL STRAIN **********" << endln;
			SmearedSteel[0]->setTrialStrain(strain);
			ORAC[0]->setTrialStrain(strain);
			RCMembraneSection.setTrialSectionDeformation(strain);
			opserr << "\n ********** GET STRAIN **********" << endln;
			const Vector& epsConcND = ORAC[0]->getStrain();
			const Vector& epsSteelND = SmearedSteel[0]->getStrain();
			const Vector& epsRCsection = RCMembraneSection.getSectionDeformation();
			opserr << "epsConcND = " << epsConcND << endln;
			opserr << "epsSteelND = " << epsSteelND << endln;
			opserr << "epsRCsection = " << epsRCsection << endln;
			opserr << "\n ********** GET STRESS **********" << endln;
			const Vector& sigConcND = ORAC[0]->getStress();
			const Vector& sigSteelND = SmearedSteel[0]->getStress();
			const Vector& sigRCsection = RCMembraneSection.getStressResultant();
			opserr << "sigConcND = " << sigConcND << endln;
			opserr << "sigSteelND = " << sigSteelND << endln;
			opserr << "sigRCsection = " << sigRCsection << endln;
			opserr << "\n ********** GET TANGENT **********" << endln;
			const Matrix& E_ConcND = ORAC[0]->getTangent();
			const Matrix& E_SteelND = SmearedSteel[0]->getTangent();
			const Matrix& E_RCsection = RCMembraneSection.getSectionTangent();
			opserr << "E_ConcND = " << E_ConcND << endln;
			opserr << "E_SteelND = " << E_SteelND << endln;
			opserr << "E_RCsection = " << E_RCsection << endln;
			opserr << "\n ********** GET INITIAL TANGENT **********" << endln;
			const Matrix& E0_ConcND = ORAC[0]->getInitialTangent();
			const Matrix& E0_SteelND = SmearedSteel[0]->getInitialTangent();
			const Matrix& E0_RCsection = RCMembraneSection.getInitialTangent();
			opserr << "E0_ConcND = " << E0_ConcND << endln;
			opserr << "E0_SteelND = " << E0_SteelND << endln;
			opserr << "E0_RCsection = " << E0_RCsection << endln;
			opserr << "\n ********** COMMIT STATE **********" << endln;
			SmearedSteel[0]->commitState();
			ORAC[0]->commitState();
			RCMembraneSection.commitState();
			// ======================= PRUEBA ===================================
			opserr << "\n *********** PRUEBA **************" << endln;
			double E0_CU = concreteUniaxialMat1->getInitialTangent();
			double E0_CC = concreteUniaxialMat2->getInitialTangent();

			double rhoRCsection = RCMembraneSection.getRho();
			/*double thicknessRCsection = RCMembraneSection.getSectionThickness();
			double EcAvgRCsection = RCMembraneSection.getEcAvg();
			Vector IP = RCMembraneSection.getInputParameters();*/

			//double rhoRCsectionFSAM = RCMembraneSectionFSAM.getRho();
			/*double thicknessRCsectionFSAM = RCMembraneSectionFSAM.getSectionThickness();
			double EcAvgRCsectionFSAM = RCMembraneSectionFSAM.getEcAvg();*/

			opserr << "E0_CU = " << E0_CU << endln;
			opserr << "E0_CC = " << E0_CC << endln;

			opserr << "rhoRCsection = " << rhoRCsection << endln;
			/*opserr << "thicknessRCsection = " << thicknessRCsection << endln;
			opserr << "EcAvgRCsection = " << EcAvgRCsection << endln;

			opserr << "rhoRCsectionFSAM = " << rhoRCsectionFSAM << endln;
			opserr << "thicknessRCsectionFSAM = " << thicknessRCsectionFSAM << endln;
			opserr << "EcAvgRCsectionFSAM = " << EcAvgRCsectionFSAM << endln;*/

			int flag = 1;
			RCMembraneSection.Print(sserr,flag);
			// ===================== FIN PRUEBA =================================

		}

		// Se cierra el archivo
		archivo_entrada.close();
	}
	else {
		cout << "No se pudo abrir el archivo." << endl;
	}

	delete concreteUniaxialMat1;
	delete concreteUniaxialMat2;
	delete steelUniaxialMat1;
	delete steelUniaxialMat2;
	//delete FSAMmaterial;
	delete[] ORAC;
	delete[] SmearedSteel;

	return 0;
}