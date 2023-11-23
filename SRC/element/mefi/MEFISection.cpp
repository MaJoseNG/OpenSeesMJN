// Code written/implemented by:	Carlos López Olea (carlos.lopez.o@ug.uchile.cl)
//								Fabian Rojas B.
//								Leonardo M. Massone
//
// User documentation available at: https://github.com/carloslopezolea/MEFI
//
// Created: 06/2022
//
// Adapted by: Maria Jose Nunez G.
//             05/2023 
// 
// Description: The Membrane Fiber (MEFI) element, is described by four nodes, each containing three degrees of freedom (DOFs), two translations, and one in-plane rotation 
// (drilling) DOF, which incorporates a blended interpolation function for the displacements over the element (Fig. 1b-c). The element formulation accommodates the quadrature 
// points and weights of the classical finite element formulation of membrane elements to resemble strips (fibers), similarly to macroscopic elements.
//
// Reference:
// 1.- López, C. N., Rojas, F., & Massone, L. M. (2022). Membrane fiber element for reinforced concrete walls – the benefits of macro and micro modeling approaches. Engineering Structures, 254, 113819.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mefi/MEFISection.cpp
//
// Rev: 1.0                                                       

#include <MEFISection.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>

void* OPS_MEFISection()
{
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 3) {
		opserr << "WARNING model dimensions and/or nodal DOF not compatible with MEFISection element\n";
	return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 8) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element MEFISection eleTag iNode jNode kNode lNode m -width -section\n";
		return 0;
    }

	int iData[6];

	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for element MEFISection" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
		opserr << "WARNING iNode jNode kNode lNode or m for element MEFISection" << iData[0] << endln;
		return 0;
	}

	int numFibers = iData[5];
	const char* str = 0;

	double* theWidth = new double[numFibers];
	int* sectionTags = new int[numFibers];

	SectionForceDeformation** theSection = new SectionForceDeformation * [numFibers];

	int numArgs = OPS_GetNumRemainingInputArgs();
	while (numArgs > 0) {
		str = OPS_GetString();
		if (strcmp(str, "-width") == 0) {
			numData = numFibers;
			if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
				opserr << "Invalid width value for MEFISection  " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-section") == 0) {
			numData = numFibers;
			if (OPS_GetIntInput(&numData, sectionTags) != 0) {
				opserr << "Invalid section tags for MEFISection  " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < numFibers; i++) {
				theSection[i] = 0;
				theSection[i] = OPS_getSectionForceDeformation(sectionTags[i]);
				if (theSection[i] == 0) {
					opserr << "Invalid section tag " << sectionTags[i] << "  for MEFISection  " << iData[0] << endln;
					return 0;
				}
			}
		}

		numArgs = OPS_GetNumRemainingInputArgs();

	}
	
    return new MEFISection(iData[0],iData[1],iData[2],iData[3],iData[4], iData[5], theSection, theWidth);
}

MEFISection::MEFISection(int tag, 
	int nd1, int nd2, int nd3, int nd4, 
	int numFibers,
	SectionForceDeformation** section,
	double* width)

	:Element (tag, ELE_TAG_MEFISection), 
	connectedExternalNodes(4), nd1Crds(2), nd2Crds(2), nd3Crds(2), nd4Crds(2),
	theSection(0), rho(0), NodeMass(0),
	MEFI_K(12, 12), MEFI_R(12), MEFI_D(12,12), MEFI_M(12,12),
	h(0), lw(0), x(numFibers), b(numFibers), nFibers(numFibers),
	BSS(3, 12), detJacobian(0.0), qdtLocations(numFibers, 2), qdtWeights(numFibers), qdtWeight(numFibers)
	
{

	// Check material width input
	if (width == 0) {
		opserr << "MEFISection::MEFISection() - Null width array passed.\n";
		exit(-1);
	}

	// Calculate locations of concrete macro-fibers in the cross-section (centerline - x = 0.0)

	for (int i = 0; i < nFibers; i++) {
		b(i) = width[i];
		lw += b(i);		// Total length of the wall
	}

	for (int i = 0; i < nFibers; i++)
		x(i) = 0.0;

	for (int i = 0; i < nFibers; i++) {
		double sumb_i = 0.0;
		for (int j = 0; j < i + 1; j++)
			sumb_i += b(j);

		x(i) = (sumb_i - b(i) / 2.0) - lw / 2.0;
	}

	// Calculate locations of quadrature points and weights
	for (int i = 0; i < nFibers; i++) {
		qdtLocations(i, 0) = x(i) * (2 / lw);
		qdtLocations(i, 1) = 0.001;
		qdtWeights(i) = b(i) * (2 / lw) * 2;
		// opserr << quadratureLocations(i, 0) << "--"<< quadratureLocations(i, 1) << "--" << quadratureWeights(i) << "\n";
	}

    // Allocate arrays of pointers to sections
	theSection = new SectionForceDeformation * [nFibers];
    
    if (theSection == 0) {
      opserr << "MEFISection::MEFISection() - failed allocate section model pointer\n";
      exit(-1);
    }

	// Get copies of the sections
	for (int i = 0; i < nFibers; i++) {
		if (section[i] == 0) {
			opserr << i;
			opserr << "MEFISection::MEFISection() - Null section pointer passed\n";
			exit(-1);
		}

		theSection[i] = section[i]->getCopy();

		if (theSection[i] == 0) {
			opserr << "MEFISection::MEFISection() - Failed to copy section\n";
			exit(-1);
		}
	}

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;

    for (int i=0; i<4; i++)
      theNodes[i] = 0;

}

MEFISection::MEFISection()
	:Element(0, ELE_TAG_MEFISection),
	connectedExternalNodes(4), nd1Crds(2), nd2Crds(2), nd3Crds(2), nd4Crds(2),
	theSection(0), rho(0), NodeMass(0),
	MEFI_K(12, 12), MEFI_R(12), MEFI_D(12, 12), MEFI_M(12, 12),
	h(0), lw(0), x(0), b(0), nFibers(0),
	BSS(3, 12), detJacobian(0.0), qdtLocations(1, 2), qdtWeights(1), qdtWeight(1)
{

    for (int i=0; i<4; i++)
      theNodes[i] = 0;
}

MEFISection::~MEFISection()
{    
  for (int i = 0; i < nFibers; i++) {
    if (theSection[i])
      delete theSection[i];
  }

  // Delete the array of pointers to Section pointer arrays
  if (theSection)
    delete [] theSection;

}

int
MEFISection::getNumExternalNodes() const
{
    return 4;
}

const ID&
MEFISection::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
MEFISection::getNodePtrs(void)
{
  return theNodes;
}

int
MEFISection::getNumDOF()
{
    return 12;
}

void
MEFISection::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
	//opserr << "FATAL ERROR FourNodeQuad (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3) {
		opserr << "MEFISection::setDomain(): 3 dof required at all nodes. " << dofNd1 << " provided at node 1, " << dofNd2 << " provided at node 2, "
			<< dofNd3 << " provided at node 4, " << dofNd4 << " provided at node 3";
	
		return;
    }

	// Get coordinates of end nodes
	nd1Crds = theNodes[0]->getCrds();
	nd2Crds = theNodes[1]->getCrds();
	nd3Crds = theNodes[2]->getCrds();
	nd4Crds = theNodes[3]->getCrds();

	// Calculate locations of quadrature points and weights
	for (int i = 0; i < nFibers; i++) {
		this->shapeFunction(qdtLocations(i, 0), qdtLocations(i, 1));
		qdtWeight(i) = qdtWeights(i) * detJacobian;
	}


    this->DomainComponent::setDomain(theDomain);
}

int
MEFISection::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "MEFISection::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the Section states
    for (int i = 0; i < nFibers; i++)
      retVal += theSection[i]->commitState();

    return retVal;
}

int
MEFISection::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < nFibers; i++)
		retVal += theSection[i]->revertToLastCommit();

    return retVal;
}

int
MEFISection::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < nFibers; i++)
		retVal += theSection[i]->revertToStart();

    return retVal;
}


int
MEFISection::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	// in-plane displacements for membrane behavior
	Vector dispMembrane(12); 
	dispMembrane.Zero();
	for (int i = 0; i < 3; i++) {
		dispMembrane(i) = disp1(i);
		dispMembrane(i + 3) = disp2(i);
		dispMembrane(i + 6) = disp3(i);
		dispMembrane(i + 9) = disp4(i);
	}

	// Strains at each fiber
	Vector strainAtQuadraturePoint(3);
	strainAtQuadraturePoint.Zero();

	int ret = 0;
	// Loop over the integration points
	for (int i = 0; i < nFibers; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(qdtLocations(i,0), qdtLocations(i, 1));

		// Interpolate strains
		strainAtQuadraturePoint = BSS * dispMembrane;

		// Set the material strain
		ret += theSection[i]->setTrialSectionDeformation(strainAtQuadraturePoint);
	}

	return ret;
}

const Matrix&
MEFISection::getTangentStiff()
{
	MEFI_K.Zero();
	Matrix Ki(12, 12); 
	Matrix KiSection(3, 3); 

	// Loop over the integration points
	for (int i = 0; i < nFibers; i++) {
		Ki.Zero(); KiSection.Zero();
		const Matrix& D = theSection[i]->getSectionTangent();
		this->shapeFunction(qdtLocations(i, 0), qdtLocations(i, 1));
		KiSection = D;
		Ki.addMatrixTripleProduct(0.0, BSS, KiSection, 1.0);
		MEFI_K = MEFI_K + Ki * qdtWeight(i);
	}
	
	return MEFI_K;
}


const Matrix&
MEFISection::getInitialStiff()
{
	MEFI_K.Zero();
	Matrix Ki(12, 12);
	Matrix KiSection(3, 3);

	// Loop over the integration points
	for (int i = 0; i < nFibers; i++) {
		Ki.Zero(); KiSection.Zero();
		const Matrix& D = theSection[i]->getInitialTangent();
		this->shapeFunction(qdtLocations(i, 0), qdtLocations(i, 1));
		KiSection = D;
		Ki.addMatrixTripleProduct(0.0, BSS, KiSection, 1.0);
		MEFI_K = MEFI_K + Ki * qdtWeight(i);
	}
	return MEFI_K;
}

const Matrix&
MEFISection::getMass()
{
	MEFI_M.Zero();
	return MEFI_M;
}

// N/A to this model - no element loads
void MEFISection::zeroLoad(void)
{
  	return;
}

// N/A to this model - no element loads
int MEFISection::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int 
MEFISection::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector&
MEFISection::getResistingForce()
{
	MEFI_R.Zero();

	// Get Trial Displacements
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector& disp3 = theNodes[2]->getTrialDisp();
	const Vector& disp4 = theNodes[3]->getTrialDisp();

	// in-plane displacements for membrane behavior
	Vector dispMembrane(12);
	dispMembrane.Zero();
	for (int i = 0; i < 3; i++) {
		dispMembrane(i) = disp1(i);
		dispMembrane(i + 3) = disp2(i);
		dispMembrane(i + 6) = disp3(i);
		dispMembrane(i + 9) = disp4(i);
	}

	//opserr << dispMembrane(0) << dispMembrane(1) << dispMembrane(2) << dispMembrane(3) << dispMembrane(4) << dispMembrane(5) << dispMembrane(6) << dispMembrane(7) << dispMembrane(8) << dispMembrane(9) << dispMembrane(10) << dispMembrane(11) << endln;


	Matrix BSST(12, 3);
	// Loop over the integration points
	for (int i = 0; i < nFibers; i++) {
		BSST.Zero();
		const Vector& Stress = theSection[i]->getStressResultant();
		this->shapeFunction(qdtLocations(i, 0), qdtLocations(i, 1));
		BSST = this->transpose(BSS);
		MEFI_R = MEFI_R + BSST * Stress * qdtWeight(i);
		
	}

	//opserr << MEFI_R(0) << MEFI_R(1) << MEFI_R(2) << MEFI_R(3) << MEFI_R(4) << MEFI_R(5) << MEFI_R(6) << MEFI_R(7) << MEFI_R(8) << MEFI_R(9) << MEFI_R(10) << MEFI_R(11) << endln;

	return MEFI_R;
}

const Vector&
MEFISection::getResistingForceIncInertia()
{
	MEFI_R.Zero();
	return MEFI_R;
}

int MEFISection::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int MEFISection::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
	return -1;
}

void MEFISection::Print(OPS_Stream &s, int flag)
{
  
}

int MEFISection::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	return -1;
}

Response* MEFISection::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  // Added by MJN
  output.tag("ElementOutput");
  output.attr("eleType", "MEFISection");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);
  output.attr("node3", connectedExternalNodes[2]);
  output.attr("node4", connectedExternalNodes[3]);

  // Material output
  if (strcmp(argv[0], "RCpanel") == 0 || strcmp(argv[0], "RCPanel")
	  || strcmp(argv[0], "RC_panel") || strcmp(argv[0], "RC_Panel") == 0)
  {
	  //Check if correct # of arguments passed
	  if (argc != 3) {
		  opserr << "WARNING: Number of recorder input for RC Panel is: " << argc - 1 << "; should be 3: panTag (one panel only: 1 to m) and $Response_Type.\n";
		  return 0;
	  }

	  int secNum = atoi(argv[1]);

	  output.tag("Material");
	  output.attr("number", secNum);

	  return theResponse = theSection[secNum - 1]->setResponse(&argv[argc - 1], argc - 2, output);
  }

  output.endTag();

  return 0;
}

int  MEFISection::getResponse(int responseID, Information &eleInfo)
{

    return 0;
}

int MEFISection::setParameter(const char **argv, int argc, Parameter &param)
{ 
  return -1;
}
    
int MEFISection::updateParameter(int parameterID, Information &info)
{	
	return -1;
}

void MEFISection::shapeFunction(double xi, double eta)
{

	double x1 = nd1Crds(0); double y1 = nd1Crds(1);
	double x2 = nd2Crds(0); double y2 = nd2Crds(1);
	double x3 = nd3Crds(0); double y3 = nd3Crds(1); 
	double x4 = nd4Crds(0); double y4 = nd4Crds(1);

	// dShapeFunction Matrix
	Matrix dShapeFunction(2, 8); dShapeFunction.Zero();
	dShapeFunction(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
	dShapeFunction(0, 1) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
	dShapeFunction(0, 2) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
	dShapeFunction(0, 3) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
	dShapeFunction(0, 4) = 0.5 * (-2.0 * xi) * (1.0 - eta);
	dShapeFunction(0, 5) = 0.5 * (1.0 - pow(eta, 2));
	dShapeFunction(0, 6) = 0.5 * (-2.0 * xi) * (1.0 + eta);
	dShapeFunction(0, 7) = 0.5 * (-1.0) * (1.0 - pow(eta, 2));
	dShapeFunction(1, 0) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
	dShapeFunction(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
	dShapeFunction(1, 2) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
	dShapeFunction(1, 3) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
	dShapeFunction(1, 4) = 0.5 * (-1.0) * (1.0 - pow(xi, 2));
	dShapeFunction(1, 5) = 0.5 * (-2.0 * eta) * (1.0 + xi);
	dShapeFunction(1, 6) = 0.5 * (1.0 - pow(xi, 2));
	dShapeFunction(1, 7) = 0.5 * (-2.0 * eta) * (1.0 - xi);

	// localCoord8Nodes Matrix
	Matrix localCoord8Nodes(8, 2); localCoord8Nodes.Zero();
	localCoord8Nodes(0, 0) = x1;              localCoord8Nodes(0, 1) = y1;
	localCoord8Nodes(1, 0) = x2;              localCoord8Nodes(1, 1) = y2;
	localCoord8Nodes(2, 0) = x3;              localCoord8Nodes(2, 1) = y3;
	localCoord8Nodes(3, 0) = x4;              localCoord8Nodes(3, 1) = y4;
	localCoord8Nodes(4, 0) = 0.5 * (x1 + x2); localCoord8Nodes(4, 1) = 0.5 * (y1 + y2);
	localCoord8Nodes(5, 0) = 0.5 * (x2 + x3); localCoord8Nodes(5, 1) = 0.5 * (y2 + y3);
	localCoord8Nodes(6, 0) = 0.5 * (x3 + x4); localCoord8Nodes(6, 1) = 0.5 * (y3 + y4);
	localCoord8Nodes(7, 0) = 0.5 * (x4 + x1); localCoord8Nodes(7, 1) = 0.5 * (y4 + y1);

	// Jacobian Matrix
	Matrix JacobianMatrix(2, 2); JacobianMatrix.Zero();
	JacobianMatrix = dShapeFunction * localCoord8Nodes;

	// Jacobian determinant
	detJacobian = JacobianMatrix(0, 0) * JacobianMatrix(1, 1) - JacobianMatrix(0, 1) * JacobianMatrix(1, 0);

	// Jacobian Matrix Inverse
	Matrix inverseJacobianMatrix(2, 2); inverseJacobianMatrix.Zero();
	inverseJacobianMatrix(0, 0) = JacobianMatrix(1, 1) / detJacobian;
	inverseJacobianMatrix(1, 0) = -JacobianMatrix(0, 1) / detJacobian;
	inverseJacobianMatrix(0, 1) = -JacobianMatrix(1, 0) / detJacobian;
	inverseJacobianMatrix(1, 1) = JacobianMatrix(0, 0) / detJacobian;

	// J Matrix
	Matrix JMatrix(4, 4); JMatrix.Zero();
	JMatrix(0, 0) = inverseJacobianMatrix(0, 0);
	JMatrix(0, 1) = inverseJacobianMatrix(0, 1);
	JMatrix(1, 0) = inverseJacobianMatrix(1, 0);
	JMatrix(1, 1) = inverseJacobianMatrix(1, 1);
	JMatrix(2, 2) = inverseJacobianMatrix(0, 0);
	JMatrix(2, 3) = inverseJacobianMatrix(0, 1);
	JMatrix(3, 2) = inverseJacobianMatrix(1, 0);
	JMatrix(3, 3) = inverseJacobianMatrix(1, 1);

	// C Matrix
	Matrix CMatrix(4, 16); CMatrix.Zero();
	CMatrix(0, 0) = (0.5) * (-(0.5) + (0.75) * (eta)-0.25 * (pow(eta, 3)));
	CMatrix(0, 1) = 0.0;
	CMatrix(0, 2) = (0.5) * (0.25 - 0.25 * eta - 0.25 * (pow(eta, 2)) + 0.25 * (pow(eta, 3)));
	CMatrix(0, 3) = 0.0;
	CMatrix(0, 4) = (0.5) * ((0.5) - (0.75) * (eta)+0.25 * (pow(eta, 3)));
	CMatrix(0, 5) = 0.0;
	CMatrix(0, 6) = (0.5) * (-(0.25) + 0.25 * eta + 0.25 * (pow(eta, 2)) - 0.25 * (pow(eta, 3)));
	CMatrix(0, 7) = 0.0;
	CMatrix(0, 8) = (0.5) * ((0.5) + (0.75) * (eta)-0.25 * (pow(eta, 3)));
	CMatrix(0, 9) = 0.0;
	CMatrix(0, 10) = (0.5) * (0.25 + 0.25 * eta - 0.25 * (pow(eta, 2)) - 0.25 * (pow(eta, 3)));
	CMatrix(0, 11) = 0.0;
	CMatrix(0, 12) = (0.5) * (-(0.5) - (0.75) * (eta)+0.25 * (pow(eta, 3)));
	CMatrix(0, 13) = 0.0;
	CMatrix(0, 14) = (0.5) * (-(0.25) - 0.25 * eta + 0.25 * (pow(eta, 2)) + 0.25 * (pow(eta, 3)));
	CMatrix(0, 15) = 0.0;

	CMatrix(1, 0) = (0.5) * (-(0.75) + (0.75) * (pow(eta, 2))) * (1.0 - xi);
	CMatrix(1, 1) = 0.0;
	CMatrix(1, 2) = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1.0 + xi);
	CMatrix(1, 3) = 0.0;
	CMatrix(1, 4) = (0.5) * (-(0.75) + (0.75) * (pow(eta, 2))) * (1.0 + xi);
	CMatrix(1, 5) = 0.0;
	CMatrix(1, 6) = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1.0 - xi);
	CMatrix(1, 7) = 0.0;
	CMatrix(1, 8) = (0.5) * ((0.75) - (0.75) * (pow(eta, 2))) * (1 + xi);
	CMatrix(1, 9) = 0.0;
	CMatrix(1, 10) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1 - xi);
	CMatrix(1, 11) = 0.0;
	CMatrix(1, 12) = (0.5) * ((0.75) - (0.75) * (pow(eta, 2))) * (1 - xi);
	CMatrix(1, 13) = 0.0;
	CMatrix(1, 14) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1 + xi);
	CMatrix(1, 15) = 0.0;

	CMatrix(2, 0) = 0.0;
	CMatrix(2, 1) = (0.5) * (1 - eta) * (-(0.75) + (0.75) * (pow(xi, 2)));
	CMatrix(2, 2) = 0.0;
	CMatrix(2, 3) = (0.5) * (1 - eta) * (-(0.25) - 0.5 * xi + (0.75) * (pow(xi, 2)));
	CMatrix(2, 4) = 0.0;
	CMatrix(2, 5) = (0.5) * (1 - eta) * ((0.75) - (0.75) * (pow(xi, 2)));
	CMatrix(2, 6) = 0.0;
	CMatrix(2, 7) = (0.5) * (1 - eta) * (-(0.25) + 0.5 * xi + (0.75) * (pow(xi, 2)));
	CMatrix(2, 8) = 0.0;
	CMatrix(2, 9) = (0.5) * (1 + eta) * ((0.75) - (0.75) * (pow(xi, 2)));
	CMatrix(2, 10) = 0.0;
	CMatrix(2, 11) = (0.5) * (1 + eta) * (-(0.25) + 0.5 * xi + (0.75) * (pow(xi, 2)));
	CMatrix(2, 12) = 0.0;
	CMatrix(2, 13) = (0.5) * (1 + eta) * (-(0.75) + (0.75) * (pow(xi, 2)));
	CMatrix(2, 14) = 0.0;
	CMatrix(2, 15) = (0.5) * (1 + eta) * (-0.25 - 0.5 * xi + (0.75) * (pow(xi, 2)));

	CMatrix(3, 0) = 0.0;
	CMatrix(3, 1) = (0.5) * (-0.5 + (0.75) * (xi)-(0.25) * (pow(xi, 3)));
	CMatrix(3, 2) = 0.0;
	CMatrix(3, 3) = (0.5) * (-(0.25) + (0.25) * xi + (0.25) * (pow(xi, 2)) - (0.25) * (pow(xi, 3)));
	CMatrix(3, 4) = 0.0;
	CMatrix(3, 5) = (0.5) * (-(0.5) - (0.75) * (xi)+(0.25) * (pow(xi, 3)));
	CMatrix(3, 6) = 0.0;
	CMatrix(3, 7) = (0.5) * (0.25 + 0.25 * xi - (0.25) * (pow(xi, 2)) - (0.25) * (pow(xi, 3)));
	CMatrix(3, 8) = 0.0;
	CMatrix(3, 9) = (0.5) * ((0.5) + (0.75) * (xi)-(0.25) * (pow(xi, 3)));
	CMatrix(3, 10) = 0.0;
	CMatrix(3, 11) = (0.5) * (-(0.25) - (0.25) * xi + (0.25) * (pow(xi, 2)) + (0.25) * (pow(xi, 3)));
	CMatrix(3, 12) = 0.0;
	CMatrix(3, 13) = (0.5) * ((0.5) - (0.75) * (xi)+(0.25) * (pow(xi, 3)));
	CMatrix(3, 14) = 0.0;
	CMatrix(3, 15) = (0.5) * ((0.25) - (0.25) * xi - (0.25) * (pow(xi, 2)) + (0.25) * (pow(xi, 3)));
	
	// Tr Matrix
	double x21 = 0.5 * (x2 - x1);
	double x34 = 0.5 * (x3 - x4);
	double y41 = 0.5 * (y4 - y1);
	double y32 = 0.5 * (y3 - y2);

	Matrix TrMatrix(16, 12); TrMatrix.Zero();
	TrMatrix(0, 0) = 1.0;
	TrMatrix(1, 1) = 1.0;
	TrMatrix(2, 2) = y41;
	TrMatrix(3, 2) = x21;
	TrMatrix(4, 3) = 1.0;
	TrMatrix(5, 4) = 1.0;
	TrMatrix(6, 5) = y32;
	TrMatrix(7, 5) = x21;
	TrMatrix(8, 6) = 1.0;
	TrMatrix(9, 7) = 1.0;
	TrMatrix(10, 8) = y32;
	TrMatrix(11, 8) = x34;
	TrMatrix(12, 9) = 1.0;
	TrMatrix(13, 10) = 1.0;
	TrMatrix(14, 11) = y41;
	TrMatrix(15, 11) = x34;

	// Bss Matrix
	Matrix AMatrix(3, 4); AMatrix.Zero();
	AMatrix(0, 0) = 1.0;
	AMatrix(1, 3) = 1.0;
	AMatrix(2, 1) = 1.0;
	AMatrix(2, 2) = 1.0;

	BSS.Zero();
	BSS = AMatrix * (JMatrix * CMatrix * TrMatrix);
}


const Matrix&
MEFISection::transpose(const Matrix& M)
{
	int i;
	int j;

	//we're always transposing 3x12 matrices for this element,
	//so always return a 12x3 .

	static int dim1 = 12;
	static int dim2 = 3;
	static Matrix Mtran(dim1, dim2);

	for (i = 0; i < dim1; i++) {
		for (j = 0; j < dim2; j++)
			Mtran(i, j) = M(j, i);
	} // end for i

	return Mtran;
}