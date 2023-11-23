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
// Source: /usr/local/cvs/OpenSees/SRC/element/mefi/MEFISection.h
//
// Rev: 1.0         

#ifndef MEFISection_h
#define MEFISection_h

#ifndef _bool_h
#include "bool.h"
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class SectionForceDeformation;
class Response;

class MEFISection : public Element
{
  public:
      MEFISection(int tag,                          // element tag
          int nd1, int nd2, int nd3, int nd4,       // node tags
          int numFibers,                            // number of fiber elements
          SectionForceDeformation** Sections,		// array of section tags
          double* Width);							// array of macro-fiber widths;
          
      MEFISection();
      ~MEFISection();

    const char *getClassType(void) const {return "FourNodeQuad";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getMass(void);    

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);

    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);


  protected:
    
  private:
    // private attributes - a copy for each object of the class
    Node* theNodes[4];                      // external node pointers 
    SectionForceDeformation **theSection;   // pointer to the section objects   
    ID connectedExternalNodes;              // Tags of nodes

    
    //static double matrixData[64];  // array data for matrix
    //static Matrix K;		// Element stiffness, damping, and mass Matrix
    //static Vector P;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads
	
    double rho;


    // private member functions - only objects of this class can call these
    void shapeFunction(double xi, double eta);


    // Nodal coordinates
    Vector nd1Crds;
    Vector nd2Crds;
    Vector nd3Crds;
    Vector nd4Crds;

    // calculated element parameters
    double NodeMass;					// nodal mass
    double nFibers;					    // number of fibers of MEFI element
    double h;							// initial height of MEFI element
    double lw;							// initial length of MEFI element
    Vector x;				        	// fiber locations
    Vector b;	        				// fiber widths

    Matrix BSS;
    double detJacobian;
    Matrix qdtLocations;
    Vector qdtWeights;
    Vector qdtWeight;

    Matrix MEFI_K;		// stiffness
    Matrix MEFI_D;		// damping
    Matrix MEFI_M;		// mass 
    Vector MEFI_R;		// force

    const Matrix& transpose(const Matrix& M);
};

#endif

