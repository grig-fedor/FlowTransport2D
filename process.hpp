#pragma once
#include "inmost.h"
#include "fvDiffusion2D.hpp"
class Process
{
protected:
	INMOST::Mesh *m;
	bool steady;
public:
	Process(INMOST::Mesh *mesh) {m = mesh;}
	virtual ~Process(){}
	virtual void fillResidual(INMOST::Residual &R) = 0;
	void setSteady(bool b) { steady = b;}
};


class process_saturatedFlow2D: public Process
{
private:
	FV_Diffusion2D * difDiscr;
	INMOST::dynamic_variable varH;
	ProcessData_SaturatedFlow2D* processdata;
	UserData* UD;
	INMOST::Tag oldH;
	bool fluxes_built;
public:
	~process_saturatedFlow2D(){
		if (processdata != NULL) delete processdata;
		if (difDiscr != NULL) delete difDiscr;
	};
	process_saturatedFlow2D(UserData *ud, std::vector<INMOST::dynamic_variable> &dvars);
	void fillResidual(INMOST::Residual & R);
	INMOST::variable getFlux(const INMOST::Face &f);
	void storeFluxes(INMOST::Tag);
};


