#include "processdata.hpp"

ProcessData_SaturatedFlow2D::ProcessData_SaturatedFlow2D(UserData *ud) :
	ProcessData_Diffusion2D(ud)
{
	setName("Flow");
}

void ProcessData_SaturatedFlow2D::multA(const INMOST::Cell& cell, double *v, double*res){
	UD->multK(cell,v,res);
}

const BCType ProcessData_SaturatedFlow2D::getBCType(const INMOST::Face &face){
	return UD->getBCPressureType(face.Integer(UD->tag_bc_pressure));
}

double ProcessData_SaturatedFlow2D::getBCValue(const INMOST::Face &face){
	return UD->getBCPressureValue(face);
}
