#include "process.hpp"

using namespace INMOST;

process_saturatedFlow2D::process_saturatedFlow2D(UserData *ud, std::vector<INMOST::dynamic_variable> &dvars):
	Process(ud->get_mesh()), UD(ud)
{
	//готовим данные для дискретизации
	processdata = new ProcessData_SaturatedFlow2D(UD);
	//конструируем дискретизацию, передавая необходимые данные
	difDiscr = new FV_Diffusion2D(processdata);
	difDiscr->build();
	fluxes_built = true;
	varH = dvars[0];
	oldH = ud->get_tag("PRESSURE");
}

void process_saturatedFlow2D::fillResidual(INMOST::Residual & R){
	if (!fluxes_built) difDiscr->build();
	for (auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
		Face face = iface->getAsFace();
		Cell BCell = face.BackCell(), FCell = face.FrontCell();
		variable q = -difDiscr->getDgradH(face,varH);
		R[varH.Index(BCell)] -= q;
		if (FCell.isValid()) R[varH.Index(FCell)] += q;
	}
	//цикл по скважинам
	for (const auto& w : UD->wm.getWells()) {
		if (!w.isActive()) continue;
		INMOST::Cell c = m->CellByLocalID(w.getLocalCellID());
		double bar[3];
		c.Barycenter(bar);
		std::cout << "xWell = "<< bar[0] << " yWell = "<< bar[1] <<  " type = " << int(w.getType())<<std::endl;
		if (w.getType() == WellType::BHP){
			double WI = w.computeWellIndex();
			double pw = w.getBHP();
			std::cout << "pw = "<<pw<<std::endl;
			std::cout << "WI = "<<WI<<std::endl;
			R[varH.Index(c)] -= WI * (varH(c) - pw);
		}
		else if (w.getType() == WellType::RATE){
			R[varH.Index(c)] += w.getRate();
		}
	}
	//производная по времени
	if (!steady){
		for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
			Cell cell = icell->getAsCell();
			//std::cout << "V = "<<cell.Volume()<<std::endl;
			R[varH.Index(cell)] -= UD->getStorage(cell) * (varH(cell) - cell.Real(oldH)) / UD->dt * cell.Volume() * UD->getThickness(cell);
			//R[varH.Index(cell)] -= UD->getStorage(cell) * (varH(cell) - cell.Real(oldH)) / UD->dt * UD->get_dx() * UD->get_dy() * UD->getThickness(cell);
		}
	}
}

variable process_saturatedFlow2D::getFlux(const Face &f){
	return -difDiscr->getDgradH(f,varH);
}

void process_saturatedFlow2D::storeFluxes(INMOST::Tag flux_tag){
	for (auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
		iface->Real(flux_tag) = getFlux(iface->getAsFace()).GetValue();
	}
}
