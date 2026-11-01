#include "fvDiffusion2D.hpp"

using namespace INMOST;

FV_Diffusion2D::FV_Diffusion2D(ProcessData_Diffusion2D *pd) :
	mesh(pd->get_mesh()), processdata(pd) {

}

void FV_Diffusion2D::build(){
	//Создаем или получаем тег проводимости
	std::string tagname = "TRANSMISSIBILITY";
	tagname += "_" + processdata->getName();
	if (!mesh->HaveTag(tagname))
		tagT = mesh->CreateTag(tagname, INMOST::DATA_REAL, INMOST::FACE, INMOST::NONE, 1);
	else
		tagT = mesh->GetTag(tagname);
	for (auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
		Face face = iface->getAsFace();
		double xFace[3];
		face.Barycenter(xFace);

		if (face.Boundary()){
			// для Неймановских граней проводимость не нужна
			if (processdata->getBCType(face) == BCType::Neumann) continue;
			Cell bCell = face.BackCell();
			double xBCell[3];
			bCell.Barycenter(xBCell);
			rMatrix nFace(3,1),kFC(3,1);
			//нормаль к грани
			face.UnitNormal(nFace.data());
			//нормированный вектор от грани к ячейке
			for (unsigned k = 0; k<3; k++)
				kFC(k,0) = xFace[k] - xBCell[k];
			double modkFC = kFC.FrobeniusNorm();
			kFC /= modkFC*modkFC;
			//сюда будем записывать результат умножения вектора kFC на тензор
			rMatrix D_kFC(3,1);
			processdata->multA(bCell,kFC.data(),D_kFC.data());
			//умножаем на нормаль
			double trans = D_kFC.DotProduct(nFace);
			face.Real(tagT) = trans;
		}
		else{//обработка внутренних граней
			Cell bCell = face.BackCell();
			Cell fCell = face.FrontCell();
			double xBCell[3],xFCell[3];
			bCell.Barycenter(xBCell);
			fCell.Barycenter(xFCell);
			rMatrix nFace(3,1),kFB(3,1),kFF(3,1);
			face.UnitNormal(nFace.data());
			for (unsigned k = 0; k<3; k++){
				kFB(k,0) = xFace[k] - xBCell[k];
				kFF(k,0) = xFace[k] - xFCell[k];
			}
			double modkFB = kFB.FrobeniusNorm();
			double modkFF = kFF.FrobeniusNorm();
			kFB /= modkFB*modkFB;
			kFF /= modkFF*modkFF;
			rMatrix D_kFB(3,1),D_kFF(3,1);
			processdata->multA(bCell,kFB.data(),D_kFB.data());
			processdata->multA(fCell,kFF.data(),D_kFF.data());
			double transB = D_kFB.DotProduct(nFace);
			double transF = D_kFF.DotProduct(nFace);
			//минус, потому что для FrontFace нормаль - внутренняя
			double trans = - transB * transF / (transB - transF);
			face.Real(tagT) = trans;
		}
	}
	std::cout << "TPFA is built"<<std::endl;
}

variable FV_Diffusion2D::getDgradH(const INMOST::Face &face, INMOST::dynamic_variable &H){
	variable flux;
	double length = face.Area();
	//минусы в выражениях, потому что потока положительный - градиент отрицательный
	if (face.Boundary()){
		//в Неймановских поток записан непосредственно
		if (processdata->getBCType(face)== BCType::Neumann){
			return -processdata->getBCValue(face)*length*processdata->getThickness(face.BackCell());
		}
		//условие Дирихле
		Cell BCell = face.BackCell();
		double value = processdata->getBCValue(face);
		flux = -face.Real(tagT) * (H(BCell) - value)*length;//умножение на мощность пласта сидит в проводимости
	}
	else{
		Cell BCell = face.BackCell(), FCell = face.FrontCell();
		flux = -face.Real(tagT) * (H(BCell) - H(FCell))*length;//умножение на мощность пласта сидит в проводимости
	}
	return flux;
}
