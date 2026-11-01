#include "WellManager.hpp"
#include <stdexcept>
#include <cmath>

void WellManager::addWell(const WellData& well) {
	wells.push_back(well);
}

void WellManager::bindToMesh(INMOST::Mesh* mesh, double dx, double dy) {
	for (auto& well : wells) {
		double xyz[3] = {well.getX(), well.getY(),0.0};
		well.set_dx(dx);
		well.set_dy(dx);
		for (auto icell = mesh->BeginCell(); icell != mesh->EndCell(); ++icell) {
			INMOST::Cell cell = icell->getAsCell();
			double bar[3];
			cell.Barycenter(bar);

			if (std::abs(bar[0]-well.getX())<dx &&
				std::abs(bar[1]-well.getY())<dy){
				well.markActive(cell.LocalID());
				std::cout << "xCell = "<< bar[0] << " yCell = "<< bar[1] << std::endl;
				break;
			}
			// Почему-то вариант ниже с Inside работает неправильно!
//			if (cell.Inside(xyz)) {
//				std::cout << "xWell = "<< well.getX() << " yWell = "<< well.getY() << std::endl;
//				double bar[3];
//				cell.Barycenter(bar);
//				std::cout << "xCell = "<< bar[0] << " yCell = "<< bar[1] << std::endl;
//				well.markActive(cell.LocalID());
//				break;
//			}
		}
	}
}

const std::vector<WellData>& WellManager::getWells() const {
	return wells;
}

std::vector<WellData>& WellManager::getWellsMutable() {
	return wells;
}
