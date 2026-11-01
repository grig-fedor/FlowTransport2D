#pragma once
#include "WellData.hpp"
#include "inmost.h"
#include <vector>

class WellManager {
public:
	void addWell(const WellData& well);
	void bindToMesh(INMOST::Mesh* mesh, double dx, double dy); // поиск ячейки для каждой скважины

	const std::vector<WellData>& getWells() const;
	std::vector<WellData>& getWellsMutable(); // если нужно менять

private:
	std::vector<WellData> wells;
};
