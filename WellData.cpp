#include "WellData.hpp"
#include <cmath>
#define M_PI           3.14159265358979323846

WellData::WellData(int id_,
				   double x_, double y_,
				   WellType type_,
				   double pw_,//забойное давление
				   double rw_,//радиус скважины
				   double skin_,
				   double rate_)
	: id(id_), x(x_), y(y_), type(type_), pw_const(pw_), rw(rw_), skin(skin_),
	  permeability(0), thickness(0), fix_rate(rate_){}

int WellData::getID() const { return id; }
double WellData::getX() const { return x; }
double WellData::getY() const { return y; }

WellType WellData::getType() const { return type; }

bool WellData::isActive() const { return active; }

void WellData::markActive(int local_id) {
	active = true;
	local_cell_id = local_id;
}

int WellData::getLocalCellID() const { return local_cell_id; }

double WellData::getBHP() const {
	 return pw_const;
}

//расчет well_index по формуле Писмана
double WellData::computeWellIndex() const {
	double re = 0.14 * sqrt(dx*dx + dy*dy); // эквивалентный радиус дренажа
	return (2.0 * M_PI * permeability * thickness) / (log(re / rw) + skin);
}

void WellData::SetComponentConcentration(const std::string &name, double value)
{
	component_concentrations[name] = value;
}

double WellData::GetComponentConcentration(const std::string &name) const
{
	auto it = component_concentrations.find(name);
	if (it == component_concentrations.end())
		throw std::runtime_error("No concentration set in well for component: " + name);
	return it->second;
}
