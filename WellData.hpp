#pragma once

#include "inmost.h"
#include <functional>
#include <string>



enum class WellType { BHP,//с фиксированным забойным давлением
					  RATE //с заданным дебитом
					}; // возможно расширение

class WellData {
public:
	WellData(int id,
			 double x_,
			 double y_,
			 WellType type_,
			 double pw_,
			 double rw_,
			 double skin_,
			 double fix_rate);

	int getID() const;
	double getX() const;
	double getY() const;
	void setPermeability(const double k){permeability = k;}
	void setThickness(const double b){thickness = b;}
	void set_dx(double dx_){dx = dx_;}
	void set_dy(double dy_){dy = dy_;}
	const double get_dx() const{return dx;}
	const double get_dy() const{return dy;}

	WellType getType() const;
	double getBHP() const;
	bool isActive() const;
	void markActive(int global_id);
	int getLocalCellID() const;
	double computeWellIndex() const;
	const double getRate() const{return fix_rate;}

	// Установить концентрацию для компоненты
	void SetComponentConcentration(const std::string& name, double value);
	// Получить концентрацию
	double GetComponentConcentration(const std::string& name) const;
	// Получить всю карту концентраций
	const std::unordered_map<std::string, double>& GetComponentConcentrations() const;

private:
	int id;
	double x, y;
	WellType type;
	double pw_const;
	double rw;
	double skin;
	double permeability;
	double thickness;
	double dx,dy;
	double fix_rate;
	std::unordered_map<std::string, double> component_concentrations;

	bool active = false;
	int local_cell_id = -1;
};
