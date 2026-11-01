#pragma once

#include "inmost.h"
#include "Material.hpp"
#include <string>
#include <unordered_map>
#include <variant>
#include <functional>
#include "WellManager.hpp"

#include "json.hpp"
using json = nlohmann::json;

//умный обработчик считывания из json
template<typename T>
T require_path(const nlohmann::json& j, const std::string& path) {
	std::istringstream ss(path);
	std::string token;
	const nlohmann::json* current = &j;

	while (std::getline(ss, token, '.')) {
		if (!current->contains(token)) {
			std::cerr << "[MODEL READING ERROR] Missing key: \"" << token << "\" in path \"" << path << "\"" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		current = &(*current)[token];
	}

	try {
		return current->get<T>();
	} catch (const std::exception& e) {
		std::cerr << "[MODEL READING ERROR] Can't read key \"" << path << "\": " << e.what() << std::endl;
		std::exit(EXIT_FAILURE);
	}
}


enum class BCType { Dirichlet, Neumann };

using BCValue = std::variant<double, std::function<double(double, double)>>;

struct BoundaryConditionGW {
	BCType type;
	BCValue value;

	double Evaluate(double x, double y) const {
		if (std::holds_alternative<double>(value)) {
			return std::get<double>(value);
		} else {
			return std::get<std::function<double(double, double)>>(value)(x, y);
		}
	}
};

struct BoundaryConditionConc {
	BCType type;
	std::unordered_map<std::string, double> values;

	double Evaluate(const std::string& name) const {
		auto it = values.find(name);
		if (it == values.end())
		{
			std::cout<<("No BC value for component: " + name)<<std::endl;
			throw std::runtime_error("No BC value for component: " + name);
		}
		return it->second;
	}
};


class UserData{
public:
	UserData(INMOST::Mesh* m);
	void RegisterMaterial(int mat_id, const Material& mat);
	const Material& GetMaterial(int mat_id) const;

	void RegisterBCPressure(int bc_id, const BoundaryConditionGW& bc);
	void RegisterBCTransport(int bc_id, const BoundaryConditionConc& bc);
	const BoundaryConditionGW& GetBCPressure(int bc_id) const;
	const BoundaryConditionConc& GetBCTransport(int bc_id) const;
	const std::string& GetFlowSolverType() const {return flowSolverType;};
	void setFlowSolverType(const std::string type) {flowSolverType = type;};
	void multK(const INMOST::Cell&, double *v, double*res);
	void Init(json& config);
	const BCType getBCPressureType(int i) {return bc_pressure_map[i].type;}
	const BCType getBCPressureType(INMOST::Face face) {return getBCPressureType(face.Integer(tag_bc_pressure));}
	const BCType getBCTransportType(int i) {return bc_transport_map[i].type;};
	const BCType getBCTransportType(INMOST::Face face){return getBCTransportType(face.Integer(tag_bc_transport));}
	const BCValue getBCPressureValue(int i) {return bc_pressure_map[i].value;}
	int getBCPressureNum(INMOST::Face face){return face.Integer(tag_bc_pressure);}
	double getBCPressureValue(INMOST::Face face);
	const std::unordered_map<int, BoundaryConditionGW>& getBCPressureMap(){return bc_pressure_map;}
	const double getBCTransportValue(int i, const std::string& name) {return bc_transport_map[i].Evaluate(name);};
	double getBCTransportValue(INMOST::Face face, const std::string& name);
	double getBCTransportValue(INMOST::Face face, const int i) {return getBCTransportValue(face,mobile_components[i]);}
	const std::unordered_map<int, BoundaryConditionConc>& getBCTransportMap(){return bc_transport_map;};
	const double getPermeability(INMOST::Cell cell);
	const double getThickness(INMOST::Cell cell) {return GetMaterial(cell.Integer(tag_mat_id)).thickness;};
	const double getStorage(INMOST::Cell cell) {return GetMaterial(cell.Integer(tag_mat_id)).storage;};
	const double getPorosity(INMOST::Cell cell);
	
	
	INMOST::Tag tag_mat_id;
	INMOST::Tag tag_bc_pressure, tag_bc_transport;
	INMOST::Tag FluxTag;
	INMOST::Mesh* get_mesh() {return mesh;};
	INMOST::Tag get_tag(std::string tagname);
	WellManager wm;
	void init_wells();
	void set_dx(double dx_){dx = dx_;}
	void set_dy(double dy_){dy = dy_;}
	const double get_dx() const{return dx;}
	const double get_dy() const{return dy;}
	void set_Lx(double Lx_){Lx = Lx_;}
	void set_Ly(double Ly_){Ly = Ly_;}
	const double get_Lx() const{return Lx;}
	const double get_Ly() const{return Ly;}
	

	bool steady;
	bool transport;
	bool warnCFL;
	double p0;
	double dt;
	unsigned Nt;
	unsigned save_intensity;
	double Tend;


	void RegisterMobComponent(const std::string& name);

	// Доступ
	int    GetMobIndex(const std::string& name) const;


	// Получить «каталоги» имён
	const std::vector<std::string>& GetMobNames() const { return mobile_components; }

	//метод для задания НУ переноса
	void SetInitialConcentration(const std::string& name, double value, bool is_mobile);
	//НУ для определнной компоненты
	double GetTransportInitialConditions(const std::string& name, bool is_mobile);

	
private:
	INMOST::Mesh* mesh;
	//Пул материалов
	std::unordered_map<int, Material> material_map;
	//Пул граничных условий для фильтрации
	std::unordered_map<int, BoundaryConditionGW> bc_pressure_map;
	//Пул граничных условий для переноса
	std::unordered_map<int, BoundaryConditionConc> bc_transport_map;
	//Тип линейного решателя для фильтрации
	std::string flowSolverType;
	//имена мобильных компонент
	std::vector<std::string> mobile_components;


	//НУ мобильных компонент
	std::map<std::string, double> initial_mobile_concentrations;
	double dx,dy,Lx,Ly;
};
