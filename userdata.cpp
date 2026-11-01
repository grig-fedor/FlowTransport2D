#include "userdata.hpp"
#include <string>

using namespace  INMOST;
using namespace std;

UserData::UserData(Mesh* m) : mesh(m){};

// Задание параметров задачи, разметка граничных условий, задание параметров численной модели
void UserData::Init(json& config){
	std::cout << "[Initialization] Saturated flow initialization..." << std::endl;

	Material mat_default;
	const auto& material_dict = config["materials"]["Default"];
	mat_default.thickness    = require_path<double>(material_dict, "thickness");
	mat_default.permeability = require_path<double>(material_dict, "permeability");
	mat_default.storage      = require_path<double>(material_dict, "storage");
	mat_default.porosity     = require_path<double>(material_dict, "porosity");
	mat_default.diffusivity  = require_path<double>(material_dict, "diffusivity");
	mat_default.alphaL       = require_path<double>(material_dict, "alphaL");
	mat_default.alphaT       = require_path<double>(material_dict, "alphaT");
	RegisterMaterial(0, mat_default);
	tag_mat_id = mesh->CreateTag("MATERIAL_ID", DATA_INTEGER, CELL, NONE, 1);
	//задаем материал по умолчанию во всей области
	for (Mesh::iteratorCell icell = mesh->BeginCell(); icell != mesh->EndCell(); ++icell) {
		icell->Integer(tag_mat_id) = 0;
	}

	//Обработка параметров фильтрационной задачи
	const auto& filtr_config = config["flow"];
	// Регистрация ГУ
	if (filtr_config.contains("bc_pressure") && filtr_config["bc_pressure"].contains("list")) {
		const auto& bc_list = filtr_config["bc_pressure"]["list"];

		for (auto& [key_str, bc_data] : bc_list.items()) {
			int id = std::stoi(key_str);
			std::string type_str = bc_data["type"];
			BCType type;

			if (type_str == "Neumann") type = BCType::Neumann;
			else if (type_str == "Dirichlet") type = BCType::Dirichlet;
			else {
				std::cerr << "Unknown BC type: " << type_str << std::endl;
				std::exit(EXIT_FAILURE);
			}

			double value = bc_data["value"];
			RegisterBCPressure(id, {type, value});
		}
	}
	// Пример ГУ, заданного значением
	RegisterBCPressure(100, {BCType::Neumann, 0.0});
	// Пример ГУ, заданного функцией
	RegisterBCPressure(99, {
		BCType::Dirichlet,
		[](double x, double y) { return x; }
	});

	tag_bc_pressure = mesh->CreateTag("BC_PRESSURE", DATA_INTEGER, FACE, NONE, 1);
	std::unordered_map<std::string, int> pressure_bc_sides;
	if (filtr_config["bc_pressure"].contains("sides")) {
		for (auto& [side, bc_id] : filtr_config["bc_pressure"]["sides"].items()) {
			pressure_bc_sides[side] = bc_id;
		}
	}

	for (Mesh::iteratorFace iface = mesh->BeginFace(); iface != mesh->EndFace(); ++iface) {
		Face face = iface->getAsFace();
		if (face.Boundary()) {
			double xFace[3];
			face.Barycenter(xFace);
			double x0 = require_path<double>(config, "mesh.x0");
			double y0 = require_path<double>(config, "mesh.y0");
			double x = xFace[0] - x0, y = xFace[1] - y0;
			int bc_id = -1;
			double eps = require_path<double>(config, "num_parameters.epsilon");
			if (std::abs(y - Ly) < dy*eps) bc_id = pressure_bc_sides["north"];
			else if (std::abs(x - Lx) < dx*eps) bc_id = pressure_bc_sides["east"];
			else if (std::abs(y - 0.0) < dy*eps) bc_id = pressure_bc_sides["south"];
			else if (std::abs(x - 0.0) < dx*eps) bc_id = pressure_bc_sides["west"];
			if (bc_id >=0 ) face.Integer(tag_bc_pressure) = bc_id;
		}
	}

	setFlowSolverType(require_path<std::string>(config, "num_parameters.linear_solver_type"));

	//стационарность задачи
	steady = require_path<bool>(config, "flow.steady");
	//начальное пластовое давление для нестационарной задачи
	if (!steady) p0 = require_path<double>(config, "flow.IC");
	//включаем перенос
	transport = config.contains("transport");
	if (!steady || transport){
		//шаг по врмени для нестационарной и транспортной задачи
		dt = require_path<double>(config, "num_parameters.time_step");
		//число шагов по времени
		Nt = require_path<int>(config, "num_parameters.time_step_number");
		//как часто сохранять решение
		save_intensity = require_path<int>(config, "num_parameters.save_intensity");
		Tend = Nt*dt;
		if (config["num_parameters"].contains("Tend")) Tend = require_path<double>(config, "num_parameters.Tend");
		warnCFL = false;
		if (config["num_parameters"].contains("warnCFL")) warnCFL = require_path<bool>(config, "num_parameters.warnCFL");
	}

	

	//заводим тег под поток
	FluxTag = mesh->CreateTag("FLUX", DATA_REAL, FACE, NONE, 1);

	if (transport){
		const auto& tr = config["transport"];
		// Подвижные компоненты
		std::vector<std::string> mob_comps = require_path<std::vector<std::string>>(tr, "mob_components");
		for (auto name : mob_comps )
		{
			RegisterMobComponent(name);
			double IC = 0;
			if (tr.contains("initial_conditions") &&
				tr["initial_conditions"].contains("mobile") &&
				tr["initial_conditions"]["mobile"].contains(name))
					IC = require_path<double>(tr, "initial_conditions.mobile."+name);
			SetInitialConcentration(name, IC, true);
		}

		//Обрабатываем ГУ переноса
		if (tr.contains("bc_transport") && tr["bc_transport"].contains("list")) {
			const auto& bc_list = tr["bc_transport"]["list"];
			for (auto& [key_str, bc_data] : bc_list.items()) {
				int id = std::stoi(key_str);
				std::string type_str = bc_data["type"];
				BCType type;

				if (type_str == "Neumann") type = BCType::Neumann;
				else if (type_str == "Dirichlet") type = BCType::Dirichlet;
				else {
					std::cerr << "Unknown transport BC type: " << type_str << std::endl;
					std::exit(EXIT_FAILURE);
				}

				std::unordered_map<std::string, double> values;
                for (const auto& comp : mobile_components) {
                    values[comp] = require_path<double>(bc_data, "values." + comp);
                }
				RegisterBCTransport(id, BoundaryConditionConc{type, values});
			}
		}

		// Расстановка ГУ переноса для задачи поршневого вытеснения в колонке слева направо
		tag_bc_transport = mesh->CreateTag("BC_TRANSPORT", DATA_INTEGER, FACE, NONE, 1);
		std::unordered_map<std::string, int> transport_bc_sides;
		if (tr["bc_transport"].contains("sides")) {
			for (auto& [side, bc_id] : tr["bc_transport"]["sides"].items()) {
				transport_bc_sides[side] = bc_id;
			}
		}

		for (Mesh::iteratorFace iface = mesh->BeginFace(); iface != mesh->EndFace(); ++iface) {
			Face face = iface->getAsFace();
			if (face.Boundary()) {
				double xFace[3];
				face.Barycenter(xFace);
				double x0 = require_path<double>(config, "mesh.x0");
				double y0 = require_path<double>(config, "mesh.y0");
				double x = xFace[0] - x0, y = xFace[1] - y0;
				int bc_id = -1;
				double eps = require_path<double>(config, "num_parameters.epsilon");
				if (std::abs(y - Ly) < dy*eps) bc_id = transport_bc_sides["north"];
				else if (std::abs(x - Lx) < dx*eps) bc_id = transport_bc_sides["east"];
				else if (std::abs(y - 0.0) < dy*eps) bc_id = transport_bc_sides["south"];
				else if (std::abs(x - 0.0) < dx*eps) bc_id = transport_bc_sides["west"];
				if (bc_id >=0 ) face.Integer(tag_bc_transport) = bc_id;
			}
		}

	}

	//добавляем скважины
	if (config.contains("wells")){
		const auto& wells = config["wells"];
		for (const auto& w : wells)
		{
			WellType type = WellType::RATE;
			std::string t = w["type"];
			if (t == "BHP") type = WellType::BHP;
			else if (t == "RATE") type = WellType::RATE;
			else {
				std::cerr << "Unknown well type: " << t << std::endl;
				std::exit(EXIT_FAILURE);
			}
			WellData well(
						require_path<int>(w, "id"),
						require_path<double>(w, "x"),
						require_path<double>(w, "y"),
						type,
						require_path<double>(w, "pw"),
						require_path<double>(w, "radius"),
						require_path<double>(w, "radius"),
						require_path<double>(w, "rate")
					);
			if (transport){
                if (w.contains("concentrations")) {
                    const auto& concs = w["concentrations"];
                    for (const auto& comp : mobile_components)
                        well.SetComponentConcentration(comp, require_path<double>(concs, comp));
                } else {
                    std::cerr << "Can't read concentrations for well: " << w["id"] << std::endl;
                    std::exit(EXIT_FAILURE);
                }
			}
			wm.addWell(well);
		}
		//init_wells();
	}
}



// Задание параметров задачи, разметка граничных условий, задание параметров численной модели


void UserData::RegisterMaterial(int mat_id, const Material& mat) {
	material_map[mat_id] = mat;
}

const Material& UserData::GetMaterial(int mat_id) const {
	auto it = material_map.find(mat_id);
	if (it == material_map.end()) throw std::runtime_error("Unknown material id");
	return it->second;
}

void UserData::RegisterBCPressure(int bc_id, const BoundaryConditionGW& bc) {
	bc_pressure_map[bc_id] = bc;
}

void UserData::RegisterBCTransport(int bc_id, const BoundaryConditionConc& bc) {
	bc_transport_map[bc_id] = bc;
}

const BoundaryConditionGW& UserData::GetBCPressure(int bc_id) const {
	auto it = bc_pressure_map.find(bc_id);
	if (it == bc_pressure_map.end()) throw std::runtime_error("Unknown pressure BC id");
	return it->second;
}

const BoundaryConditionConc& UserData::GetBCTransport(int bc_id) const {
	auto it = bc_transport_map.find(bc_id);
	if (it == bc_transport_map.end()) throw std::runtime_error("Unknown transport BC id");
	return it->second;
}

void UserData::multK(const INMOST::Cell & cell, double *v, double *res){
	double Kf = getPermeability(cell);
	double b = GetMaterial(cell.Integer(tag_mat_id)).thickness;
	for (int i = 0; i<2; i++) res[i] = v[i]*Kf*b;
}


double UserData::getBCPressureValue(INMOST::Face face){
	// Получаем ID граничного условия из тега
	int bc_id = face.Integer(get_tag("BC_PRESSURE"));

	// Получаем само граничное условие
	const BoundaryConditionGW& bc = bc_pressure_map.at(bc_id);

	// Получаем координаты центра грани
	double xyz[3];
	face.Centroid(xyz);
	double x = xyz[0];
	double y = xyz[1];

	// Вычисляем значение (константа или функция)
	return bc.Evaluate(x, y);
}

 double UserData::getBCTransportValue(INMOST::Face face, const std::string& name){
	int bc_id = face.Integer(get_tag("BC_TRANSPORT"));
	return getBCTransportValue(bc_id,name);
}

void UserData::init_wells(){
	wm.bindToMesh(get_mesh(),dx,dy);
	std::vector<WellData>& wells = wm.getWellsMutable();
	for (auto& w : wells) {
		if (!w.isActive()) continue;
		INMOST::Cell c = mesh->CellByLocalID(w.getLocalCellID());
		w.setPermeability(getPermeability(c));
		w.setThickness(getThickness(c));
	}
}

void UserData::RegisterMobComponent(const std::string& name) {
	// Проверка: уже зарегистрирована?
	if (std::find(mobile_components.begin(), mobile_components.end(), name) != mobile_components.end())
		return;

	mobile_components.push_back(name);
}





int UserData::GetMobIndex(const std::string& name) const {
	auto it = std::find(mobile_components.begin(), mobile_components.end(), name);
	if (it == mobile_components.end()){
		std::cerr << "Unknown mobile component "<< name << std::endl;
		std::exit(EXIT_FAILURE);
	}
	return static_cast<int>(std::distance(mobile_components.begin(), it));
}







void UserData::SetInitialConcentration(const std::string& name, double value, bool is_mobile) {
    if (is_mobile) initial_mobile_concentrations[name] = value;
}

double UserData::GetTransportInitialConditions(const std::string& name, bool is_mobile) {
    return is_mobile ? initial_mobile_concentrations.at(name) : 0.0;
}

const double UserData::getPorosity(INMOST::Cell cell) {
    return GetMaterial(cell.Integer(tag_mat_id)).porosity;
}

const double UserData::getPermeability(INMOST::Cell cell) {
    return GetMaterial(cell.Integer(tag_mat_id)).permeability;
}

INMOST::Tag UserData::get_tag(string tagname){
	if (mesh->HaveTag(tagname)) return mesh->GetTag(tagname);
	else{
		std::cerr << "There is no "<< tagname << " tag on the mesh!"<< std::endl;
		std::exit(EXIT_FAILURE);
	}
}
