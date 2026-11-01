#include "problem.hpp"

using namespace INMOST;
using namespace std;

//#define SAVE_RUN_DATA

Problem::Problem(UserData* ud_) : mesh(ud_->get_mesh()),ud(ud_) {
	// Создание тега давления, если не создан
	if (!mesh->HaveTag("PRESSURE")) {
		pressure_tag = mesh->CreateTag("PRESSURE", DATA_REAL, CELL, NONE, 1);
	} else {
		pressure_tag = mesh->GetTag("PRESSURE");
	}

	//для нестационарной задачи инициализируем пластовое давление
	if (!ud->steady){
		for (auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
			icell->Real(pressure_tag) = ud->p0;
		}
	}

	flux_tag = ud->get_tag("FLUX");
	if (ud->transport){

		if (!mesh->HaveTag("MOB_CONC"))
			mob_conc_tag = mesh->CreateTag("MOB_CONC",DATA_REAL,CELL,NONE,ud->GetMobNames().size());
		else
			mob_conc_tag = mesh->GetTag("MOB_CONC");
		//расстановка НУ транспорта на сетке
		SetTranIC();
		mob_conc_old_tag = mob_conc_tag;
		
	}
	if (ud->transport || !ud->steady) Save("step_0.vtk");
	ud->init_wells();
}



void Problem::Save(const std::string& filename) {
    std::vector<INMOST::Tag> created_tags;
    if (ud->transport && mesh->HaveTag("MOB_CONC")) {
        INMOST::Tag concentrations = mesh->GetTag("MOB_CONC");
        const auto& names = ud->GetMobNames();
        for (int i = 0; i < static_cast<int>(names.size()); ++i) {
            INMOST::Tag component;
            if (mesh->HaveTag(names[i])) component = mesh->GetTag(names[i]);
            else {
                component = mesh->CreateTag(names[i], INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
                created_tags.push_back(component);
            }
            for (auto cell = mesh->BeginCell(); cell != mesh->EndCell(); ++cell)
                cell->Real(component) = cell->RealArray(concentrations)[i];
        }
    }
    mesh->Save(filename);
    for (auto tag : created_tags) mesh->DeleteTag(tag);
}

void Problem::RunSteadyFlow(bool store_fluxes) {
	std::cout << "[Run] Confined flow running..." << std::endl;


	//Указываю текущий автоматизатор
	Automatizator::MakeCurrent(&aut);
	//Регестрирую переменную давления, получаю ее индекс
	auto indH = aut.RegisterTag(pressure_tag, CELL);
	//Завожу переменную под давление
	dynamic_variable varH(aut, indH);
	//Назначем индексы (положения в глобальном векторе) зарегестрированным переменным
	aut.EnumerateEntries();
	residual = Residual("R", aut.GetFirstIndex(), aut.GetLastIndex());

	//Загружаю переменные для решения задачи фильтрации
	std::vector<dynamic_variable> varsFlow;
	varsFlow.push_back(varH);

	process_saturatedFlow2D Flow(ud,varsFlow);

	Flow.setSteady(ud->steady);

	INMOST::Solver S(ud->GetFlowSolverType(), "flow");
	S.SetParameter("relative_tolerance", "1e-12");
	S.SetParameter("absolute_tolerance", "1e-15");

	Sparse::Vector sol("sol", aut.GetFirstIndex(), aut.GetLastIndex());

	residual.Clear();
	Flow.fillResidual(residual);
	double norm2 = 0.0;
	norm2 = residual.Norm();
	std::cout << "Norm after fillResidual "<< norm2 << std::endl;
	S.SetMatrix(residual.GetJacobian());
	bool solved = S.Solve(residual.GetResidual(),sol);
	if(!solved){
		std::cout << "Linear solver failed: " << S.GetReason() << std::endl;
		std::cout << "Residual: " << S.Residual() << std::endl;
		exit(1);
	}
	for (auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++)
		icell->Real(pressure_tag) = -sol[varH.Index(icell->getAsCell())];
	if (store_fluxes){
		Flow.storeFluxes(flux_tag);
	}
	std::stringstream filename;
	filename << "Pressure" << ".vtk";
	Save(filename.str());
	cout << "Max dt = " << compute_max_dT()<<endl;
	std::cout << "Flow is done" << std::endl;
}


void Problem::RunUnsteadyFlow() {
	std::cout << "[Run] Confined flow running..." << std::endl;


	//Указываю текущий автоматизатор
	Automatizator::MakeCurrent(&aut);
	//Регестрирую переменную давления, получаю ее индекс
	auto indH = aut.RegisterTag(pressure_tag, CELL);
	//Завожу переменную под давление
	dynamic_variable varH(aut, indH);
	//Назначем индексы (положения в глобальном векторе) зарегестрированным переменным
	aut.EnumerateEntries();
	residual = Residual("R", aut.GetFirstIndex(), aut.GetLastIndex());

	//Загружаю переменные для решения задачи фильтрации
	std::vector<dynamic_variable> varsFlow;
	varsFlow.push_back(varH);

	process_saturatedFlow2D Flow(ud,varsFlow);

	Flow.setSteady(ud->steady);

	INMOST::Solver S(ud->GetFlowSolverType(), "flow");
	S.SetParameter("relative_tolerance", "1e-12");
	S.SetParameter("absolute_tolerance", "1e-15");

	Sparse::Vector sol("sol", aut.GetFirstIndex(), aut.GetLastIndex());

	for (unsigned time_step = 0; time_step<ud->Nt; time_step++){
		AssembleSolveUpdateFlow(Flow, varH, S, sol);
		std::cout << "Flow on step "<< time_step<< " is done" << std::endl;
	}
	std::cout << "Flow is done" << std::endl;
}


void Problem::RunSteadyFlowTransport(){
	RunSteadyFlow(true);
	for (unsigned time_step = 0; time_step<ud->Nt; time_step++){
		make_explicit_advection_step();
		
		
		std::cout << "Transport on step "<< time_step<< " is done" << std::endl;
		if (time_step % ud->save_intensity == 0){
			std::stringstream filename;
			filename << "step_" << (time_step/ud->save_intensity + 1) << ".vtk";
			Save(filename.str());
			std::cout << "Saved: " << filename.str() << std::endl;
		}
	}
}

void Problem::make_explicit_advection_step(){
	int comp_number = mob_conc_tag.GetSize();
	int cell_number = mesh->NumberOfCells();
	double **ConcNew = new double*[cell_number];
	double **source = new double*[cell_number];
	for (int i = 0; i<cell_number; i++) ConcNew[i] = new double[comp_number];
	for (int i = 0; i<cell_number; i++) source[i] = new double[comp_number];
	bool *isSource = new bool[cell_number];
	//расчет перетока концентрации через скважины
	run_sources(source,mob_conc_tag,isSource);

	for (Mesh::iteratorCell icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
		ElementArray<Face> faces = icell->getFaces();
		Cell cell = icell->getAsCell();
		double *MassFlux = new double[comp_number];
		for (int i=0; i<comp_number; i++) MassFlux[i] = 0;
		for (ElementArray<Face>::iterator iface = faces.begin(); iface != faces.end(); iface++){
			Face face = iface->getAsFace();
			double Flux = (face->BackCell() == cell) ? face.Real(flux_tag) :
													   - face.Real(flux_tag);
			if (Flux >= 0)
				for (int i=0; i<comp_number; i++)
					MassFlux[i] -= cell.RealArray(mob_conc_old_tag)[i] * Flux;
			else{
				//определеяем противопотоковое значение
				if (face->Boundary()){
					if (ud->getBCTransportType(face) == BCType::Dirichlet)
						for (int i=0; i<comp_number; i++)
							MassFlux[i] -= ud->getBCTransportValue(face,i) * Flux;
				}
				else
				if (face->BackCell() == cell)
					for (int i=0; i<comp_number; i++)
						MassFlux[i]	-= face->FrontCell()->RealArray(mob_conc_old_tag)[i] * Flux;
				else{
					for (int i=0; i<comp_number; i++)
						MassFlux[i]	-= face->BackCell()->RealArray(mob_conc_old_tag)[i] * Flux;
				}
			}
		}
		int cellIndex = cell->DataLocalID();
		double V = cell->Volume()*ud->getThickness(cell);
		double coef = ud->dt/V/ud->getPorosity(cell);
		for (int i=0; i<comp_number; i++)
			ConcNew[cellIndex][i] = cell->RealArray(mob_conc_old_tag)[i] + MassFlux[i] * coef;
		delete[] MassFlux;
		if (isSource[cellIndex]){
			for (int i=0; i<comp_number; i++){
				double temp  = source[cellIndex][i]* coef;
				ConcNew[cellIndex][i] += source[cellIndex][i]* coef;
			}
		}
	}
	for (Mesh::iteratorCell icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
		for (int i = 0; i<comp_number; i++)
			icell->RealArray(mob_conc_tag)[i] = ConcNew[icell->DataLocalID()][i];
	}
	for (int i = 0; i < cell_number; i++) delete[] ConcNew[i];
	delete[] ConcNew;
	for (int i = 0; i < cell_number; i++) delete[] source[i];
	delete[] source;
	delete[] isSource;
}

void Problem::SetTranIC()
{
    // Initial concentrations of transported components.
    for (const auto& name : ud->GetMobNames()) {
        int id = ud->GetMobIndex(name);
        double value = ud->GetTransportInitialConditions(name, true);
        for (auto cell = mesh->BeginCell(); cell != mesh->EndCell(); ++cell)
            cell->RealArray(mob_conc_tag)[id] = value;
    }
}

double Problem::get_well_rate(const WellData& w) const {
	if (!w.isActive()) return 0.0;
	if (w.getType() == WellType::RATE) return w.getRate();
	else if (w.getType() == WellType::BHP){
		const INMOST::Cell c = mesh->CellByLocalID(w.getLocalCellID());
		double WI = w.computeWellIndex();
		double pw = w.getBHP();
		double p = c.Real(pressure_tag);
		return WI*(pw - p);
	}
	else{
		throw std::runtime_error("Unknown well type");
	}
}
void Problem::run_sources(double **source, INMOST::Tag conc_prev, bool *isSource){
	unsigned int cell_number = mesh->NumberOfCells();
	unsigned int comp_number = conc_prev.GetSize();
	for (unsigned int i = 0; i<cell_number; i++){
		for (unsigned int j = 0; j < comp_number; j++)
			source[i][j] = 0;
		isSource[i] = false;
	}
	//цикл по скважинам
	for (const auto& w : ud->wm.getWells()) {
		if (!w.isActive()) continue;
		INMOST::Cell c = mesh->CellByLocalID(w.getLocalCellID());
		unsigned int id_cell = w.getLocalCellID();
		isSource[id_cell] = true;
		double rate = get_well_rate(w);
		if (rate >= 0){
			std::vector<std::string> comp_names = ud->GetMobNames();
			for (unsigned int i = 0; i < comp_number; i++)
				source[id_cell][i] += rate*w.GetComponentConcentration(comp_names[i]);
		}
		else{
			for (unsigned int i = 0; i < comp_number; i++)
				source[id_cell][i] += rate*c.RealArray(conc_prev)[i];
		}
	}
}



void Problem::AssembleSolveUpdateFlow(process_saturatedFlow2D& Flow,
									 INMOST::dynamic_variable& varH,
									 INMOST::Solver& S,
									 INMOST::Sparse::Vector& sol, bool store_fluxes){
	residual.Clear();
	Flow.fillResidual(residual);

	double norm2 = residual.Norm();
	std::cout << "Norm after fillResidual " << norm2 << std::endl;
	S.SetMatrix(residual.GetJacobian());
	bool solved = S.Solve(residual.GetResidual(), sol);
	if (!solved) {
		std::cout << "Linear solver failed: " << S.GetReason() << std::endl;
		std::cout << "Residual: " << S.Residual() << std::endl;
		std::exit(1);
	}
	for (auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++)
		icell->Real(pressure_tag) -= sol[varH.Index(icell->getAsCell())];
	if (store_fluxes) Flow.storeFluxes(flux_tag);
}

void Problem::RunUnsteadyFlowTransport(){
	std::cout << "[Run] Unsteady flow + transport..." << std::endl;
	// Подготовка автоматизатора и неизвестных потока
	Automatizator::MakeCurrent(&aut);
	auto indH = aut.RegisterTag(pressure_tag, CELL);
	dynamic_variable varH(aut, indH);
	aut.EnumerateEntries();
	residual = Residual("R", aut.GetFirstIndex(), aut.GetLastIndex());
	std::vector<dynamic_variable> varsFlow;
	varsFlow.push_back(varH);
	process_saturatedFlow2D Flow(ud, varsFlow);
	Flow.setSteady(false); // именно НЕстационарный расчёт потока
	INMOST::Solver S(ud->GetFlowSolverType(), "flow");
	S.SetParameter("relative_tolerance", "1e-12");
	S.SetParameter("absolute_tolerance", "1e-15");
	INMOST::Sparse::Vector sol("sol", aut.GetFirstIndex(), aut.GetLastIndex());
	double t = 0;
	for (unsigned time_step = 0; t < ud->Tend; ++time_step){
		AssembleSolveUpdateFlow(Flow, varH, S, sol, true);
		std::cout << "Flow on step "<< time_step<< " is done" << std::endl;
		WriteBalances(t,ud->dt);
		make_explicit_advection_step();
		std::cout << "Transport on step "<< time_step<< " is done" << std::endl;
		if (time_step % ud->save_intensity == 0){
			std::stringstream filename;
			filename << "step_" << (time_step/ud->save_intensity + 1) << ".vtk";
			Save(filename.str());
			std::cout << "Saved: " << filename.str() << std::endl;
		}
		t+=ud->dt;
		double recom_dT = compute_max_dT();
		cout << "Recommended dt = "<< recom_dT << endl;
		if (ud->warnCFL) ud->dt = recom_dT;
	}
	std::cout << "Unsteady flow + transport is done" << std::endl;
}



double Problem::CellBulkVolume(const INMOST::Cell & c) const{
	// В 2D: Volume() = площадь; объём слоя = площадь * мощность
	return c.Volume() * ud->getThickness(c);
}

double Problem::CellPoreVolume(const INMOST::Cell & c) const{
	return CellBulkVolume(c) * ud->getPorosity(c);
}

void Problem::WriteBalances(double time_day, double dt){
	// --- Секции, которые будем наполнять и писать в JSON ---
	json j;
	j["time_days"] = time_day;
	j["dt_days"]   = dt;

	// ===== 1) Массы растворённых компонент =====
	json mob_mass = json::object();
	const auto & mobnames = ud->GetMobNames();
	const int n = (int)mobnames.size();

	std::vector<double> totals(n, 0.0);
	for (auto ic = mesh->BeginCell(); ic != mesh->EndCell(); ++ic) {
		const INMOST::Cell c = ic->getAsCell();
		const double Vp = CellPoreVolume(c);
		for (int i = 0; i < n; ++i) totals[i] += c->RealArray(mob_conc_tag)[i] * Vp; // масса = C * Vp
	}
	for (int i = 0; i < n; ++i) mob_mass[mobnames[i]] = totals[i];
	j["mobile_mass"] = mob_mass;

	// ===== 3) Потоки через границы по ГУ переноса =====
	json bc_flux = json::object();
	std::unordered_map<int, std::unordered_map<std::string, double>> acc_out; // из модели
	std::unordered_map<int, std::unordered_map<std::string, double>> acc_in;  // в модель

	Tag bc_transport = ud->get_tag("BC_TRANSPORT");
	Tag tag_flux = ud->get_tag("FLUX"); // объёмный расход (для Dirichlet)
	for (auto ifa = mesh->BeginFace(); ifa != mesh->EndFace(); ++ifa) {
		const INMOST::Face f = ifa->getAsFace();
		if (!f.Boundary()) continue;
		if (!tag_flux.isValid()) continue; // нет расхода — не считаем
		const int bc_id = f.Integer(bc_transport);
		const double Qf = f.Real(tag_flux); // [м^3/с], + : из модели наружу
		if (Qf >= 0){
			INMOST::Cell c = f.BackCell();
			for (size_t i = 0; i < mobnames.size(); ++i) {
				const std::string & comp = mobnames[i];
				acc_out[bc_id][comp] += Qf * c->RealArray(mob_conc_tag)[i];
			}
		}
		else{
			if (ud->getBCTransportType(bc_id) == BCType::Neumann){
				//TO DO: обработка неймановских граней будет важна для диффузии
				continue;
			}
			else if (ud->getBCTransportType(bc_id) == BCType::Dirichlet){
				for (size_t i = 0; i < mobnames.size(); ++i) {
					const std::string & comp = mobnames[i];
					acc_in[bc_id][comp] -= Qf * ud->getBCTransportValue(bc_id,comp);
				}
			}
		}
	}
	// Сформировать JSON-узел
	for (const auto & id_acc : ud->getBCTransportMap()) {//проходим по всем зарегистрированным id ГУ переноса
		int id = id_acc.first;//
		json jentry;
		json jout = json::object(), jin = json::object();
		auto it_out = acc_out.find(id);
		if (it_out != acc_out.end())
			for (const auto & name : mobnames) jout[name] = it_out->second[name];
		else
			for (const auto & name : mobnames) jout[name] = 0;
		auto it_in = acc_in.find(id);
		if (it_in != acc_in.end())
			for (const auto & name : mobnames) jin[name] = it_in->second[name];
		else
			for (const auto & name : mobnames) jin[name] = 0;
		jentry["out"] = jout;
		jentry["in"]  = jin;
		bc_flux[std::to_string(id)] = jentry;
	}
	j["boundary_flux"] = bc_flux;


	// ===== 3b) ВОДА: Потоки через границы по ГУ фильтрации =====
	// Для Dirichlet используем объёмный расход на грани из тега FLUX (м^3/с)
	// Для Neumann берём значение из BC_PRESSURE как объёмный поток на единицу длины
	// (наружу положительный) и умножаем на длина*мощность.
	json water_bc_flux = json::object();
	{
		std::unordered_map<int, double> acc_out; // из модели (+)
		std::unordered_map<int, double> acc_in;  // в модель   (+)
		for (auto ifa = mesh->BeginFace(); ifa != mesh->EndFace(); ++ifa){
			const INMOST::Face f = ifa->getAsFace();
			if (!f.Boundary()) continue;

			// геометрия (2D): длина грани * мощность слоя
			const double L = f.Area();
			double h = ud->getThickness(f.BackCell());
			const double geom = L * h;
			double Qf = 0;
			if (ud->getBCPressureType(f)==BCType::Dirichlet){
				// Нужен объёмный расход через грань
			   if (!tag_flux.isValid()) continue;
			   Qf = f.Real(tag_flux); // м^3/с, знак: + наружу

			}
			else if (ud->getBCPressureType(f)==BCType::Neumann){
				Qf = ud->getBCPressureValue(f)*geom;
			}
			int bc_id = ud->getBCPressureNum(f);
			if (Qf >= 0.0) acc_out[bc_id] += Qf; else acc_in[bc_id] += (-Qf);

		}
		// Упакуем в JSON
		for (const auto & id_acc : ud->getBCPressureMap()) {//проходим по всем зарегистрированным id ГУ фильтрации
			int id = id_acc.first;
			json jentry;
			auto it_out = acc_out.find(id);
			if (it_out != acc_out.end())
				jentry["out"] = it_out->second;
			else
				jentry["out"] = 0;
			auto it_in = acc_in.find(id);
			if (it_in != acc_in.end())
				jentry["in"] = it_in->second;
			else
				jentry["in"] = 0;
			water_bc_flux[std::to_string(id)] = jentry;
		}
	}
	j["water_boundary_flux"] = water_bc_flux;


	// ===== 4) Потоки по скважинам =====
	json wells_flux = json::object();
	const auto & wells = ud->wm.getWells();
	for (const auto & w : wells) {
		if (!w.isActive()) continue;
		json jout = json::object(), jin = json::object();
		double Q = get_well_rate(w); // + закачка, - отбор
		if (Q>=0){
			for (const auto & name : mobnames) {
				jin[name]  = Q * w.GetComponentConcentration(name);
				jout[name] = 0.0;
			}
		}
		else{
			for (const auto & name : mobnames) {
				jin[name]  = 0.0;
				Cell c = mesh->CellByLocalID(w.getLocalCellID());
				double conc = c.RealArray(mob_conc_tag)[ud->GetMobIndex(name)];
				jout[name] = -conc*Q;
			}
		}
		json jentry; jentry["in"] = jin; jentry["out"] = jout;
		wells_flux[std::to_string(w.getID())] = jentry;
	}
	j["wells_flux"] = wells_flux;

	// ===== 4b) ВОДА: Потоки по скважинам =====
	// Положительный get_well_rate(w) — закачка в модель (in), отрицательный — отбор (out).
	json water_wells_flux = json::object();
	{
		const auto & wells = ud->wm.getWells();
		for (const auto & w : wells)
		{
			if (!w.isActive()) continue;
			const double Q = get_well_rate(w); // м^3/с
			json entry;
			if (Q >= 0.0) { entry["in"] =  Q; entry["out"] = 0.0; }
			else          { entry["in"] = 0.0; entry["out"] = -Q; }
			water_wells_flux[ std::to_string(w.getID()) ] = entry;
		}
	}
	j["water_wells_flux"] = water_wells_flux;

	// ===== 5) Запись в ЕДИНЫЙ файл (JSON Lines) =====
	std::ios_base::openmode mode = std::ios::out;
	mode |= balances_initialized ? std::ios::app : std::ios::trunc; // первый раз очищаем
	std::ofstream ff(balances_file, mode);
	ff << j.dump() << '\n'; // одна запись = одна строка
	ff.close();
	balances_initialized = true;
}

double Problem::compute_max_dT(){
	double max_dT = 1e10;
	for (auto ic = mesh->BeginCell(); ic != mesh->EndCell(); ic++){
		ElementArray<Face> c_faces = ic->getFaces();
		double posflux=0e0;
		double negflux=0e0;
		Tag Flux = ud->get_tag("FLUX");
		for (ElementArray<Face>::iterator iface = c_faces.begin(); iface != c_faces.end(); ++iface)
		{
			if (iface->Real(Flux) > 0e0)
				posflux += iface->Real(Flux);
			else
				negflux += fabs(iface->Real(Flux));
		}
		double porosity = ud->getPorosity(ic->getAsCell());
		double fl = max(posflux,negflux);
		if (fl>0.0)
			max_dT = min(max_dT, ic->Volume()*ud->getThickness(ic->getAsCell())*porosity/fl);
	}
	double CFLmax = 0.5;
	return CFLmax*max_dT;
}
