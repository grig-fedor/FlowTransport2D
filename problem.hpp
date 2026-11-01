#pragma once

#include "process.hpp"
#include "userdata.hpp"




class Problem {
public:
	explicit Problem(UserData* ud_);
	~Problem() = default;
	void Save(const std::string& filename);


	void RunSteadyFlow(bool store_fluxes = false);
	void RunUnsteadyFlow();
	void RunSteadyFlowTransport();
	void RunUnsteadyFlowTransport();
	void Run(){
		if (!ud->transport)
			if (ud->steady)
				RunSteadyFlow(true);
			else
				RunUnsteadyFlow();
		else if (ud->steady)
			RunSteadyFlowTransport();
		else
			RunUnsteadyFlowTransport();
	}
	void make_explicit_advection_step();//шагатель транспорта в старомодном стиле

	// Записать балансы за текущий шаг времени:
	//  - time_day: физическое время (сут)
	//  - dt      : шаг по времени (сут)
	void WriteBalances(double time_day, double dt);
	double compute_max_dT();

private:
	INMOST::Mesh* mesh;
	UserData* ud;
	INMOST::Automatizator aut;
	INMOST::Residual residual;
	INMOST::Sparse::Vector solution;
	INMOST::Tag pressure_tag;
	INMOST::Tag pressure_old_tag; //напор на предыдущем шаге по времени
	INMOST::Tag mob_conc_tag; // подвижные концентрации
	INMOST::Tag mob_conc_old_tag; // подвижные концентрации на предыдущем шаге по времени
	INMOST::Tag flux_tag;
	void SetTranIC();//Установка НУ для задачи переноса
	double get_well_rate(const WellData& w) const; //расчет или извлечение дебита скважины
	void run_sources(double** source, INMOST::Tag conc_prev, bool *isSource);
	void AssembleSolveUpdateFlow(process_saturatedFlow2D& Flow,
								 INMOST::dynamic_variable& varH,
								 INMOST::Solver& S,
								 INMOST::Sparse::Vector& sol, bool store_fluxes = false);
	// --- сервисные хелперы для балансов ---
	double CellBulkVolume(const INMOST::Cell & c) const;     // [м^3] = area*thickness
	double CellPoreVolume(const INMOST::Cell & c) const;     // [м^3] = bulk * porosity
	bool        balances_initialized = false;                // для очистки файла при первом вызове
	std::string balances_file        = "balances.jsonl";     // единый файл-лента для хранения балансов
};
