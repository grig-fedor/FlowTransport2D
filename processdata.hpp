#pragma once

#include "userdata.hpp"

//Класс-прокладка между UserData и дискретизацией

class ProcessData
{
protected:
	UserData *UD;
	//Название процесса
	std::string name;
public:
	ProcessData(UserData *ud) : UD(ud){};
	~ProcessData() {};
	const std::string getName() const {return name;}
	void setName(const std::string Name) {name = Name;}
	INMOST::Mesh * get_mesh() {return UD->get_mesh();};
	double getThickness(INMOST::Cell cell){return UD->getThickness(cell);}
};

//Абстрактный класс для диффузиеподобных процессов

class ProcessData_Diffusion2D : public ProcessData
{
public:
	ProcessData_Diffusion2D(UserData *ud) : ProcessData(ud) {};
	virtual ~ProcessData_Diffusion2D(){};
	//умножение на матрицу
	virtual void multA(const INMOST::Cell&, double *v, double *res) = 0;
	// тип ГУ
	virtual const BCType getBCType(const INMOST::Face &) = 0;
	// значение ГУ
	virtual double getBCValue(const INMOST::Face &) = 0;
};

//для напорной фильтрации
class ProcessData_SaturatedFlow2D : public ProcessData_Diffusion2D
{
public:
	ProcessData_SaturatedFlow2D(UserData *ud);
	~ProcessData_SaturatedFlow2D(){};
	//умножение на матрицу
	void multA(const INMOST::Cell&, double *v, double *res);
	// тип ГУ
	const BCType getBCType(const INMOST::Face &);
	// значение ГУ
	double getBCValue(const INMOST::Face &);
};
