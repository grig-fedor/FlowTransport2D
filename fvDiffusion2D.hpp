#pragma once
#include "inmost.h"
#include "processdata.hpp"



class FV_Diffusion2D
{
protected:
	INMOST::Mesh *mesh;
	INMOST::Tag tagT; //тег, хранящий проводимость грани
	ProcessData_Diffusion2D *processdata;
public:
	FV_Diffusion2D(ProcessData_Diffusion2D *pd);
	void build();
	//получение потока через грань
	//--Разница между dynamic_variable и variable:
	//---dynamic_variable нужен как источник информации о переменной (в т.ч. ее индекс в глобальном векторе)
	//---variable нужен как участник выражения
	INMOST::variable getDgradH(const INMOST::Face &f, INMOST::dynamic_variable &H);

};
