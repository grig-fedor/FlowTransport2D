#include "RectangularMeshBuilder.hpp"
#include "problem.hpp"
#include "inmost.h"





//#include "analytical.hpp"


using namespace INMOST;

int main(int argc, char** argv) {
#ifdef USE_MPI
	Mesh::Initialize(&argc, &argv);
#else
	Mesh::Initialize(nullptr, nullptr);
#endif
	{
		const std::string config_file =
			(argc > 1) ? argv[1] : "model.json";

		std::ifstream f(config_file);

		if (!f.is_open()) {
			std::cerr << "Can't open configuration file: "
					  << config_file << std::endl;
			std::exit(EXIT_FAILURE);
		}
		json config = json::parse(f, nullptr, true, true); // поддержка комментариев
		int rank = 0;
#ifdef USE_MPI
		MPI_Comm_rank(INMOST_MPI_COMM_WORLD, &rank);
#endif
		if (rank == 0)
			std::cout << "[ULA] Starting mesh generation..." << std::endl;

		double Lx = require_path<double>(config,"mesh.Lx"), Ly = require_path<double>(config,"mesh.Ly");;
		int Nx = require_path<int>(config,"mesh.Nx");
		int Ny = require_path<int>(config,"mesh.Ny");
		double dx = Lx/Nx;
		double dy = Ly/Ny;

		double x0 = require_path<double>(config,"mesh.x0");
		double y0 = require_path<double>(config,"mesh.y0");

		Mesh* m = RectangularMeshBuilder::Generate(Lx, Ly, Nx, Ny, x0, y0);
		//обеспечиваем, чтобы нормали смотрели наружу BackCell
		for (Mesh::iteratorFace iface = m->BeginFace(); iface != m->EndFace(); ++iface )
			iface->FixNormalOrientation();


		if (rank == 0)
			std::cout << "[ULA] Mesh generated. Creating problem..." << std::endl;


		// распределяет сетку между процессами: каждая ячейка копируется туда, где нужна
		// Передай на соседние процессы все ячейки, которые граничат по граням (bridge = FACE), и создай призрачные копии этих ячеек (mask = CELL)
#ifdef USE_MPI
		m->ExchangeGhost(1, FACE, CELL);
#endif
		// назначает уникальные глобальные идентификаторы (для вывода и синхронизации)
		m->AssignGlobalID(CELL | FACE | NODE);

		UserData UD(m);
		UD.set_dx(dx);
		UD.set_dy(dy);
		UD.set_Lx(Lx);
		UD.set_Ly(Ly);
		UD.Init(config);

		Problem* prob = new Problem(&UD);

		prob->Save("initial.vtk");
		prob->Run();
		//Аналитика для задачи Тейса
//		Tag tag_a = m->CreateTag("Analytical", DATA_REAL, CELL, NONE, 1);
//		for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
//					Cell cell = icell->getAsCell();
//					double bar[3];
//					cell.Barycenter(bar);
//					double corX = 0.5*(Lx - dx);
//					double corY = 0.5*(Ly - dy);
//					bar[0] -= corX;
//					bar[1] -= corY;
//					cell.Real(tag_a) = theis(bar,6);
//		}
		prob->Save("result.vtk");
		delete prob;
		delete m;

		if (rank == 0)
			std::cout << "[ULA] Finished." << std::endl;
	}
	Mesh::Finalize(); // завершение INMOST (и MPI, если нужно)
	return 0;
}
