#include "RectangularMeshBuilder.hpp"
#ifdef USE_MPI
#include <mpi.h>
#endif
using namespace INMOST;

Mesh* RectangularMeshBuilder::Generate(double Lx, double Ly, int Nx, int Ny, double x0, double y0) {
	Mesh* m = new Mesh();

	//Привязывает сетку m к текущему MPI-контексту, чтобы в будущем можно было делать параллельные вычисления (если нужно)
	m->SetCommunicator(INMOST_MPI_COMM_WORLD);

	//устанавливает опцию вывода (verbosity) для операций ввода/вывода файла (Load, Save) в INMOST 0 = молчать
	m->SetFileOption("VERBOSITY", "0");

	m->SetFileOption("GRID_DIMENSION", "2");

	int rank = 0;
#ifdef USE_MPI
	MPI_Comm_rank(INMOST_MPI_COMM_WORLD, &rank);
#endif

	if (rank == 0) {
		double dx = Lx / Nx;
		double dy = Ly / Ny;

		std::vector<Node> all_nodes;
		all_nodes.reserve((Nx + 1) * (Ny + 1));

		// генерируем узлы
		for (int j = 0; j <= Ny; ++j) {
			for (int i = 0; i <= Nx; ++i) {
				double coords[3] = {x0 + i * dx, y0 + j * dy, 0.0};
				all_nodes.push_back(m->CreateNode(coords));
			}
		}

		// создаем ячейки (только на 0-м процессе)
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				int n0 = j * (Nx + 1) + i;
				int n1 = n0 + 1;
				int n2 = n1 + (Nx + 1);
				int n3 = n0 + (Nx + 1);

				ElementArray<Node> face_nodes(m);
				ElementArray<Face> cell_faces(m);

				ElementArray<Node> cell_nodes(m);
				cell_nodes.push_back(all_nodes[n0]);
				cell_nodes.push_back(all_nodes[n1]);
				cell_nodes.push_back(all_nodes[n2]);
				cell_nodes.push_back(all_nodes[n3]);

				// левая грань
				face_nodes.clear();
				face_nodes.push_back(all_nodes[n0]);
				face_nodes.push_back(all_nodes[n3]);
				cell_faces.push_back(m->CreateFace(face_nodes).first);

				// верхняя грань
				face_nodes.clear();
				face_nodes.push_back(all_nodes[n3]);
				face_nodes.push_back(all_nodes[n2]);
				cell_faces.push_back(m->CreateFace(face_nodes).first);

				// правая грань
				face_nodes.clear();
				face_nodes.push_back(all_nodes[n2]);
				face_nodes.push_back(all_nodes[n1]);
				cell_faces.push_back(m->CreateFace(face_nodes).first);

				// нижняя грань
				face_nodes.clear();
				face_nodes.push_back(all_nodes[n1]);
				face_nodes.push_back(all_nodes[n0]);
				cell_faces.push_back(m->CreateFace(face_nodes).first);

				m->CreateCell(cell_faces,cell_nodes);
			}
		}
	}
	//m->ReorderEmpty(NODE | EDGE | FACE | CELL);
//---!db
//	for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
//			Cell cell = icell->getAsCell();
//			std::cout << "V = "<<cell.Volume()<<std::endl;
//			}

	return m;
}
