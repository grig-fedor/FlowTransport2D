cmake_minimum_required(VERSION 3.10)
project(FlowTransport2D)

set(CMAKE_CXX_STANDARD 17)

# --- Опция для управления параллельной сборкой
option(PARALLEL "Enable parallel build with MPI support" OFF)

# --- Управляем INMOST и препроцессором
if(PARALLEL)
	add_definitions(-DUSE_MPI)
	set(USE_MPI ON CACHE BOOL "Enable MPI in INMOST" FORCE)
	set(USE_MPI_P2P ON CACHE BOOL "Enable MPI P2P in INMOST" FORCE)
	set(USE_MPI_FILE ON CACHE BOOL "Enable MPI file IO in INMOST" FORCE)
	set(USE_MPI2 ON CACHE BOOL "Enable MPI-2 features in INMOST" FORCE)

	if(WIN32)
		# === Указание путей к MS-MPI ===
		set(MSMPI_INCLUDE "C:/Program Files (x86)/Microsoft SDKs/MPI/Include")
		set(MSMPI_LIB "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64")

		include_directories(${MSMPI_INCLUDE})
		link_directories(${MSMPI_LIB})
		link_libraries(msmpi)
	endif()
else()
	set(USE_MPI OFF CACHE BOOL "Disable MPI in INMOST" FORCE)
	set(USE_MPI_P2P OFF CACHE BOOL "Disable MPI P2P in INMOST" FORCE)
	set(USE_MPI_FILE OFF CACHE BOOL "Disable MPI file IO in INMOST" FORCE)
	set(USE_MPI2 OFF CACHE BOOL "Disable MPI-2 features in INMOST" FORCE)
endif()

# === Подключение INMOST ===
add_subdirectory(external/inmost)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/inmost/Source/Headers)

# === JSON header-only library ===
include_directories(${CMAKE_SOURCE_DIR}/external)

# === Подключение собственного кода ===
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(FlowTransport2D
  main.cpp
  problem.cpp
  RectangularMeshBuilder.cpp
  process.cpp
  fvDiffusion2D.cpp
  userdata.cpp
  processdata.cpp
  WellData.cpp
  WellManager.cpp
# analytical.cpp
)

target_link_libraries(FlowTransport2D inmost)
