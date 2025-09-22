cmake_minimum_required(VERSION 3.10)
project(Uranio VERSION 1.0 LANGUAGES CXX)

# Estándar C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Incluir headers públicos
include_directories(${PROJECT_SOURCE_DIR}/include)

# Agregar todos los .cpp de src automáticamente
file(GLOB SOURCES "src/*.cpp")

# Crear biblioteca estática con todos los módulos
add_library(uranio_lib STATIC ${SOURCES})

# Crear ejecutable principal
add_executable(uranio_exec main.cpp)

# Linkear el ejecutable con la librería
target_link_libraries(uranio_exec uranio_lib)

# Opcional: definir dónde va a quedar el binario
set_target_properties(uranio_exec PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Opcional: ejecutables de tests
#add_executable(test tests/test1.cpp)
#target_link_libraries(test_old PRIVATE uranio)