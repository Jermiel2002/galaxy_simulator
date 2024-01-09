# CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

# Nom du projet
project(MyProject)

# Version du langage C++
set(CMAKE_CXX_STANDARD 11)

# Liste des fichiers sources
set(SOURCES
    src/Boite.cpp
    src/Octree.cpp
    src/Particule3D.cpp
    src/main.cpp
)

# Liste des fichiers headers
set(HEADERS
    include/Boite.hpp
    include/Octree.hpp
    include/Particule3D.hpp
)

# Ajouter l'exécutable
add_executable(my_program ${SOURCES} ${HEADERS})

# Spécifier le répertoire d'inclusion des headers
target_include_directories(my_program PRIVATE include)
