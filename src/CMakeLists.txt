# Add lab 1 executable target
add_executable(
    lab1
    "${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab1.cxx"
)

# Add lab 2 executable target
add_executable(
    lab2
    "${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/SearchAlgorithms.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab2.cxx"
)

# Add lab 3 executable target
add_executable(
    lab3
    "${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/SearchAlgorithms.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/GradientAlgorithms.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab3.cxx"
)

# Add lab 4 executable target
add_executable(
    lab4
    "${CMAKE_CURRENT_SOURCE_DIR}/GeneticAlgorithm.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab4.cxx"
)

# Add lab 4 (3) executable target
add_executable(
    lab4_3
    "${CMAKE_CURRENT_SOURCE_DIR}/GeneticAlgorithm.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab4_3.cxx"
)

# Add lab 4 (4) executable target
add_executable(
    lab4_4
    "${CMAKE_CURRENT_SOURCE_DIR}/GeneticAlgorithm.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab4_4.cxx"
)

# Add lab 4 (5) executable target
add_executable(
    lab4_5
    "${CMAKE_CURRENT_SOURCE_DIR}/GeneticAlgorithm.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab4_5.cxx"
)

# Add lab 5 executable target
add_executable(
    lab5
    "${CMAKE_CURRENT_SOURCE_DIR}/Matrix.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/NumInt.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/Lab5.cxx"
)
