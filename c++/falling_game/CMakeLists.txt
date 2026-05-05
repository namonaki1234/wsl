cmake_minimum_required(VERSION 3.22)
project(FallingGame)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# COMPONENTSの後に、使う機能を具体的に指定します
find_package(SFML COMPONENTS Graphics Window System CONFIG REQUIRED)

add_executable(FallingGame main.cpp)

# リンクする際もコンポーネントを指定します
target_link_libraries(FallingGame PRIVATE SFML::Graphics SFML::Window SFML::System)