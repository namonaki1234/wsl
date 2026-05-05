cd ~/wsl/c++/falling_game/build
rm -rf *  # 一旦きれいにリセット
cmake .. -DCMAKE_TOOLCHAIN_FILE=/root/vcpkg/scripts/buildsystems/vcpkg.cmake
make
./FallingGame