cmake -S cpp -B cpp/build
cmake --build cpp/build -j
cd cpp/build && ./main
cd ../../ && python3 python/graphs.py