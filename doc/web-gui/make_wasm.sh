emcc src/ruckig/wasm.cpp src/ruckig/brake.cpp \
    src/ruckig/position_first_step1.cpp src/ruckig/position_first_step2.cpp \
    src/ruckig/position_second_step1.cpp src/ruckig/position_second_step2.cpp \
    src/ruckig/position_third_step1.cpp src/ruckig/position_third_step2.cpp \
    src/ruckig/velocity_second_step1.cpp src/ruckig/velocity_second_step2.cpp \
    src/ruckig/velocity_third_step1.cpp src/ruckig/velocity_third_step2.cpp \
    -Iinclude -Ithird_party \
    -std=c++17 -lembind -Os \
    -s MODULARIZE -s EXPORT_ES6 -s EXPORT_NAME='RuckigModule' -s ENVIRONMENT='web' -s EXPORTED_RUNTIME_METHODS=ccall,cwrap \
    -o doc/web-gui/src/ruckig.js