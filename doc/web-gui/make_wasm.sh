emcc src/ruckig/wasm.cpp src/ruckig/brake.cpp \
    src/ruckig/position-first-step1.cpp src/ruckig/position-first-step2.cpp src/ruckig/position-second-step1.cpp src/ruckig/position-second-step2.cpp src/ruckig/position-third-step1.cpp src/ruckig/position-third-step2.cpp \
    src/ruckig/velocity-second-step1.cpp src/ruckig/velocity-second-step2.cpp src/ruckig/velocity-third-step1.cpp src/ruckig/velocity-third-step2.cpp \
    -Iinclude -Ithird_party \
    -std=c++17 -lembind -Os \
    -s MODULARIZE -s EXPORT_ES6 -s EXPORT_NAME='RuckigModule' -s ENVIRONMENT='web' -s EXPORTED_RUNTIME_METHODS=ccall,cwrap \
    -o doc/web-gui/src/ruckig.js