# If there is no argument, create subdirectory for patch and copy all relevant files
if [ "$#" -eq 0 ]; then
    mkdir patch-cpp11
    cp -r include src examples cmake test CMakeLists.txt patch-cpp11/
    cd patch-cpp11
fi

# Download optional polyfill
mkdir -p include/nonstd/
wget -O include/nonstd/optional.hpp https://raw.githubusercontent.com/martinmoene/optional-lite/master/include/nonstd/optional.hpp
wget -O include/nonstd/variant.hpp https://raw.githubusercontent.com/martinmoene/variant-lite/master/include/nonstd/variant.hpp

# Replace optional, if constexpr, and C++ version
sed -i -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g' include/ruckig/*.hpp
sed -i -e 's|std::nullopt|nonstd::nullopt|g' test/otg-*.cpp
sed -i -e 's|std::variant|nonstd::variant|g' -e 's|std::visit|nonstd::visit|g' -e 's|std::holds_alternative|nonstd::holds_alternative|g' -e 's|<variant>|<nonstd/variant.hpp>|g' include/ruckig/*.hpp
sed -i -e 's|cxx_std_17|cxx_std_11|g' -e 's|if(BUILD_PYTHON_MODULE)|if(FALSE)|g' -e 's|if(BUILD_BENCHMARK)|if(FALSE)|g' -e 's| /WX||g' -e 's| -Werror||g' CMakeLists.txt
