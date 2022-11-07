# If there is no argument, create subdirectory for patch and copy all relevant files
if [ "$#" -eq 0 ]; then
    mkdir patch-cpp11
    cp -r include src examples cmake test third_party patch-cpp11/
    cd patch-cpp11
fi

# Download optional polyfill
mkdir -p include/nonstd/
wget -O include/nonstd/optional.hpp https://raw.githubusercontent.com/martinmoene/optional-lite/master/include/nonstd/optional.hpp

# Replace optional, if constexpr, and C++ version
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g' ../include/ruckig/block.hpp > include/ruckig/block.hpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g' ../include/ruckig/calculator_target.hpp > include/ruckig/calculator_target.hpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' ../include/ruckig/input_parameter.hpp > include/ruckig/input_parameter.hpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g' ../include/ruckig/profile.hpp > include/ruckig/profile.hpp
sed -e 's|if constexpr|if|g' ../include/ruckig/roots.hpp > include/ruckig/roots.hpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g' ../include/ruckig/ruckig.hpp > include/ruckig/ruckig.hpp
sed -e 's|if constexpr|if|g' ../include/ruckig/trajectory.hpp > include/ruckig/trajectory.hpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|<optional>|<nonstd/optional.hpp>|g' ../test/test-target.cpp > test/test-target.cpp
sed -e 's|std::nullopt|nonstd::nullopt|g' -e 's|<optional>|<nonstd/optional.hpp>|g' ../test/test-target-known.cpp > test/test-target-known.cpp

sed -e 's|cxx_std_17|cxx_std_11|g' -e 's|if(BUILD_BENCHMARK)|if(FALSE)|g' ../CMakeLists.txt > CMakeLists.txt
