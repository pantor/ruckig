# If there is no argument, create subdirectory for patch and copy all relevant files
if [ "$#" -eq 0 ]; then
    mkdir patch-cpp11
    cp -r include src examples cmake test third_party CMakeLists.txt patch-cpp11/
    cd patch-cpp11
fi

# Download optional polyfill
mkdir -p include/nonstd/
wget -O include/nonstd/optional.hpp https://raw.githubusercontent.com/martinmoene/optional-lite/master/include/nonstd/optional.hpp

# Replace function, be careful with overwriting a file
replace () {
    sed "${@:2}" "$1" > "$1.2" && mv "$1.2" "$1"
}

# Replace optional, if constexpr, and C++ version
replace include/ruckig/block.hpp -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g'
replace include/ruckig/calculator_target.hpp -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g'
replace include/ruckig/input_parameter.hpp -e 's|std::nullopt|nonstd::nullopt|g' -e 's|std::optional|nonstd::optional|g' -e 's|<optional>|<nonstd/optional.hpp>|g'
replace include/ruckig/position.hpp -e 's|<optional>|<nonstd/optional.hpp>|g'
replace include/ruckig/profile.hpp -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g'
replace include/ruckig/roots.hpp -e 's|if constexpr|if|g'
replace include/ruckig/ruckig.hpp -e 's|<optional>|<nonstd/optional.hpp>|g' -e 's|if constexpr|if|g'
replace include/ruckig/trajectory.hpp -e 's|if constexpr|if|g'
replace include/ruckig/velocity.hpp -e 's|<optional>|<nonstd/optional.hpp>|g'
replace test/test-target.cpp -e 's|std::nullopt|nonstd::nullopt|g' -e 's|<optional>|<nonstd/optional.hpp>|g'
replace test/test-target-known.cpp -e 's|std::nullopt|nonstd::nullopt|g' -e 's|<optional>|<nonstd/optional.hpp>|g'

replace CMakeLists.txt -e 's|cxx_std_17|cxx_std_11|g' -e 's|if(BUILD_BENCHMARK)|if(FALSE)|g'
