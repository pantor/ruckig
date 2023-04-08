#pragma once

#include <exception>
#include <stdexcept>
#include <string>


namespace ruckig {

struct RuckigError: public std::runtime_error {
    explicit RuckigError(const std::string& message): std::runtime_error("\n[ruckig] " + message + "\n") {}
};

}
