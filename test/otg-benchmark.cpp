#include <random>
#include "randomizer.hpp"

#include <ruckig/parameter.hpp>
#include <ruckig/ruckig.hpp>

#ifdef WITH_REFLEXXES
#include <ruckig/alternative/reflexxes.hpp>
#endif


using namespace ruckig;


template<size_t DOFs, class OTGType>
double check_calculation(OTGType& otg, InputParameter<DOFs>& input) {
    OutputParameter<DOFs> output;
    auto result = otg.update(input, output);
    return output.calculation_duration;
}


int main() {
    std::normal_distribution<double> position_dist {0.0, 4.0};
    std::normal_distribution<double> dynamic_dist {0.0, 0.8};
    std::uniform_real_distribution<double> limit_dist {0.1, 12.0};

    {
        constexpr size_t DOFs {7};
        Ruckig<DOFs, true> otg {0.005};
        // Reflexxes<DOFs> otg {0.005};
        InputParameter<DOFs> input;
        
        srand(43);
        Randomizer<DOFs, decltype(position_dist)> p { position_dist };
        Randomizer<DOFs, decltype(dynamic_dist)> d { dynamic_dist };
        Randomizer<DOFs, decltype(limit_dist)> l { limit_dist };

        double moving_average {0.0};
        double worst {0.0};
        size_t n {1};

        for (size_t i = 0; i < 32 * 1024; ++i) {
            p.fill(input.current_position);
            d.fill_or_zero(input.current_velocity, 0.9);
            d.fill_or_zero(input.current_acceleration, 0.8);
            p.fill(input.target_position);
            d.fill_or_zero(input.target_velocity, 0.7);
            d.fill_or_zero(input.target_acceleration, 0.6);
            l.fill(input.max_velocity, input.target_velocity);
            l.fill(input.max_acceleration, input.target_acceleration);
            l.fill(input.max_jerk);

            if (!otg.validate_input(input)) {
                continue;
            }

            double time = check_calculation(otg, input);
            moving_average = moving_average + (time - moving_average) / n;
            worst = std::max(worst, time);
            ++n;
        }

        std::cout << "Average Calculation Duration 3 Dof: " << moving_average << std::endl;
        std::cout << "Worst Calculation Duration 3 Dof: " << worst << std::endl;
    }
}