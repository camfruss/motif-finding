#pragma once

#include "data.hpp"

template <typename T>
class Serial {
    public: 
        Serial(const Data& data);
        [[nodiscard]] std::vector<std::size_t> find_motifs(std::size_t k);

    private:
        const std::array<T, 4> m_background;
        const Data m_data;

        std::array<T, 4> calculate_noise(int sample_size = 100);
        std::vector<std::size_t> init_positions(std::size_t end_buffer);
        void init_pwm(int k);
        void update_pwm();
};

