#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include "StateSpace.hpp"
#include "ConditionedProbability.hpp"

namespace UE {

/**
 * CostLandscape - Paisaje de Costo Condicionado F_n^R
 *
 * Implementa el paisaje de costo:
 * F_n^R[ψ] = -ε_n log p_n(ψ|R)
 *
 * Propiedades:
 * - Mínimos locales → Clases metaestables A_k
 * - Barreras → Costos de transición
 * - ε_n → Temperatura efectiva (resolución)
 *
 * Ref: UE Framework Sección 4.4 - Paisaje condicionado
 */
template<typename ConfigType = std::vector<double>>
class CostLandscape {
public:
    using Config = ConfigType;

    struct LandscapePoint {
        Config config;
        double probability;  // p_n(ψ|R)
        double cost;         // F_n^R[ψ]
        bool is_minimum;     // ¿Es mínimo local?
        bool is_maximum;     // ¿Es máximo local (barrera)?
    };

private:
    const StateSpace<Config>* state_space_;
    const ConditionedProbability<Config>* cond_prob_;
    int level_idx_;
    double epsilon_n_;

    // Paisaje discretizado
    std::vector<LandscapePoint> landscape_;

    // Regularización para evitar log(0)
    double min_probability_ = 1e-10;

public:
    CostLandscape(
        const StateSpace<Config>& space,
        const ConditionedProbability<Config>& prob,
        int level = 0)
        : state_space_(&space)
        , cond_prob_(&prob)
        , level_idx_(level)
        , epsilon_n_(space.getLevel(level).epsilon_n)
    {}

    // Calcular F_n^R para una configuración
    double calculate(const Config& psi) const {
        auto result = cond_prob_->calculate(psi);
        double p = std::max(result.probability, min_probability_);
        return -epsilon_n_ * std::log(p);
    }

    // Construir paisaje sobre un grid 1D
    void buildGrid1D(int num_points) {
        if (state_space_->dimension() != 1) return;

        landscape_.clear();
        double lower = state_space_->lowerBounds()[0];
        double upper = state_space_->upperBounds()[0];
        double step = (upper - lower) / (num_points - 1);

        for (int i = 0; i < num_points; ++i) {
            Config psi(1);
            psi[0] = lower + i * step;

            auto result = cond_prob_->calculate(psi);
            double p = std::max(result.probability, min_probability_);
            double F = -epsilon_n_ * std::log(p);

            landscape_.push_back({psi, p, F, false, false});
        }

        // Identificar mínimos y máximos locales
        identifyExtrema();
    }

    // Construir paisaje sobre un grid 2D
    void buildGrid2D(int points_per_dim) {
        if (state_space_->dimension() != 2) return;

        landscape_.clear();
        double lower_x = state_space_->lowerBounds()[0];
        double upper_x = state_space_->upperBounds()[0];
        double lower_y = state_space_->lowerBounds()[1];
        double upper_y = state_space_->upperBounds()[1];

        double step_x = (upper_x - lower_x) / (points_per_dim - 1);
        double step_y = (upper_y - lower_y) / (points_per_dim - 1);

        for (int i = 0; i < points_per_dim; ++i) {
            for (int j = 0; j < points_per_dim; ++j) {
                Config psi(2);
                psi[0] = lower_x + i * step_x;
                psi[1] = lower_y + j * step_y;

                auto result = cond_prob_->calculate(psi);
                double p = std::max(result.probability, min_probability_);
                double F = -epsilon_n_ * std::log(p);

                landscape_.push_back({psi, p, F, false, false});
            }
        }
    }

    // Identificar extremos locales (mínimos = cuencas, máximos = barreras)
    void identifyExtrema() {
        if (landscape_.size() < 3) return;

        // Para 1D, simple comparación con vecinos
        if (state_space_->dimension() == 1) {
            for (size_t i = 1; i < landscape_.size() - 1; ++i) {
                double F = landscape_[i].cost;
                double F_left = landscape_[i - 1].cost;
                double F_right = landscape_[i + 1].cost;

                landscape_[i].is_minimum = (F < F_left && F < F_right);
                landscape_[i].is_maximum = (F > F_left && F > F_right);
            }
        }
    }

    // Obtener mínimos locales (candidatos a clases metaestables)
    std::vector<LandscapePoint> getMinima() const {
        std::vector<LandscapePoint> minima;
        for (const auto& point : landscape_) {
            if (point.is_minimum) {
                minima.push_back(point);
            }
        }
        // Ordenar por costo (menor primero)
        std::sort(minima.begin(), minima.end(),
            [](const LandscapePoint& a, const LandscapePoint& b) {
                return a.cost < b.cost;
            });
        return minima;
    }

    // Obtener máximos locales (barreras)
    std::vector<LandscapePoint> getMaxima() const {
        std::vector<LandscapePoint> maxima;
        for (const auto& point : landscape_) {
            if (point.is_maximum) {
                maxima.push_back(point);
            }
        }
        return maxima;
    }

    // Calcular barrera entre dos mínimos
    double calculateBarrier(const Config& min1, const Config& min2) const {
        // Encontrar el máximo de F en el camino entre min1 y min2
        // (Para 1D, es el máximo entre las posiciones)

        if (state_space_->dimension() != 1 || landscape_.empty()) {
            return std::numeric_limits<double>::infinity();
        }

        double x1 = min1[0];
        double x2 = min2[0];
        if (x1 > x2) std::swap(x1, x2);

        double max_cost = 0.0;
        double min_cost_1 = std::numeric_limits<double>::infinity();
        double min_cost_2 = std::numeric_limits<double>::infinity();

        for (const auto& point : landscape_) {
            double x = point.config[0];
            if (x >= x1 && x <= x2) {
                max_cost = std::max(max_cost, point.cost);
            }
            // Encontrar costo en los mínimos
            if (std::abs(x - min1[0]) < 0.1) {
                min_cost_1 = std::min(min_cost_1, point.cost);
            }
            if (std::abs(x - min2[0]) < 0.1) {
                min_cost_2 = std::min(min_cost_2, point.cost);
            }
        }

        // Barrera desde el mínimo más bajo
        double base_cost = std::min(min_cost_1, min_cost_2);
        return max_cost - base_cost;
    }

    // Calcular dominancia D_H = p_1/p_2
    double calculateDominanceRatio() const {
        auto minima = getMinima();
        if (minima.size() < 2) return std::numeric_limits<double>::infinity();

        double p1 = minima[0].probability;
        double p2 = minima[1].probability;

        return p1 / std::max(p2, min_probability_);
    }

    // Calcular dominancia Δ = p_1 - p_2
    double calculateDominanceDifference() const {
        auto minima = getMinima();
        if (minima.size() < 2) return 1.0;

        double p1 = minima[0].probability;
        double p2 = minima[1].probability;

        return p1 - p2;
    }

    // Getters
    const std::vector<LandscapePoint>& getLandscape() const { return landscape_; }
    double epsilon() const { return epsilon_n_; }

    void setEpsilon(double eps) { epsilon_n_ = eps; }
    void setMinProbability(double p) { min_probability_ = p; }
};

} // namespace UE
