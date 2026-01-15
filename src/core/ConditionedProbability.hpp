#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "StateSpace.hpp"
#include "ERSCondition.hpp"

namespace UE {

/**
 * ConditionedProbability - Calculador de Probabilidad Condicionada p_n(ψ|R)
 *
 * Calcula la masa de probabilidad condicionada:
 * p_n(ψ|R) = P(d_n(φ_n, ψ) ≤ δ_n | R)
 *
 * Es decir, la probabilidad de estar cerca de ψ dado que R se sostiene.
 *
 * Métodos de estimación:
 * 1. Muestreo directo de trayectorias condicionadas
 * 2. Histograma en espacio de configuraciones
 * 3. Kernel density estimation
 *
 * Ref: UE Framework Sección 4.3 - Masa de probabilidad condicionada
 */
template<typename ConfigType = std::vector<double>>
class ConditionedProbability {
public:
    using Config = ConfigType;
    using Trajectory = std::vector<Config>;

    struct EstimationResult {
        double probability;       // p_n(ψ|R)
        double standard_error;    // Error estándar
        int samples_used;         // Muestras usadas
        int samples_satisfying_R; // Muestras que cumplen R
    };

private:
    const StateSpace<Config>* state_space_;
    const ERSCondition<Config>* ers_condition_;
    int coarse_grain_level_;

    // Trayectorias condicionadas almacenadas
    std::vector<Trajectory> conditioned_trajectories_;

    // Configuraciones visitadas bajo R
    std::vector<Config> conditioned_samples_;

public:
    ConditionedProbability(
        const StateSpace<Config>& space,
        const ERSCondition<Config>& ers,
        int level = 0)
        : state_space_(&space)
        , ers_condition_(&ers)
        , coarse_grain_level_(level)
    {}

    // Agregar trayectoria condicionada (que cumple R)
    void addConditionedTrajectory(const Trajectory& traj) {
        if (ers_condition_->checkTrajectory(traj)) {
            conditioned_trajectories_.push_back(traj);
            // Agregar configuraciones individuales
            for (const auto& config : traj) {
                conditioned_samples_.push_back(config);
            }
        }
    }

    // Agregar múltiples trayectorias, filtrando las que cumplen R
    void addTrajectories(const std::vector<Trajectory>& trajectories) {
        for (const auto& traj : trajectories) {
            addConditionedTrajectory(traj);
        }
    }

    // Calcular p_n(ψ|R) para una configuración específica
    EstimationResult calculate(const Config& psi) const {
        if (conditioned_samples_.empty()) {
            return {0.0, 0.0, 0, 0};
        }

        const auto& level = state_space_->getLevel(coarse_grain_level_);
        double delta_n = level.delta_n;

        // Contar muestras dentro de la bola B(ψ, δ_n)
        int count = 0;
        for (const auto& sample : conditioned_samples_) {
            if (state_space_->distance(sample, psi, coarse_grain_level_) <= delta_n) {
                ++count;
            }
        }

        int N = static_cast<int>(conditioned_samples_.size());
        double p = static_cast<double>(count) / N;

        // Error estándar binomial
        double se = std::sqrt(p * (1.0 - p) / N);

        return {p, se, N, static_cast<int>(conditioned_trajectories_.size())};
    }

    // Calcular p_n(ψ|R) para un grid de configuraciones (1D o 2D)
    std::vector<std::pair<Config, double>> calculateGrid(
        const std::vector<Config>& grid_points) const
    {
        std::vector<std::pair<Config, double>> result;
        for (const auto& psi : grid_points) {
            auto est = calculate(psi);
            result.emplace_back(psi, est.probability);
        }
        return result;
    }

    // Calcular histograma de probabilidades
    std::map<int, double> calculateHistogram(int num_bins) const {
        if (conditioned_samples_.empty() || state_space_->dimension() != 1) {
            return {};
        }

        double lower = state_space_->lowerBounds()[0];
        double upper = state_space_->upperBounds()[0];
        double bin_width = (upper - lower) / num_bins;

        std::map<int, double> histogram;
        for (int i = 0; i < num_bins; ++i) {
            histogram[i] = 0.0;
        }

        for (const auto& sample : conditioned_samples_) {
            int bin = static_cast<int>((sample[0] - lower) / bin_width);
            bin = std::max(0, std::min(num_bins - 1, bin));
            histogram[bin] += 1.0;
        }

        // Normalizar
        double total = static_cast<double>(conditioned_samples_.size());
        for (auto& [bin, count] : histogram) {
            count /= total;
        }

        return histogram;
    }

    // Encontrar los máximos de p_n(ψ|R) - candidatos a clases metaestables
    std::vector<Config> findModes(int num_bins = 50) const {
        std::vector<Config> modes;

        if (conditioned_samples_.empty()) return modes;

        // Para 1D, usar histograma
        if (state_space_->dimension() == 1) {
            auto hist = calculateHistogram(num_bins);
            double lower = state_space_->lowerBounds()[0];
            double upper = state_space_->upperBounds()[0];
            double bin_width = (upper - lower) / num_bins;

            // Encontrar máximos locales
            std::vector<int> bins_sorted;
            for (const auto& [bin, _] : hist) {
                bins_sorted.push_back(bin);
            }

            for (int bin : bins_sorted) {
                double p = hist[bin];
                double p_left = (bin > 0) ? hist[bin - 1] : 0.0;
                double p_right = (bin < num_bins - 1) ? hist[bin + 1] : 0.0;

                if (p > p_left && p > p_right && p > 0.01) {
                    Config mode(1);
                    mode[0] = lower + (bin + 0.5) * bin_width;
                    modes.push_back(mode);
                }
            }
        }
        // Para dimensiones mayores, usar clustering (simplificado)
        else {
            // K-means simplificado con los datos
            // Por ahora, retornar el centroide como modo único
            Config centroid(state_space_->dimension(), 0.0);
            for (const auto& sample : conditioned_samples_) {
                for (int d = 0; d < state_space_->dimension(); ++d) {
                    centroid[d] += sample[d];
                }
            }
            for (int d = 0; d < state_space_->dimension(); ++d) {
                centroid[d] /= conditioned_samples_.size();
            }
            modes.push_back(centroid);
        }

        return modes;
    }

    // Calcular probabilidades de las clases (p_k)
    std::vector<double> calculateClassProbabilities(
        const std::vector<Config>& class_centers) const
    {
        std::vector<double> p_k(class_centers.size(), 0.0);

        if (conditioned_samples_.empty()) return p_k;

        // Asignar cada muestra a la clase más cercana
        for (const auto& sample : conditioned_samples_) {
            int closest = 0;
            double min_dist = state_space_->distance(sample, class_centers[0]);

            for (size_t k = 1; k < class_centers.size(); ++k) {
                double d = state_space_->distance(sample, class_centers[k]);
                if (d < min_dist) {
                    min_dist = d;
                    closest = static_cast<int>(k);
                }
            }
            p_k[closest] += 1.0;
        }

        // Normalizar
        double total = static_cast<double>(conditioned_samples_.size());
        for (auto& p : p_k) {
            p /= total;
        }

        return p_k;
    }

    // Getters
    size_t numConditionedSamples() const { return conditioned_samples_.size(); }
    size_t numConditionedTrajectories() const { return conditioned_trajectories_.size(); }
    const std::vector<Config>& getSamples() const { return conditioned_samples_; }

    // Limpiar datos
    void clear() {
        conditioned_trajectories_.clear();
        conditioned_samples_.clear();
    }
};

} // namespace UE
