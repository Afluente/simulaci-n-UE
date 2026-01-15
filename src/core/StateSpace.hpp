#pragma once
#include <vector>
#include <cmath>
#include <functional>
#include <random>

namespace UE {

/**
 * StateSpace - Espacio de Configuraciones con Coarse-Graining
 *
 * Implementa el espacio de estados φ con:
 * - Escala n (nivel de coarse-graining)
 * - Métrica d_n (distancia entre configuraciones)
 * - Tolerancia δ_n (resolución efectiva)
 *
 * Ref: UE Framework Sección 4 - Selección Condicionada
 */
template<typename ConfigType = std::vector<double>>
class StateSpace {
public:
    using Config = ConfigType;
    using DistanceFunc = std::function<double(const Config&, const Config&)>;

    struct CoarseGrainLevel {
        int n;              // Nivel de escala
        double delta_n;     // Tolerancia δ_n
        double epsilon_n;   // Parámetro de temperatura ε_n = 1/β_n
    };

private:
    int dimension_;
    std::vector<CoarseGrainLevel> levels_;
    DistanceFunc distance_func_;

    // Límites del espacio de configuraciones
    std::vector<double> lower_bounds_;
    std::vector<double> upper_bounds_;

    // Generador aleatorio
    mutable std::mt19937 rng_;

public:
    StateSpace(int dim, unsigned seed = 42)
        : dimension_(dim)
        , rng_(seed)
    {
        lower_bounds_.resize(dim, -10.0);
        upper_bounds_.resize(dim, 10.0);

        // Métrica por defecto: distancia euclídea
        distance_func_ = [](const Config& a, const Config& b) {
            double sum = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                double diff = a[i] - b[i];
                sum += diff * diff;
            }
            return std::sqrt(sum);
        };

        // Nivel de coarse-graining por defecto
        levels_.push_back({1, 0.5, 1.0});
    }

    // Configurar límites del espacio
    void setBounds(const std::vector<double>& lower, const std::vector<double>& upper) {
        lower_bounds_ = lower;
        upper_bounds_ = upper;
    }

    // Configurar métrica personalizada
    void setDistanceFunction(DistanceFunc func) {
        distance_func_ = std::move(func);
    }

    // Añadir nivel de coarse-graining
    void addCoarseGrainLevel(int n, double delta_n, double epsilon_n) {
        levels_.push_back({n, delta_n, epsilon_n});
    }

    // Obtener nivel de coarse-graining
    const CoarseGrainLevel& getLevel(size_t idx) const {
        return levels_.at(idx);
    }

    size_t numLevels() const { return levels_.size(); }
    int dimension() const { return dimension_; }

    // Distancia entre configuraciones en escala n
    double distance(const Config& a, const Config& b, int level_idx = 0) const {
        // La métrica puede depender de n (coarse-graining)
        double d = distance_func_(a, b);
        // Escalar por el nivel si es necesario
        return d;
    }

    // Verificar si dos configuraciones son equivalentes en escala n
    // d_n(φ, ψ) ≤ δ_n
    bool areEquivalent(const Config& a, const Config& b, int level_idx = 0) const {
        double d = distance(a, b, level_idx);
        return d <= levels_[level_idx].delta_n;
    }

    // Generar configuración aleatoria uniforme
    Config randomConfig() const {
        Config config(dimension_);
        std::uniform_real_distribution<double> dist;
        for (int i = 0; i < dimension_; ++i) {
            dist = std::uniform_real_distribution<double>(lower_bounds_[i], upper_bounds_[i]);
            config[i] = dist(rng_);
        }
        return config;
    }

    // Generar configuración con perturbación gaussiana
    Config perturbConfig(const Config& base, double sigma) const {
        Config config = base;
        std::normal_distribution<double> noise(0.0, sigma);
        for (int i = 0; i < dimension_; ++i) {
            config[i] += noise(rng_);
            // Aplicar límites
            config[i] = std::max(lower_bounds_[i], std::min(upper_bounds_[i], config[i]));
        }
        return config;
    }

    const std::vector<double>& lowerBounds() const { return lower_bounds_; }
    const std::vector<double>& upperBounds() const { return upper_bounds_; }

    std::mt19937& rng() const { return rng_; }
};

} // namespace UE
