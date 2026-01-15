#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "../core/StateSpace.hpp"
#include "../core/CostLandscape.hpp"

namespace UE {

/**
 * MetastableClass - Clase Metaestable A_k
 *
 * Representa una cuenca/atractor en el paisaje condicionado.
 * Validado dinámicamente por: τ_relax(A_k) << E[τ_exit(A_k)]
 *
 * Ref: UE Framework Sección 5 - Clases y Dinámica
 */
template<typename ConfigType = std::vector<double>>
struct MetastableClass {
    int id;                      // Identificador de la clase
    ConfigType center;           // Centro/atractor de la clase
    double probability;          // p_k = probabilidad de la clase
    double cost;                 // F_k = costo en el centro

    // Propiedades dinámicas (se llenan después de validación)
    double tau_relax;            // Tiempo de relajación
    double tau_exit_mean;        // Tiempo medio de salida
    double tau_exit_std;         // Desviación de τ_exit
    double metastability_ratio;  // τ_exit / τ_relax

    // Región de la clase (límites)
    ConfigType lower_bound;
    ConfigType upper_bound;

    // Estado de validación
    bool is_validated;
    std::string validation_status;

    MetastableClass()
        : id(-1)
        , probability(0)
        , cost(0)
        , tau_relax(0)
        , tau_exit_mean(0)
        , tau_exit_std(0)
        , metastability_ratio(0)
        , is_validated(false)
    {}
};

/**
 * MetastableClassDetector - Detector de Clases Metaestables
 *
 * Identifica y valida clases metaestables A_k a partir de:
 * 1. Mínimos locales del paisaje F_n^R
 * 2. Validación dinámica: τ_relax << τ_exit
 * 3. Probabilidad alta de permanencia durante τ_p
 *
 * Ref: UE Framework Sección 5.1 - Identificación de clases
 */
template<typename ConfigType = std::vector<double>>
class MetastableClassDetector {
public:
    using Config = ConfigType;
    using Class = MetastableClass<Config>;

private:
    const StateSpace<Config>* state_space_;
    const CostLandscape<Config>* landscape_;

    std::vector<Class> classes_;

    // Parámetros de detección
    double min_probability_ = 0.01;        // Probabilidad mínima para ser clase
    double min_metastability_ratio_ = 10.0; // τ_exit/τ_relax mínimo

public:
    MetastableClassDetector(
        const StateSpace<Config>& space,
        const CostLandscape<Config>& landscape)
        : state_space_(&space)
        , landscape_(&landscape)
    {}

    // Detectar clases a partir de mínimos del paisaje
    void detectFromLandscape() {
        classes_.clear();

        auto minima = landscape_->getMinima();
        int id = 0;

        for (const auto& min : minima) {
            if (min.probability >= min_probability_) {
                Class c;
                c.id = id++;
                c.center = min.config;
                c.probability = min.probability;
                c.cost = min.cost;
                c.is_validated = false;
                c.validation_status = "Pendiente validación dinámica";

                // Estimar región de la clase (simplificado)
                c.lower_bound = c.center;
                c.upper_bound = c.center;
                double delta = state_space_->getLevel(0).delta_n * 2.0;
                for (size_t d = 0; d < c.center.size(); ++d) {
                    c.lower_bound[d] -= delta;
                    c.upper_bound[d] += delta;
                }

                classes_.push_back(c);
            }
        }

        // Ordenar por probabilidad (mayor primero)
        std::sort(classes_.begin(), classes_.end(),
            [](const Class& a, const Class& b) {
                return a.probability > b.probability;
            });

        // Re-asignar IDs después de ordenar
        for (size_t i = 0; i < classes_.size(); ++i) {
            classes_[i].id = static_cast<int>(i);
        }
    }

    // Verificar si una configuración pertenece a una clase
    int classifyConfig(const Config& config) const {
        for (const auto& c : classes_) {
            bool inside = true;
            for (size_t d = 0; d < config.size(); ++d) {
                if (config[d] < c.lower_bound[d] || config[d] > c.upper_bound[d]) {
                    inside = false;
                    break;
                }
            }
            if (inside) return c.id;
        }
        return -1;  // No pertenece a ninguna clase
    }

    // Validar clase con tiempos dinámicos
    void validateClass(int class_id, double tau_relax, double tau_exit_mean, double tau_exit_std) {
        if (class_id < 0 || class_id >= static_cast<int>(classes_.size())) return;

        auto& c = classes_[class_id];
        c.tau_relax = tau_relax;
        c.tau_exit_mean = tau_exit_mean;
        c.tau_exit_std = tau_exit_std;

        if (tau_relax > 0) {
            c.metastability_ratio = tau_exit_mean / tau_relax;
        }

        // Criterio de validación: τ_relax << τ_exit
        if (c.metastability_ratio >= min_metastability_ratio_) {
            c.is_validated = true;
            c.validation_status = "Validada: τ_exit/τ_relax = " +
                std::to_string(c.metastability_ratio);
        } else {
            c.is_validated = false;
            c.validation_status = "No validada: τ_exit/τ_relax = " +
                std::to_string(c.metastability_ratio) + " < " +
                std::to_string(min_metastability_ratio_);
        }
    }

    // Obtener clases detectadas
    const std::vector<Class>& getClasses() const { return classes_; }

    // Obtener solo clases validadas
    std::vector<Class> getValidatedClasses() const {
        std::vector<Class> validated;
        for (const auto& c : classes_) {
            if (c.is_validated) {
                validated.push_back(c);
            }
        }
        return validated;
    }

    // Calcular medidas de dominancia
    struct DominanceMetrics {
        double D_H;     // p_1/p_2 (ratio de dominancia)
        double Delta;   // p_1 - p_2 (diferencia)
        double p_1;     // Probabilidad de la clase dominante
        double p_2;     // Probabilidad de la segunda clase
        int dominant_class_id;
    };

    DominanceMetrics calculateDominance() const {
        DominanceMetrics m{};

        if (classes_.size() < 2) {
            m.D_H = std::numeric_limits<double>::infinity();
            m.Delta = 1.0;
            m.p_1 = classes_.empty() ? 0.0 : classes_[0].probability;
            m.p_2 = 0.0;
            m.dominant_class_id = classes_.empty() ? -1 : classes_[0].id;
            return m;
        }

        m.p_1 = classes_[0].probability;
        m.p_2 = classes_[1].probability;
        m.D_H = m.p_1 / std::max(m.p_2, 1e-10);
        m.Delta = m.p_1 - m.p_2;
        m.dominant_class_id = classes_[0].id;

        return m;
    }

    // Número de clases
    size_t numClasses() const { return classes_.size(); }

    // Configuración
    void setMinProbability(double p) { min_probability_ = p; }
    void setMinMetastabilityRatio(double r) { min_metastability_ratio_ = r; }
};

} // namespace UE
