#pragma once
#include <vector>
#include <functional>
#include <deque>
#include <cmath>

namespace UE {

/**
 * ERSCondition - Evento Raro Sostenido (Sustained Rare Event)
 *
 * Define la condición R que caracteriza el evento raro sostenido:
 * R := {f(φ_t) > f_c durante t ∈ [0, τ_p]}
 *
 * Componentes:
 * - f: Observable que define el evento
 * - f_c: Umbral crítico
 * - τ_p: Persistencia mínima
 * - P(R): Probabilidad del evento (debe ser rara)
 *
 * Ref: UE Framework Sección 4.2 - Definición de R
 */
template<typename ConfigType = std::vector<double>>
class ERSCondition {
public:
    using Config = ConfigType;
    using Observable = std::function<double(const Config&)>;

private:
    Observable observable_;      // f(φ)
    double threshold_;           // f_c
    int persistence_;            // τ_p (en pasos discretos)
    bool above_threshold_;       // true si R es f > f_c, false si f < f_c

    // Buffer para verificar persistencia
    std::deque<double> history_;

public:
    ERSCondition(Observable obs, double threshold, int persistence, bool above = true)
        : observable_(std::move(obs))
        , threshold_(threshold)
        , persistence_(persistence)
        , above_threshold_(above)
    {}

    // Evaluar el observable en una configuración
    double evaluate(const Config& config) const {
        return observable_(config);
    }

    // Verificar si una configuración cumple la condición instantánea
    bool checkInstantaneous(const Config& config) const {
        double f = observable_(config);
        return above_threshold_ ? (f > threshold_) : (f < threshold_);
    }

    // Verificar si una trayectoria cumple R (persistencia completa)
    bool checkTrajectory(const std::vector<Config>& trajectory) const {
        if (static_cast<int>(trajectory.size()) < persistence_) {
            return false;
        }

        // Verificar los últimos τ_p pasos
        for (int i = 0; i < persistence_; ++i) {
            size_t idx = trajectory.size() - persistence_ + i;
            if (!checkInstantaneous(trajectory[idx])) {
                return false;
            }
        }
        return true;
    }

    // Actualizar historial y verificar condición R
    bool updateAndCheck(const Config& config) {
        double f = observable_(config);
        history_.push_back(f);

        if (static_cast<int>(history_.size()) > persistence_) {
            history_.pop_front();
        }

        return checkFromHistory();
    }

    // Verificar R desde el historial actual
    bool checkFromHistory() const {
        if (static_cast<int>(history_.size()) < persistence_) {
            return false;
        }

        for (double f : history_) {
            bool meets = above_threshold_ ? (f > threshold_) : (f < threshold_);
            if (!meets) return false;
        }
        return true;
    }

    // Reiniciar historial
    void reset() {
        history_.clear();
    }

    // Estimar P(R) mediante muestreo
    template<typename Sampler>
    double estimateProbability(Sampler& sampler, int num_samples) const {
        int count = 0;
        for (int i = 0; i < num_samples; ++i) {
            auto trajectory = sampler.generateTrajectory(persistence_);
            if (checkTrajectory(trajectory)) {
                ++count;
            }
        }
        return static_cast<double>(count) / num_samples;
    }

    // Getters
    double threshold() const { return threshold_; }
    int persistence() const { return persistence_; }
    bool isAboveThreshold() const { return above_threshold_; }

    // Calcular fracción del tiempo que se cumple R en una trayectoria
    double fractionSatisfied(const std::vector<Config>& trajectory) const {
        if (trajectory.empty()) return 0.0;

        int satisfied = 0;
        for (const auto& config : trajectory) {
            if (checkInstantaneous(config)) {
                ++satisfied;
            }
        }
        return static_cast<double>(satisfied) / trajectory.size();
    }
};

/**
 * Fábrica de condiciones ERS comunes
 */
template<typename Config = std::vector<double>>
class ERSFactory {
public:
    // Condición de energía alta
    static ERSCondition<Config> highEnergy(
        std::function<double(const Config&)> energy,
        double E_c,
        int tau_p)
    {
        return ERSCondition<Config>(energy, E_c, tau_p, true);
    }

    // Condición de energía baja (estado fundamental)
    static ERSCondition<Config> lowEnergy(
        std::function<double(const Config&)> energy,
        double E_c,
        int tau_p)
    {
        return ERSCondition<Config>(energy, E_c, tau_p, false);
    }

    // Condición de norma del orden (|m| > m_c)
    static ERSCondition<Config> orderParameter(
        std::function<double(const Config&)> magnetization,
        double m_c,
        int tau_p)
    {
        auto norm = [magnetization](const Config& c) {
            return std::abs(magnetization(c));
        };
        return ERSCondition<Config>(norm, m_c, tau_p, true);
    }

    // Condición compuesta (AND de múltiples condiciones)
    static ERSCondition<Config> composite(
        const std::vector<ERSCondition<Config>>& conditions)
    {
        auto combined = [conditions](const Config& c) {
            for (const auto& cond : conditions) {
                if (!cond.checkInstantaneous(c)) {
                    return 0.0;  // No cumple alguna
                }
            }
            return 1.0;  // Cumple todas
        };
        int max_persistence = 0;
        for (const auto& cond : conditions) {
            max_persistence = std::max(max_persistence, cond.persistence());
        }
        return ERSCondition<Config>(combined, 0.5, max_persistence, true);
    }
};

} // namespace UE
