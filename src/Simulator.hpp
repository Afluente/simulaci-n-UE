#pragma once
#include <vector>
#include <functional>
#include <random>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>

#include "core/StateSpace.hpp"
#include "core/ERSCondition.hpp"
#include "core/ConditionedProbability.hpp"
#include "core/CostLandscape.hpp"
#include "analysis/MetastableClassDetector.hpp"
#include "analysis/DynamicsAnalyzer.hpp"

namespace UE {

/**
 * Simulator - Simulador Principal del Framework UE
 *
 * Orquesta la simulación completa:
 * 1. Define espacio de estados y condición R
 * 2. Genera trayectorias (condicionadas y no condicionadas)
 * 3. Calcula probabilidades condicionadas p_n(ψ|R)
 * 4. Construye paisaje F_n^R
 * 5. Detecta y valida clases metaestables
 * 6. Analiza dinámica y transiciones
 *
 * Ref: UE Framework - Protocolo Completo (Sección 7)
 */
template<typename ConfigType = std::vector<double>>
class Simulator {
public:
    using Config = ConfigType;
    using Trajectory = std::vector<Config>;
    using Dynamics = std::function<Config(const Config&)>;

    struct SimulationConfig {
        int num_trajectories = 1000;      // Trayectorias a generar
        int trajectory_length = 500;       // Longitud de cada trayectoria
        int landscape_resolution = 100;    // Puntos para el paisaje
        double noise_sigma = 0.1;          // Ruido en la dinámica
        unsigned seed = 42;                // Semilla aleatoria
    };

    struct SimulationResults {
        // Estadísticas básicas
        int total_trajectories;
        int conditioned_trajectories;
        double fraction_satisfying_R;

        // Clases detectadas
        int num_classes;
        int num_validated_classes;

        // Dominancia
        double D_H;
        double Delta;
        int dominant_class_id;

        // Tiempos característicos por clase
        std::vector<double> tau_relax;
        std::vector<double> tau_exit;

        // Matriz de transición
        std::vector<std::vector<double>> transition_matrix;
    };

private:
    SimulationConfig config_;
    StateSpace<Config> state_space_;
    std::unique_ptr<ERSCondition<Config>> ers_condition_;
    std::unique_ptr<ConditionedProbability<Config>> cond_prob_;
    std::unique_ptr<CostLandscape<Config>> landscape_;
    std::unique_ptr<MetastableClassDetector<Config>> class_detector_;
    std::unique_ptr<DynamicsAnalyzer<Config>> dynamics_;

    Dynamics dynamics_func_;
    std::vector<Trajectory> all_trajectories_;

    std::mt19937 rng_;

public:
    Simulator(int dimension, const SimulationConfig& cfg = SimulationConfig())
        : config_(cfg)
        , state_space_(dimension, cfg.seed)
        , rng_(cfg.seed)
    {}

    // Configurar límites del espacio
    void setBounds(const std::vector<double>& lower, const std::vector<double>& upper) {
        state_space_.setBounds(lower, upper);
    }

    // Configurar la condición ERS (R)
    void setERSCondition(
        std::function<double(const Config&)> observable,
        double threshold,
        int persistence,
        bool above = true)
    {
        ers_condition_ = std::make_unique<ERSCondition<Config>>(
            observable, threshold, persistence, above
        );
    }

    // Configurar la dinámica del sistema
    void setDynamics(Dynamics func) {
        dynamics_func_ = std::move(func);
    }

    // Configurar nivel de coarse-graining
    void addCoarseGrainLevel(int n, double delta_n, double epsilon_n) {
        state_space_.addCoarseGrainLevel(n, delta_n, epsilon_n);
    }

    // Generar una trayectoria
    Trajectory generateTrajectory(const Config& initial) {
        Trajectory traj;
        traj.reserve(config_.trajectory_length);

        Config current = initial;
        traj.push_back(current);

        std::normal_distribution<double> noise(0.0, config_.noise_sigma);

        for (int t = 1; t < config_.trajectory_length; ++t) {
            // Aplicar dinámica
            current = dynamics_func_(current);

            // Añadir ruido
            for (size_t d = 0; d < current.size(); ++d) {
                current[d] += noise(rng_);
                // Aplicar límites
                current[d] = std::max(state_space_.lowerBounds()[d],
                             std::min(state_space_.upperBounds()[d], current[d]));
            }

            traj.push_back(current);
        }

        return traj;
    }

    // Ejecutar simulación completa
    SimulationResults run() {
        SimulationResults results{};

        std::cout << "=== SIMULACION UNIVERSO EMERGENTE ===" << std::endl;
        std::cout << std::endl;

        // Paso 1: Generar trayectorias
        std::cout << "[1] Generando " << config_.num_trajectories << " trayectorias..." << std::endl;
        all_trajectories_.clear();
        all_trajectories_.reserve(config_.num_trajectories);

        for (int i = 0; i < config_.num_trajectories; ++i) {
            Config initial = state_space_.randomConfig();
            auto traj = generateTrajectory(initial);
            all_trajectories_.push_back(std::move(traj));
        }
        results.total_trajectories = config_.num_trajectories;

        // Paso 2: Filtrar trayectorias que cumplen R
        std::cout << "[2] Filtrando trayectorias bajo condicion R..." << std::endl;
        cond_prob_ = std::make_unique<ConditionedProbability<Config>>(
            state_space_, *ers_condition_
        );

        int count_R = 0;
        for (const auto& traj : all_trajectories_) {
            if (ers_condition_->checkTrajectory(traj)) {
                cond_prob_->addConditionedTrajectory(traj);
                ++count_R;
            }
        }
        results.conditioned_trajectories = count_R;
        results.fraction_satisfying_R = static_cast<double>(count_R) / config_.num_trajectories;

        std::cout << "    Trayectorias bajo R: " << count_R
                  << " (" << std::fixed << std::setprecision(1)
                  << (results.fraction_satisfying_R * 100) << "%)" << std::endl;

        if (count_R == 0) {
            std::cout << "    [!] No hay trayectorias que cumplan R. Ajustar parametros." << std::endl;
            return results;
        }

        // Paso 3: Construir paisaje F_n^R
        std::cout << "[3] Construyendo paisaje condicionado F_n^R..." << std::endl;
        landscape_ = std::make_unique<CostLandscape<Config>>(
            state_space_, *cond_prob_
        );

        if (state_space_.dimension() == 1) {
            landscape_->buildGrid1D(config_.landscape_resolution);
        } else if (state_space_.dimension() == 2) {
            landscape_->buildGrid2D(static_cast<int>(std::sqrt(config_.landscape_resolution)));
        }

        auto minima = landscape_->getMinima();
        std::cout << "    Minimos locales encontrados: " << minima.size() << std::endl;

        // Paso 4: Detectar clases metaestables
        std::cout << "[4] Detectando clases metaestables..." << std::endl;
        class_detector_ = std::make_unique<MetastableClassDetector<Config>>(
            state_space_, *landscape_
        );
        class_detector_->detectFromLandscape();

        results.num_classes = static_cast<int>(class_detector_->numClasses());
        std::cout << "    Clases candidatas: " << results.num_classes << std::endl;

        // Paso 5: Validar dinámicamente
        std::cout << "[5] Validando clases con dinamica..." << std::endl;
        dynamics_ = std::make_unique<DynamicsAnalyzer<Config>>(
            state_space_, *class_detector_
        );

        // Convertir trayectorias condicionadas a vector
        std::vector<Trajectory> cond_trajs;
        for (const auto& traj : all_trajectories_) {
            if (ers_condition_->checkTrajectory(traj)) {
                cond_trajs.push_back(traj);
            }
        }

        dynamics_->validateAllClasses(cond_trajs);

        auto validated = class_detector_->getValidatedClasses();
        results.num_validated_classes = static_cast<int>(validated.size());
        std::cout << "    Clases validadas: " << results.num_validated_classes << std::endl;

        // Recopilar tiempos
        for (const auto& c : class_detector_->getClasses()) {
            results.tau_relax.push_back(c.tau_relax);
            results.tau_exit.push_back(c.tau_exit_mean);
        }

        // Paso 6: Calcular dominancia
        std::cout << "[6] Calculando metricas de dominancia..." << std::endl;
        auto dom = class_detector_->calculateDominance();
        results.D_H = dom.D_H;
        results.Delta = dom.Delta;
        results.dominant_class_id = dom.dominant_class_id;

        std::cout << "    D_H (p1/p2) = " << std::fixed << std::setprecision(2)
                  << results.D_H << std::endl;
        std::cout << "    Delta (p1-p2) = " << std::fixed << std::setprecision(4)
                  << results.Delta << std::endl;

        // Paso 7: Matriz de transición
        std::cout << "[7] Calculando matriz de transicion..." << std::endl;
        auto trans = dynamics_->calculateTransitionMatrix(cond_trajs);
        results.transition_matrix = trans.P;

        // Imprimir resumen
        printSummary(results);

        return results;
    }

    // Imprimir resumen de resultados
    void printSummary(const SimulationResults& results) const {
        std::cout << std::endl;
        std::cout << "=== RESUMEN DE RESULTADOS ===" << std::endl;
        std::cout << std::endl;

        std::cout << "Trayectorias totales: " << results.total_trajectories << std::endl;
        std::cout << "Trayectorias bajo R: " << results.conditioned_trajectories
                  << " (" << std::fixed << std::setprecision(1)
                  << (results.fraction_satisfying_R * 100) << "%)" << std::endl;
        std::cout << std::endl;

        std::cout << "Clases metaestables detectadas: " << results.num_classes << std::endl;
        std::cout << "Clases validadas (τ_exit >> τ_relax): " << results.num_validated_classes << std::endl;
        std::cout << std::endl;

        if (!results.tau_relax.empty()) {
            std::cout << "Tiempos caracteristicos por clase:" << std::endl;
            for (size_t k = 0; k < results.tau_relax.size(); ++k) {
                std::cout << "  Clase " << k << ": τ_relax=" << std::fixed << std::setprecision(1)
                          << results.tau_relax[k] << ", τ_exit=" << results.tau_exit[k] << std::endl;
            }
            std::cout << std::endl;
        }

        std::cout << "Dominancia:" << std::endl;
        std::cout << "  Clase dominante: " << results.dominant_class_id << std::endl;
        std::cout << "  D_H = " << std::fixed << std::setprecision(2) << results.D_H << std::endl;
        std::cout << "  Δ = " << std::fixed << std::setprecision(4) << results.Delta << std::endl;
        std::cout << std::endl;

        if (!results.transition_matrix.empty()) {
            std::cout << "Matriz de transicion P[i→j]:" << std::endl;
            std::cout << "     ";
            for (size_t j = 0; j < results.transition_matrix.size(); ++j) {
                std::cout << std::setw(8) << j;
            }
            std::cout << std::endl;
            for (size_t i = 0; i < results.transition_matrix.size(); ++i) {
                std::cout << "  " << i << ": ";
                for (size_t j = 0; j < results.transition_matrix[i].size(); ++j) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(3)
                              << results.transition_matrix[i][j];
                }
                std::cout << std::endl;
            }
        }
    }

    // Acceso a componentes
    const StateSpace<Config>& stateSpace() const { return state_space_; }
    const std::vector<Trajectory>& trajectories() const { return all_trajectories_; }
    const CostLandscape<Config>* landscape() const { return landscape_.get(); }
    const MetastableClassDetector<Config>* classDetector() const { return class_detector_.get(); }
    const ConditionedProbability<Config>* conditionedProb() const { return cond_prob_.get(); }
    const ERSCondition<Config>* ersCondition() const { return ers_condition_.get(); }

    // Datos para visualización

    // Obtener datos de probabilidad p_n(ψ|R) en un grid
    std::vector<std::pair<Config, double>> getProbabilityData(int num_points = 100) const {
        std::vector<std::pair<Config, double>> data;
        if (!cond_prob_ || state_space_.dimension() != 1) return data;

        double lower = state_space_.lowerBounds()[0];
        double upper = state_space_.upperBounds()[0];
        double step = (upper - lower) / (num_points - 1);

        for (int i = 0; i < num_points; ++i) {
            Config psi(1);
            psi[0] = lower + i * step;
            auto result = cond_prob_->calculate(psi);
            data.emplace_back(psi, result.probability);
        }
        return data;
    }

    // Obtener datos del paisaje F_n^R
    std::vector<std::pair<Config, double>> getCostData() const {
        std::vector<std::pair<Config, double>> data;
        if (!landscape_) return data;

        for (const auto& point : landscape_->getLandscape()) {
            data.emplace_back(point.config, point.cost);
        }
        return data;
    }

    // Obtener centros de clases
    std::vector<Config> getClassCenters() const {
        std::vector<Config> centers;
        if (!class_detector_) return centers;

        for (const auto& c : class_detector_->getClasses()) {
            centers.push_back(c.center);
        }
        return centers;
    }

    // Obtener trayectorias condicionadas
    std::vector<Trajectory> getConditionedTrajectories() const {
        std::vector<Trajectory> cond_trajs;
        if (!ers_condition_) return cond_trajs;

        for (const auto& traj : all_trajectories_) {
            if (ers_condition_->checkTrajectory(traj)) {
                cond_trajs.push_back(traj);
            }
        }
        return cond_trajs;
    }

    // Obtener valores x de configuraciones bajo R (para histograma)
    std::vector<double> getConditionedXValues() const {
        std::vector<double> values;
        if (!cond_prob_) return values;

        for (const auto& sample : cond_prob_->getSamples()) {
            if (!sample.empty()) {
                values.push_back(sample[0]);
            }
        }
        return values;
    }
};

} // namespace UE
