#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "../core/StateSpace.hpp"
#include "MetastableClassDetector.hpp"

namespace UE {

/**
 * DynamicsAnalyzer - Analizador de Dinámica bajo R
 *
 * Mide propiedades dinámicas clave:
 * - τ_relax(A_k): Tiempo de relajación dentro de la clase
 * - τ_exit(A_k): Tiempo de primera salida de la clase
 * - P(A_k → A_j): Matriz de transición entre clases
 * - Histéresis y dependencia de condiciones iniciales
 *
 * Ref: UE Framework Sección 5.2 - Dinámica y Transiciones
 */
template<typename ConfigType = std::vector<double>>
class DynamicsAnalyzer {
public:
    using Config = ConfigType;
    using Trajectory = std::vector<Config>;
    using Class = MetastableClass<Config>;

    struct RelaxationResult {
        double tau_relax;         // Tiempo característico
        double autocorr_decay;    // Decaimiento de autocorrelación
        std::vector<double> autocorrelation;
    };

    struct ExitTimeResult {
        double mean;              // E[τ_exit]
        double std_dev;           // Desviación estándar
        double median;            // Mediana
        std::vector<int> exit_times;  // Tiempos individuales
        int num_samples;
    };

    struct TransitionMatrix {
        std::vector<std::vector<double>> P;  // P[i][j] = prob de ir de i a j
        std::vector<int> class_ids;          // IDs de las clases
    };

private:
    const StateSpace<Config>* state_space_;
    MetastableClassDetector<Config>* class_detector_;

public:
    DynamicsAnalyzer(
        const StateSpace<Config>& space,
        MetastableClassDetector<Config>& detector)
        : state_space_(&space)
        , class_detector_(&detector)
    {}

    // Calcular tiempo de relajación desde autocorrelación
    RelaxationResult calculateRelaxation(
        const Trajectory& trajectory,
        int class_id,
        int max_lag = 100) const
    {
        RelaxationResult result;
        result.tau_relax = 0;
        result.autocorr_decay = 0;

        // Filtrar solo configuraciones dentro de la clase
        std::vector<Config> in_class;
        for (const auto& config : trajectory) {
            if (class_detector_->classifyConfig(config) == class_id) {
                in_class.push_back(config);
            }
        }

        if (in_class.size() < static_cast<size_t>(max_lag)) {
            return result;
        }

        // Calcular autocorrelación de la primera componente
        // (simplificado para demostración)
        std::vector<double> x(in_class.size());
        for (size_t i = 0; i < in_class.size(); ++i) {
            x[i] = in_class[i][0];
        }

        // Media
        double mean = 0;
        for (double xi : x) mean += xi;
        mean /= x.size();

        // Varianza
        double var = 0;
        for (double xi : x) var += (xi - mean) * (xi - mean);
        var /= x.size();

        if (var < 1e-10) {
            result.tau_relax = 1;
            return result;
        }

        // Autocorrelación
        result.autocorrelation.resize(max_lag);
        for (int lag = 0; lag < max_lag; ++lag) {
            double sum = 0;
            int count = 0;
            for (size_t i = 0; i < x.size() - lag; ++i) {
                sum += (x[i] - mean) * (x[i + lag] - mean);
                ++count;
            }
            result.autocorrelation[lag] = (count > 0) ? sum / (count * var) : 0;
        }

        // Estimar τ_relax como el lag donde autocorr cae a 1/e
        double threshold = 1.0 / std::exp(1.0);
        for (int lag = 0; lag < max_lag; ++lag) {
            if (result.autocorrelation[lag] < threshold) {
                result.tau_relax = lag;
                break;
            }
        }
        if (result.tau_relax == 0) {
            result.tau_relax = max_lag;  // No decayó suficiente
        }

        return result;
    }

    // Calcular tiempos de salida de una clase
    ExitTimeResult calculateExitTimes(
        const std::vector<Trajectory>& trajectories,
        int class_id) const
    {
        ExitTimeResult result;
        result.mean = 0;
        result.std_dev = 0;
        result.median = 0;
        result.num_samples = 0;

        for (const auto& traj : trajectories) {
            // Encontrar primera entrada a la clase
            int entry_time = -1;
            for (size_t t = 0; t < traj.size(); ++t) {
                if (class_detector_->classifyConfig(traj[t]) == class_id) {
                    entry_time = static_cast<int>(t);
                    break;
                }
            }

            if (entry_time < 0) continue;

            // Encontrar primera salida después de entrada
            for (size_t t = entry_time + 1; t < traj.size(); ++t) {
                if (class_detector_->classifyConfig(traj[t]) != class_id) {
                    int exit_time = static_cast<int>(t) - entry_time;
                    result.exit_times.push_back(exit_time);
                    break;
                }
            }
        }

        result.num_samples = static_cast<int>(result.exit_times.size());

        if (result.num_samples > 0) {
            // Media
            double sum = 0;
            for (int t : result.exit_times) sum += t;
            result.mean = sum / result.num_samples;

            // Desviación estándar
            double sq_sum = 0;
            for (int t : result.exit_times) {
                double diff = t - result.mean;
                sq_sum += diff * diff;
            }
            result.std_dev = std::sqrt(sq_sum / result.num_samples);

            // Mediana
            std::vector<int> sorted = result.exit_times;
            std::sort(sorted.begin(), sorted.end());
            result.median = sorted[sorted.size() / 2];
        }

        return result;
    }

    // Calcular matriz de transición entre clases
    TransitionMatrix calculateTransitionMatrix(
        const std::vector<Trajectory>& trajectories) const
    {
        TransitionMatrix result;

        const auto& classes = class_detector_->getClasses();
        int K = static_cast<int>(classes.size());

        result.P.resize(K, std::vector<double>(K, 0.0));
        for (const auto& c : classes) {
            result.class_ids.push_back(c.id);
        }

        // Contar transiciones
        std::vector<std::vector<int>> counts(K, std::vector<int>(K, 0));
        std::vector<int> total_from(K, 0);

        for (const auto& traj : trajectories) {
            int prev_class = -1;
            for (const auto& config : traj) {
                int curr_class = class_detector_->classifyConfig(config);
                if (curr_class >= 0 && curr_class < K) {
                    if (prev_class >= 0 && prev_class < K && prev_class != curr_class) {
                        counts[prev_class][curr_class]++;
                        total_from[prev_class]++;
                    }
                    prev_class = curr_class;
                }
            }
        }

        // Normalizar a probabilidades
        for (int i = 0; i < K; ++i) {
            if (total_from[i] > 0) {
                for (int j = 0; j < K; ++j) {
                    result.P[i][j] = static_cast<double>(counts[i][j]) / total_from[i];
                }
            }
        }

        return result;
    }

    // Validar todas las clases con datos dinámicos
    void validateAllClasses(const std::vector<Trajectory>& trajectories) {
        for (const auto& c : class_detector_->getClasses()) {
            // Calcular τ_relax (usando la primera trayectoria larga)
            double tau_relax = 1.0;
            for (const auto& traj : trajectories) {
                if (traj.size() > 100) {
                    auto relax_result = calculateRelaxation(traj, c.id);
                    if (relax_result.tau_relax > 0) {
                        tau_relax = relax_result.tau_relax;
                        break;
                    }
                }
            }

            // Calcular τ_exit
            auto exit_result = calculateExitTimes(trajectories, c.id);

            // Validar la clase
            class_detector_->validateClass(
                c.id,
                tau_relax,
                exit_result.mean,
                exit_result.std_dev
            );
        }
    }

    // Detectar histéresis: dependencia de condición inicial
    struct HysteresisResult {
        bool detected;
        double asymmetry;  // Diferencia en probabilidades finales
        std::vector<double> final_probs_forward;   // Comenzando desde clase 0
        std::vector<double> final_probs_backward;  // Comenzando desde clase 1
    };

    HysteresisResult detectHysteresis(
        const std::vector<Trajectory>& forward_trajectories,
        const std::vector<Trajectory>& backward_trajectories) const
    {
        HysteresisResult result;
        result.detected = false;
        result.asymmetry = 0;

        const auto& classes = class_detector_->getClasses();
        int K = static_cast<int>(classes.size());

        result.final_probs_forward.resize(K, 0);
        result.final_probs_backward.resize(K, 0);

        // Contar estados finales desde forward
        int count_f = 0;
        for (const auto& traj : forward_trajectories) {
            if (!traj.empty()) {
                int final_class = class_detector_->classifyConfig(traj.back());
                if (final_class >= 0 && final_class < K) {
                    result.final_probs_forward[final_class]++;
                    count_f++;
                }
            }
        }
        if (count_f > 0) {
            for (auto& p : result.final_probs_forward) p /= count_f;
        }

        // Contar estados finales desde backward
        int count_b = 0;
        for (const auto& traj : backward_trajectories) {
            if (!traj.empty()) {
                int final_class = class_detector_->classifyConfig(traj.back());
                if (final_class >= 0 && final_class < K) {
                    result.final_probs_backward[final_class]++;
                    count_b++;
                }
            }
        }
        if (count_b > 0) {
            for (auto& p : result.final_probs_backward) p /= count_b;
        }

        // Calcular asimetría
        double max_diff = 0;
        for (int k = 0; k < K; ++k) {
            double diff = std::abs(result.final_probs_forward[k] -
                                   result.final_probs_backward[k]);
            max_diff = std::max(max_diff, diff);
        }
        result.asymmetry = max_diff;
        result.detected = (max_diff > 0.1);  // Umbral arbitrario

        return result;
    }
};

} // namespace UE
