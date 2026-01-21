#pragma once
#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <memory>
#include <iostream>

namespace UE {

/**
 * SystemConfig - Configuración de Sistema Dinámico
 *
 * Define un sistema completo para simulación UE:
 * - Potencial/Energía
 * - Dinámica (ecuaciones de movimiento)
 * - Condición R (ERS)
 * - Espacio de configuraciones
 *
 * Permite definir sistemas arbitrarios sin modificar el código.
 */
template<typename ConfigType = std::vector<double>>
struct SystemConfig {
    using Config = ConfigType;
    using PotentialFunc = std::function<double(const Config&)>;
    using DynamicsFunc = std::function<Config(const Config&)>;
    using ObservableFunc = std::function<double(const Config&)>;

    // Identificación
    std::string name;
    std::string description;

    // Dimensión del espacio
    int dimension;

    // Límites del espacio de configuraciones
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;

    // Funciones del sistema
    PotentialFunc potential;      // V(x) - Potencial/Energía
    DynamicsFunc dynamics;        // x_{t+1} = f(x_t)
    ObservableFunc observable;    // θ(x) para condición R

    // Condición R (ERS)
    double theta_threshold;       // θ*
    int tau_persistence;          // τ*
    bool above_threshold;         // true: θ > θ*, false: θ < θ*

    // Parámetros de simulación sugeridos
    double suggested_noise;       // σ sugerido
    int suggested_trajectories;   // N sugerido
    int suggested_length;         // T sugerido

    // Interpretación de resultados
    std::string class_interpretation;  // Qué significan las clases
};

/**
 * SystemFactory - Fábrica de Sistemas Predefinidos
 */
template<typename Config = std::vector<double>>
class SystemFactory {
public:

    /**
     * Doble Pozo Potencial
     * V(x) = (x² - 1)²
     * Dos mínimos en x = ±1, barrera en x = 0
     */
    static SystemConfig<Config> doubleWell() {
        SystemConfig<Config> sys;
        sys.name = "double_well";
        sys.description = "Potencial de doble pozo: V(x) = (x² - 1)²";
        sys.dimension = 1;
        sys.lower_bounds = {-3.0};
        sys.upper_bounds = {3.0};

        sys.potential = [](const Config& x) {
            double val = x[0] * x[0] - 1.0;
            return val * val;
        };

        sys.dynamics = [](const Config& x) {
            const double dt = 0.05;
            double x_val = x[0];
            double force = -4.0 * x_val * x_val * x_val + 4.0 * x_val;
            return Config{x_val + dt * force};
        };

        sys.observable = sys.potential;  // Usar energía como observable
        sys.theta_threshold = 0.5;
        sys.tau_persistence = 10;
        sys.above_threshold = false;  // E < θ*

        sys.suggested_noise = 0.25;
        sys.suggested_trajectories = 3000;
        sys.suggested_length = 800;

        sys.class_interpretation = "Dos pozos (x = ±1) separados por barrera en x = 0";
        return sys;
    }

    /**
     * Triple Pozo Potencial
     * V(x) = x⁶ - 3x⁴ + 2x²
     * Tres mínimos
     */
    static SystemConfig<Config> tripleWell() {
        SystemConfig<Config> sys;
        sys.name = "triple_well";
        sys.description = "Potencial de triple pozo: V(x) = x⁶ - 3x⁴ + 2x²";
        sys.dimension = 1;
        sys.lower_bounds = {-2.0};
        sys.upper_bounds = {2.0};

        sys.potential = [](const Config& x) {
            double x2 = x[0] * x[0];
            double x4 = x2 * x2;
            double x6 = x4 * x2;
            return x6 - 3.0 * x4 + 2.0 * x2;
        };

        sys.dynamics = [](const Config& x) {
            const double dt = 0.02;
            double x_val = x[0];
            double x2 = x_val * x_val;
            double x3 = x2 * x_val;
            double x5 = x3 * x2;
            // F = -dV/dx = -6x⁵ + 12x³ - 4x
            double force = -6.0 * x5 + 12.0 * x3 - 4.0 * x_val;
            return Config{x_val + dt * force};
        };

        sys.observable = sys.potential;
        sys.theta_threshold = 0.3;
        sys.tau_persistence = 15;
        sys.above_threshold = false;

        sys.suggested_noise = 0.2;
        sys.suggested_trajectories = 5000;
        sys.suggested_length = 1000;

        sys.class_interpretation = "Tres pozos: central (x=0) y laterales (x≈±1)";
        return sys;
    }

    /**
     * Pozo Asimétrico
     * V(x) = (x² - 1)² + αx
     * Un pozo más profundo que el otro
     */
    static SystemConfig<Config> asymmetricWell(double alpha = 0.3) {
        SystemConfig<Config> sys;
        sys.name = "asymmetric_well";
        sys.description = "Pozo asimetrico: V(x) = (x² - 1)² + " + std::to_string(alpha) + "x";
        sys.dimension = 1;
        sys.lower_bounds = {-3.0};
        sys.upper_bounds = {3.0};

        sys.potential = [alpha](const Config& x) {
            double val = x[0] * x[0] - 1.0;
            return val * val + alpha * x[0];
        };

        sys.dynamics = [alpha](const Config& x) {
            const double dt = 0.05;
            double x_val = x[0];
            // F = -dV/dx = -4x(x² - 1) - α
            double force = -4.0 * x_val * (x_val * x_val - 1.0) - alpha;
            return Config{x_val + dt * force};
        };

        sys.observable = sys.potential;
        sys.theta_threshold = 0.5;
        sys.tau_persistence = 10;
        sys.above_threshold = false;

        sys.suggested_noise = 0.25;
        sys.suggested_trajectories = 3000;
        sys.suggested_length = 800;

        sys.class_interpretation = "Dos pozos con diferente profundidad (dominancia esperada)";
        return sys;
    }

    /**
     * Oscilador Armónico con Barrera
     * V(x,y) = (x² + y² - 1)²
     * Anillo de mínimos en el círculo unitario
     */
    static SystemConfig<Config> mexicanHat() {
        SystemConfig<Config> sys;
        sys.name = "mexican_hat";
        sys.description = "Sombrero mexicano 2D: V(x,y) = (x² + y² - 1)²";
        sys.dimension = 2;
        sys.lower_bounds = {-2.0, -2.0};
        sys.upper_bounds = {2.0, 2.0};

        sys.potential = [](const Config& x) {
            double r2 = x[0] * x[0] + x[1] * x[1];
            double val = r2 - 1.0;
            return val * val;
        };

        sys.dynamics = [](const Config& x) {
            const double dt = 0.05;
            double x_val = x[0], y_val = x[1];
            double r2 = x_val * x_val + y_val * y_val;
            double factor = -4.0 * (r2 - 1.0);
            double fx = factor * x_val;
            double fy = factor * y_val;
            return Config{x_val + dt * fx, y_val + dt * fy};
        };

        sys.observable = sys.potential;
        sys.theta_threshold = 0.3;
        sys.tau_persistence = 10;
        sys.above_threshold = false;

        sys.suggested_noise = 0.2;
        sys.suggested_trajectories = 5000;
        sys.suggested_length = 500;

        sys.class_interpretation = "Anillo de minimos en r=1 (ruptura continua de simetria)";
        return sys;
    }

    /**
     * Modelo de Ising 1D simplificado
     * Cadena de N spins con interacción de vecinos
     */
    static SystemConfig<Config> ising1D(int num_spins = 10) {
        SystemConfig<Config> sys;
        sys.name = "ising_1d";
        sys.description = "Cadena de Ising 1D con " + std::to_string(num_spins) + " spins";
        sys.dimension = num_spins;
        sys.lower_bounds = std::vector<double>(num_spins, -1.0);
        sys.upper_bounds = std::vector<double>(num_spins, 1.0);

        // Energía: E = -J Σ s_i s_{i+1}
        sys.potential = [num_spins](const Config& s) {
            double E = 0.0;
            for (int i = 0; i < num_spins - 1; ++i) {
                E -= s[i] * s[i + 1];
            }
            return E;
        };

        // Dinámica: Glauber simplificado (relajación hacia mínimo local)
        sys.dynamics = [num_spins](const Config& s) {
            Config next = s;
            const double dt = 0.1;
            for (int i = 0; i < num_spins; ++i) {
                // Campo local: h_i = s_{i-1} + s_{i+1}
                double h = 0.0;
                if (i > 0) h += s[i - 1];
                if (i < num_spins - 1) h += s[i + 1];
                // Relajar hacia sign(h)
                next[i] += dt * (h - s[i]);
                // Clamp a [-1, 1]
                next[i] = std::max(-1.0, std::min(1.0, next[i]));
            }
            return next;
        };

        // Observable: magnetización total
        sys.observable = [num_spins](const Config& s) {
            double m = 0.0;
            for (int i = 0; i < num_spins; ++i) {
                m += s[i];
            }
            return std::abs(m) / num_spins;  // |m| normalizada
        };

        sys.theta_threshold = 0.7;  // Alta magnetización
        sys.tau_persistence = 20;
        sys.above_threshold = true;  // |m| > θ*

        sys.suggested_noise = 0.15;
        sys.suggested_trajectories = 3000;
        sys.suggested_length = 500;

        sys.class_interpretation = "Dos fases ordenadas: todos +1 o todos -1";
        return sys;
    }

    /**
     * Sistema personalizado
     */
    static SystemConfig<Config> custom(
        const std::string& name,
        int dim,
        typename SystemConfig<Config>::PotentialFunc potential,
        typename SystemConfig<Config>::DynamicsFunc dynamics,
        typename SystemConfig<Config>::ObservableFunc observable,
        double theta, int tau, bool above,
        const std::vector<double>& lower,
        const std::vector<double>& upper)
    {
        SystemConfig<Config> sys;
        sys.name = name;
        sys.description = "Sistema personalizado: " + name;
        sys.dimension = dim;
        sys.lower_bounds = lower;
        sys.upper_bounds = upper;
        sys.potential = potential;
        sys.dynamics = dynamics;
        sys.observable = observable;
        sys.theta_threshold = theta;
        sys.tau_persistence = tau;
        sys.above_threshold = above;
        sys.suggested_noise = 0.2;
        sys.suggested_trajectories = 2000;
        sys.suggested_length = 500;
        sys.class_interpretation = "Clases definidas por el usuario";
        return sys;
    }

    /**
     * Listar sistemas disponibles
     */
    static void listSystems() {
        std::cout << "Sistemas disponibles:" << std::endl;
        std::cout << "  double_well    - Doble pozo: V(x) = (x² - 1)²" << std::endl;
        std::cout << "  triple_well    - Triple pozo: V(x) = x⁶ - 3x⁴ + 2x²" << std::endl;
        std::cout << "  asymmetric     - Pozo asimetrico: V(x) = (x² - 1)² + αx" << std::endl;
        std::cout << "  mexican_hat    - Sombrero mexicano 2D: V(r) = (r² - 1)²" << std::endl;
        std::cout << "  ising_1d       - Cadena de Ising 1D" << std::endl;
    }

    /**
     * Obtener sistema por nombre
     */
    static SystemConfig<Config> getSystem(const std::string& name) {
        if (name == "double_well") return doubleWell();
        if (name == "triple_well") return tripleWell();
        if (name == "asymmetric") return asymmetricWell();
        if (name == "mexican_hat") return mexicanHat();
        if (name == "ising_1d") return ising1D();

        // Default
        std::cerr << "Sistema '" << name << "' no encontrado. Usando double_well." << std::endl;
        return doubleWell();
    }
};

} // namespace UE
