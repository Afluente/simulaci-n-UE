/**
 * Simulación del Framework Universo Emergente (UE)
 *
 * Ejemplo: Sistema de Doble Pozo Potencial
 *
 * Este ejemplo demuestra los conceptos centrales del framework UE:
 * - Selección Condicionada (SC): qué es típico bajo la condición R
 * - Evento Raro Sostenido (ERS): condición R = energía baja sostenida
 * - Clases Metaestables: los dos pozos como atractores
 * - Dinámica: τ_relax << τ_exit para cada pozo
 * - Transiciones: saltos entre pozos
 *
 * Potencial: V(x) = (x² - 1)² = x⁴ - 2x² + 1
 *   - Mínimos en x = ±1 (dos pozos)
 *   - Máximo (barrera) en x = 0
 *
 * Ref: UE Framework v1.5
 */

#include <iostream>
#include <cmath>
#include <string>
#include "src/Simulator.hpp"
#include "src/Visualization.hpp"

using Config = std::vector<double>;

// Potencial de doble pozo: V(x) = (x² - 1)²
double potential(const Config& x) {
    double val = x[0] * x[0] - 1.0;
    return val * val;
}

// Fuerza: F(x) = -dV/dx = -4x(x² - 1) = -4x³ + 4x
double force(const Config& x) {
    double x_val = x[0];
    return -4.0 * x_val * x_val * x_val + 4.0 * x_val;
}

// Dinámica: descenso por gradiente con paso dt
Config dynamics(const Config& x) {
    const double dt = 0.05;  // Paso de tiempo
    Config next = x;
    next[0] += dt * force(x);
    return next;
}

// Energía total (para la condición R)
double energy(const Config& x) {
    return potential(x);
}

void printUsage(const char* prog) {
    std::cout << "Uso: " << prog << " [opciones]" << std::endl;
    std::cout << std::endl;
    std::cout << "Opciones:" << std::endl;
    std::cout << "  --no-plot    No mostrar graficos (solo texto)" << std::endl;
    std::cout << "  --help       Mostrar esta ayuda" << std::endl;
}

int main(int argc, char* argv[]) {
    bool show_plots = true;

    // Parsear argumentos
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-plot") {
            show_plots = false;
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }

    std::cout << "============================================" << std::endl;
    std::cout << "  SIMULACION UNIVERSO EMERGENTE (UE)" << std::endl;
    std::cout << "  Sistema: Doble Pozo Potencial" << std::endl;
    std::cout << "  V(x) = (x² - 1)²" << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::endl;

    // Configuración de la simulación
    UE::Simulator<Config>::SimulationConfig cfg;
    cfg.num_trajectories = 2000;     // Más trayectorias para mejor estadística
    cfg.trajectory_length = 500;     // Pasos por trayectoria
    cfg.noise_sigma = 0.3;           // Ruido térmico (permite transiciones)
    cfg.seed = 42;
    cfg.landscape_resolution = 100;

    // Crear simulador 1D
    UE::Simulator<Config> sim(1, cfg);

    // Configurar límites del espacio [-3, 3]
    sim.setBounds({-3.0}, {3.0});

    // Configurar nivel de coarse-graining
    // n=1, δ_n=0.3 (resolución), ε_n=0.5 (temperatura)
    sim.addCoarseGrainLevel(1, 0.3, 0.5);

    // Configurar la dinámica del sistema
    sim.setDynamics(dynamics);

    // Configurar la condición ERS (R):
    // R = "energía menor que 0.5 durante al menos 10 pasos"
    // Esto selecciona trayectorias que pasan tiempo en los pozos
    double E_threshold = 0.5;
    int tau_p = 10;  // Persistencia mínima

    std::cout << "Condicion R (ERS): E(x) < " << E_threshold
              << " durante τ_p = " << tau_p << " pasos" << std::endl;
    std::cout << std::endl;

    sim.setERSCondition(energy, E_threshold, tau_p, false);  // false = menor que

    // Ejecutar simulación
    auto results = sim.run();

    // Información adicional sobre el sistema
    std::cout << std::endl;
    std::cout << "=== INTERPRETACION ===" << std::endl;
    std::cout << std::endl;
    std::cout << "El sistema de doble pozo tiene:" << std::endl;
    std::cout << "  - Dos minimos (pozos) en x = +1 y x = -1" << std::endl;
    std::cout << "  - Una barrera en x = 0" << std::endl;
    std::cout << std::endl;

    if (results.num_classes >= 2) {
        std::cout << "Bajo la condicion R (energia baja sostenida):" << std::endl;
        std::cout << "  - Se detectaron " << results.num_classes << " clases metaestables" << std::endl;
        std::cout << "  - Esto corresponde a los dos pozos del potencial" << std::endl;
        std::cout << std::endl;

        if (results.D_H > 1.5) {
            std::cout << "Dominancia detectada (D_H = " << results.D_H << "):" << std::endl;
            std::cout << "  - Una clase es mas probable que la otra bajo R" << std::endl;
            std::cout << "  - Esto puede deberse a asimetria en condiciones iniciales" << std::endl;
        } else {
            std::cout << "Equilibrio aproximado (D_H = " << results.D_H << "):" << std::endl;
            std::cout << "  - Ambas clases tienen probabilidad similar bajo R" << std::endl;
            std::cout << "  - El sistema es aproximadamente simetrico" << std::endl;
        }
    } else if (results.num_classes == 1) {
        std::cout << "Solo se detecto una clase dominante." << std::endl;
        std::cout << "  - El ruido puede ser muy bajo para transiciones" << std::endl;
        std::cout << "  - O la condicion R es muy restrictiva" << std::endl;
    } else {
        std::cout << "No se detectaron clases metaestables." << std::endl;
        std::cout << "  - Verificar parametros de la condicion R" << std::endl;
    }

    // === VISUALIZACIÓN ===
    if (show_plots && results.conditioned_trajectories > 0) {
        std::cout << std::endl;
        std::cout << "=== GENERANDO GRAFICOS ===" << std::endl;
        std::cout << std::endl;

        // Crear visualizador
        UE::Visualization viz(".", true);  // Directorio actual, interactivo

        // Obtener datos para visualización
        auto prob_data = sim.getProbabilityData(100);
        auto cost_data = sim.getCostData();
        auto class_centers = sim.getClassCenters();
        auto cond_trajs = sim.getConditionedTrajectories();
        auto x_values = sim.getConditionedXValues();

        // Verificar si gnuplot está disponible
        if (UE::Visualization::isGnuplotAvailable()) {
            std::cout << "Usando gnuplot para graficos..." << std::endl;
            std::cout << "(Cierre la ventana del grafico para continuar)" << std::endl;

            // Mostrar gráfico combinado con todos los resultados
            viz.plotCombined(
                x_values,
                prob_data,
                cost_data,
                class_centers,
                cond_trajs,
                "Simulacion UE: Doble Pozo Potencial"
            );

            std::cout << std::endl;
            std::cout << "Graficos adicionales disponibles:" << std::endl;
            std::cout << "  - Histograma de configuraciones" << std::endl;
            std::cout << "  - Trayectorias individuales" << std::endl;
            std::cout << "  - Comparacion potencial vs paisaje" << std::endl;

            // Mostrar histograma
            std::cout << std::endl;
            std::cout << "Mostrando histograma de configuraciones bajo R..." << std::endl;
            viz.plotHistogram(x_values, 50, "Distribucion de Configuraciones bajo R", "x");

            // Mostrar trayectorias
            std::cout << std::endl;
            std::cout << "Mostrando trayectorias de muestra..." << std::endl;
            viz.plotTrajectories(cond_trajs, 8, "Trayectorias bajo Condicion R");

            // Comparación potencial vs paisaje condicionado
            std::cout << std::endl;
            std::cout << "Mostrando comparacion V(x) vs F_n^R..." << std::endl;
            viz.plotPotentialComparison(
                potential,
                cost_data,
                -3.0, 3.0, 100,
                "Potencial Original vs Paisaje Condicionado"
            );
        } else {
            std::cout << "gnuplot no encontrado. Usando visualizacion ASCII..." << std::endl;
            std::cout << "(Instala gnuplot para graficos interactivos: sudo dnf install gnuplot)" << std::endl;

            // Usar visualización ASCII
            viz.plotCombinedASCII(x_values, prob_data, cost_data, class_centers);
        }

        // Guardar datos para uso externo
        std::cout << std::endl;
        std::cout << "Datos guardados en archivos .dat para uso externo:" << std::endl;
        std::cout << "  - hist_combined.dat   (histograma)" << std::endl;
        std::cout << "  - prob_combined.dat   (probabilidad)" << std::endl;
        std::cout << "  - cost_combined.dat   (paisaje costo)" << std::endl;
        std::cout << "  - pot_combined.dat    (potencial original)" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "=== FIN DE SIMULACION ===" << std::endl;

    return 0;
}
