/**
 * Simulación del Framework Universo Emergente (UE)
 *
 * Simulador configurable para diferentes sistemas dinámicos.
 * Implementa el protocolo completo del marco UE v1.5.
 *
 * Uso:
 *   ./simulacion_ue [opciones]
 *
 * Sistemas disponibles:
 *   double_well  - Potencial de doble pozo (default)
 *   triple_well  - Potencial de triple pozo
 *   asymmetric   - Pozo asimétrico (con dominancia)
 *   mexican_hat  - Sombrero mexicano 2D
 *   ising_1d     - Cadena de Ising 1D
 *
 * Ref: UE Framework v1.5
 */

#include <iostream>
#include <cmath>
#include <string>
#include "src/Simulator.hpp"
#include "src/Visualization.hpp"
#include "src/analysis/ValidationSuite.hpp"
#include "src/systems/SystemConfig.hpp"

using Config = std::vector<double>;

void printUsage(const char* prog) {
    std::cout << "Uso: " << prog << " [opciones]" << std::endl;
    std::cout << std::endl;
    std::cout << "Opciones:" << std::endl;
    std::cout << "  --system <nombre>  Sistema a simular (default: double_well)" << std::endl;
    std::cout << "  --list             Listar sistemas disponibles" << std::endl;
    std::cout << "  --no-plot          No mostrar graficos (solo texto)" << std::endl;
    std::cout << "  --validate         Ejecutar protocolo de validacion UE" << std::endl;
    std::cout << "  --quick            Simulacion rapida (menos trayectorias)" << std::endl;
    std::cout << "  --seed <n>         Semilla aleatoria (default: 42)" << std::endl;
    std::cout << "  --noise <σ>        Nivel de ruido (override)" << std::endl;
    std::cout << "  --traj <n>         Numero de trayectorias (override)" << std::endl;
    std::cout << "  --help             Mostrar esta ayuda" << std::endl;
    std::cout << std::endl;
    std::cout << "Ejemplos:" << std::endl;
    std::cout << "  " << prog << " --system double_well --validate" << std::endl;
    std::cout << "  " << prog << " --system asymmetric --no-plot" << std::endl;
    std::cout << "  " << prog << " --system ising_1d --quick" << std::endl;
}

template<typename ConfigType>
void runSimulation(
    const UE::SystemConfig<ConfigType>& system,
    bool show_plots,
    bool run_validation,
    bool quick_mode,
    unsigned seed,
    double noise_override,
    int traj_override)
{
    std::cout << "============================================" << std::endl;
    std::cout << "  SIMULACION UNIVERSO EMERGENTE (UE)" << std::endl;
    std::cout << "  Sistema: " << system.name << std::endl;
    std::cout << "  " << system.description << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << std::endl;

    // Configuración de la simulación
    typename UE::Simulator<ConfigType>::SimulationConfig cfg;
    if (quick_mode) {
        cfg.num_trajectories = 500;
        cfg.trajectory_length = 200;
    } else {
        cfg.num_trajectories = (traj_override > 0) ? traj_override : system.suggested_trajectories;
        cfg.trajectory_length = system.suggested_length;
    }
    cfg.noise_sigma = (noise_override > 0) ? noise_override : system.suggested_noise;
    cfg.seed = seed;
    cfg.landscape_resolution = 100;

    // Crear simulador
    UE::Simulator<ConfigType> sim(system.dimension, cfg);

    // Configurar límites del espacio
    sim.setBounds(system.lower_bounds, system.upper_bounds);

    // Configurar nivel de coarse-graining
    sim.addCoarseGrainLevel(1, 0.3, 0.5);

    // Configurar la dinámica del sistema
    sim.setDynamics(system.dynamics);

    // Configurar la condición ERS (R)
    std::cout << "Condicion R (ERS): θ(x) " << (system.above_threshold ? ">" : "<")
              << " " << system.theta_threshold
              << " durante τ_p = " << system.tau_persistence << " pasos" << std::endl;
    std::cout << std::endl;

    sim.setERSCondition(
        system.observable,
        system.theta_threshold,
        system.tau_persistence,
        system.above_threshold
    );

    // Ejecutar simulación
    auto results = sim.run();

    // Información adicional sobre el sistema
    std::cout << std::endl;
    std::cout << "=== INTERPRETACION ===" << std::endl;
    std::cout << std::endl;
    std::cout << system.class_interpretation << std::endl;
    std::cout << std::endl;

    if (results.num_classes >= 2) {
        std::cout << "Bajo la condicion R:" << std::endl;
        std::cout << "  - Se detectaron " << results.num_classes << " clases metaestables" << std::endl;
        std::cout << std::endl;

        if (results.D_H > 1.5) {
            std::cout << "Dominancia detectada (D_H = " << results.D_H << "):" << std::endl;
            std::cout << "  - Una clase es mas probable que la otra bajo R" << std::endl;
        } else {
            std::cout << "Equilibrio aproximado (D_H = " << results.D_H << "):" << std::endl;
            std::cout << "  - Las clases tienen probabilidad similar bajo R" << std::endl;
        }
    } else if (results.num_classes == 1) {
        std::cout << "Solo se detecto una clase dominante." << std::endl;
    } else {
        std::cout << "No se detectaron clases metaestables." << std::endl;
        std::cout << "  - Verificar parametros de la condicion R" << std::endl;
    }

    // === VISUALIZACIÓN ===
    if (show_plots && results.conditioned_trajectories > 0 && system.dimension == 1) {
        std::cout << std::endl;
        std::cout << "=== GENERANDO GRAFICOS ===" << std::endl;
        std::cout << std::endl;

        UE::Visualization viz(".", true);

        auto prob_data = sim.getProbabilityData(100);
        auto cost_data = sim.getCostData();
        auto class_centers = sim.getClassCenters();
        auto cond_trajs = sim.getConditionedTrajectories();
        auto x_values = sim.getConditionedXValues();

        if (UE::Visualization::isGnuplotAvailable()) {
            std::cout << "Usando gnuplot para graficos..." << std::endl;
            std::cout << "(Cierre la ventana del grafico para continuar)" << std::endl;

            viz.plotCombined(
                x_values,
                prob_data,
                cost_data,
                class_centers,
                cond_trajs,
                "Simulacion UE: " + system.name
            );

            std::cout << std::endl;
            std::cout << "Mostrando histograma..." << std::endl;
            viz.plotHistogram(x_values, 50, "Distribucion bajo R", "x");

            std::cout << "Mostrando trayectorias..." << std::endl;
            viz.plotTrajectories(cond_trajs, 8, "Trayectorias bajo R");

            std::cout << "Mostrando comparacion potencial vs paisaje..." << std::endl;
            viz.plotPotentialComparison(
                system.potential,
                cost_data,
                system.lower_bounds[0], system.upper_bounds[0], 100,
                "V(x) vs F_n^R"
            );
        } else {
            std::cout << "gnuplot no encontrado. Usando ASCII..." << std::endl;
            viz.plotCombinedASCII(x_values, prob_data, cost_data, class_centers);
        }

        std::cout << std::endl;
        std::cout << "Datos guardados en archivos .dat" << std::endl;
    }

    // === VALIDACIÓN UE ===
    if (run_validation && results.conditioned_trajectories > 0) {
        std::cout << std::endl;

        UE::ValidationSuite<ConfigType> validation(
            system.observable,
            system.dynamics,
            system.theta_threshold,
            system.tau_persistence,
            system.above_threshold,
            seed
        );

        auto all_trajs = sim.trajectories();

        auto val_results = validation.runFullValidation(
            sim,
            all_trajs,
            results.num_classes,
            results.D_H
        );

        if (val_results.all_passed) {
            std::cout << std::endl;
            std::cout << "La simulacion ha pasado todos los tests de validacion." << std::endl;
        } else {
            std::cout << std::endl;
            std::cout << "ATENCION: Algunos tests no pasaron. Revisar parametros." << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "=== FIN DE SIMULACION ===" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string system_name = "double_well";
    bool show_plots = true;
    bool run_validation = false;
    bool quick_mode = false;
    unsigned seed = 42;
    double noise_override = -1.0;
    int traj_override = -1;

    // Parsear argumentos
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--system" && i + 1 < argc) {
            system_name = argv[++i];
        } else if (arg == "--list") {
            UE::SystemFactory<Config>::listSystems();
            return 0;
        } else if (arg == "--no-plot") {
            show_plots = false;
        } else if (arg == "--validate") {
            run_validation = true;
        } else if (arg == "--quick") {
            quick_mode = true;
        } else if (arg == "--seed" && i + 1 < argc) {
            seed = static_cast<unsigned>(std::stoi(argv[++i]));
        } else if (arg == "--noise" && i + 1 < argc) {
            noise_override = std::stod(argv[++i]);
        } else if (arg == "--traj" && i + 1 < argc) {
            traj_override = std::stoi(argv[++i]);
        } else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
    }

    // Obtener configuración del sistema
    auto system = UE::SystemFactory<Config>::getSystem(system_name);

    // Ejecutar simulación
    runSimulation(
        system,
        show_plots,
        run_validation,
        quick_mode,
        seed,
        noise_override,
        traj_override
    );

    return 0;
}
