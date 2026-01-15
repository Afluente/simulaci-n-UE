#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace UE {

/**
 * Visualization - Visualización de Resultados UE
 *
 * Genera gráficos usando gnuplot:
 * - Paisaje de costo F_n^R(ψ)
 * - Distribución de probabilidad p_n(ψ|R)
 * - Trayectorias de muestra
 * - Histograma de configuraciones
 * - Posiciones de clases metaestables
 */
class Visualization {
private:
    std::string output_dir_;
    bool interactive_;

    // Ejecutar comando gnuplot
    void runGnuplot(const std::string& script) const {
        std::string script_file = output_dir_ + "/plot_script.gp";
        std::ofstream out(script_file);
        out << script;
        out.close();

        std::string cmd = "gnuplot ";
        if (interactive_) {
            cmd += "-persist ";
        }
        cmd += script_file + " 2>/dev/null";
        std::system(cmd.c_str());
    }

public:
    Visualization(const std::string& output_dir = ".", bool interactive = true)
        : output_dir_(output_dir)
        , interactive_(interactive)
    {}

    void setOutputDir(const std::string& dir) { output_dir_ = dir; }
    void setInteractive(bool interactive) { interactive_ = interactive; }

    /**
     * Graficar paisaje de costo F_n^R y probabilidad p_n(ψ|R)
     */
    template<typename Config>
    void plotLandscape(
        const std::vector<std::pair<Config, double>>& probability_data,
        const std::vector<std::pair<Config, double>>& cost_data,
        const std::vector<Config>& class_centers,
        const std::string& title = "Paisaje Condicionado F_n^R") const
    {
        // Guardar datos de probabilidad
        std::string prob_file = output_dir_ + "/probability.dat";
        std::ofstream prob_out(prob_file);
        for (const auto& [config, prob] : probability_data) {
            prob_out << config[0] << " " << prob << "\n";
        }
        prob_out.close();

        // Guardar datos de costo
        std::string cost_file = output_dir_ + "/cost.dat";
        std::ofstream cost_out(cost_file);
        for (const auto& [config, cost] : cost_data) {
            cost_out << config[0] << " " << cost << "\n";
        }
        cost_out.close();

        // Guardar posiciones de clases
        std::string class_file = output_dir_ + "/classes.dat";
        std::ofstream class_out(class_file);
        for (size_t i = 0; i < class_centers.size(); ++i) {
            class_out << class_centers[i][0] << " " << i << "\n";
        }
        class_out.close();

        // Script de gnuplot
        std::ostringstream script;
        script << "set terminal qt size 900,600 font 'Sans,12'\n";
        script << "set title '" << title << "' font 'Sans,14'\n";
        script << "set xlabel 'Configuracion ψ' font 'Sans,12'\n";
        script << "set grid\n";
        script << "set key top right\n";
        script << "\n";
        script << "set multiplot layout 2,1\n";
        script << "\n";
        script << "# Panel superior: Probabilidad p_n(ψ|R)\n";
        script << "set ylabel 'p_n(ψ|R)' font 'Sans,12'\n";
        script << "set style fill solid 0.3\n";
        script << "plot '" << prob_file << "' with filledcurves y1=0 lc rgb '#3366cc' title 'Probabilidad', \\\n";
        script << "     '" << prob_file << "' with lines lw 2 lc rgb '#0044aa' notitle";

        if (!class_centers.empty()) {
            script << ", \\\n     '" << class_file << "' using 1:(0):(0):(0.5) with vectors nohead lw 2 lc rgb '#cc0000' title 'Clases A_k'";
        }
        script << "\n\n";

        script << "# Panel inferior: Costo F_n^R\n";
        script << "set ylabel 'F_n^R[ψ]' font 'Sans,12'\n";
        script << "plot '" << cost_file << "' with lines lw 2 lc rgb '#cc6600' title 'Costo (paisaje)'";

        if (!class_centers.empty()) {
            script << ", \\\n     '" << class_file << "' using 1:(0) with points pt 7 ps 2 lc rgb '#cc0000' title 'Minimos (clases)'";
        }
        script << "\n\n";

        script << "unset multiplot\n";

        runGnuplot(script.str());
    }

    /**
     * Graficar histograma de configuraciones bajo R
     */
    void plotHistogram(
        const std::vector<double>& values,
        int num_bins,
        const std::string& title = "Distribucion de Configuraciones bajo R",
        const std::string& xlabel = "Configuracion x") const
    {
        if (values.empty()) return;

        // Calcular histograma
        double min_val = values[0], max_val = values[0];
        for (double v : values) {
            min_val = std::min(min_val, v);
            max_val = std::max(max_val, v);
        }

        double range = max_val - min_val;
        if (range < 1e-10) range = 1.0;
        double bin_width = range / num_bins;

        std::vector<int> counts(num_bins, 0);
        for (double v : values) {
            int bin = static_cast<int>((v - min_val) / bin_width);
            bin = std::max(0, std::min(num_bins - 1, bin));
            counts[bin]++;
        }

        // Guardar datos
        std::string hist_file = output_dir_ + "/histogram.dat";
        std::ofstream out(hist_file);
        for (int i = 0; i < num_bins; ++i) {
            double x = min_val + (i + 0.5) * bin_width;
            double freq = static_cast<double>(counts[i]) / values.size();
            out << x << " " << freq << "\n";
        }
        out.close();

        // Script
        std::ostringstream script;
        script << "set terminal qt size 800,500 font 'Sans,12'\n";
        script << "set title '" << title << "' font 'Sans,14'\n";
        script << "set xlabel '" << xlabel << "' font 'Sans,12'\n";
        script << "set ylabel 'Frecuencia relativa' font 'Sans,12'\n";
        script << "set grid\n";
        script << "set style fill solid 0.6\n";
        script << "set boxwidth " << (bin_width * 0.9) << "\n";
        script << "plot '" << hist_file << "' with boxes lc rgb '#3366cc' title 'p_n(x|R)'\n";

        runGnuplot(script.str());
    }

    /**
     * Graficar trayectorias de muestra
     */
    template<typename Config>
    void plotTrajectories(
        const std::vector<std::vector<Config>>& trajectories,
        int max_trajectories = 10,
        const std::string& title = "Trayectorias bajo Condicion R") const
    {
        if (trajectories.empty()) return;

        int num_to_plot = std::min(static_cast<int>(trajectories.size()), max_trajectories);

        // Guardar datos de trayectorias
        std::vector<std::string> traj_files;
        for (int i = 0; i < num_to_plot; ++i) {
            std::string file = output_dir_ + "/traj_" + std::to_string(i) + ".dat";
            traj_files.push_back(file);

            std::ofstream out(file);
            for (size_t t = 0; t < trajectories[i].size(); ++t) {
                out << t << " " << trajectories[i][t][0] << "\n";
            }
            out.close();
        }

        // Script
        std::ostringstream script;
        script << "set terminal qt size 900,500 font 'Sans,12'\n";
        script << "set title '" << title << "' font 'Sans,14'\n";
        script << "set xlabel 'Tiempo t' font 'Sans,12'\n";
        script << "set ylabel 'Configuracion x(t)' font 'Sans,12'\n";
        script << "set grid\n";
        script << "set key outside right\n";

        // Líneas horizontales para los pozos
        script << "set arrow from graph 0, first 1 to graph 1, first 1 nohead lc rgb '#999999' dt 2\n";
        script << "set arrow from graph 0, first -1 to graph 1, first -1 nohead lc rgb '#999999' dt 2\n";
        script << "set label 'Pozo +1' at graph 0.02, first 1.2 font 'Sans,10'\n";
        script << "set label 'Pozo -1' at graph 0.02, first -0.8 font 'Sans,10'\n";

        script << "plot ";
        for (int i = 0; i < num_to_plot; ++i) {
            if (i > 0) script << ", \\\n     ";
            script << "'" << traj_files[i] << "' with lines lw 1.5 title 'Traj " << (i+1) << "'";
        }
        script << "\n";

        runGnuplot(script.str());
    }

    /**
     * Graficar potencial original vs paisaje condicionado
     */
    template<typename PotentialFunc>
    void plotPotentialComparison(
        PotentialFunc potential,
        const std::vector<std::pair<std::vector<double>, double>>& conditioned_cost,
        double x_min, double x_max, int num_points,
        const std::string& title = "Potencial Original vs Paisaje Condicionado") const
    {
        // Datos del potencial original
        std::string pot_file = output_dir_ + "/potential.dat";
        std::ofstream pot_out(pot_file);
        double dx = (x_max - x_min) / (num_points - 1);
        for (int i = 0; i < num_points; ++i) {
            double x = x_min + i * dx;
            std::vector<double> config = {x};
            pot_out << x << " " << potential(config) << "\n";
        }
        pot_out.close();

        // Datos del paisaje condicionado
        std::string cond_file = output_dir_ + "/conditioned.dat";
        std::ofstream cond_out(cond_file);
        for (const auto& [config, cost] : conditioned_cost) {
            cond_out << config[0] << " " << cost << "\n";
        }
        cond_out.close();

        // Script
        std::ostringstream script;
        script << "set terminal qt size 800,500 font 'Sans,12'\n";
        script << "set title '" << title << "' font 'Sans,14'\n";
        script << "set xlabel 'Configuracion x' font 'Sans,12'\n";
        script << "set ylabel 'Energia / Costo' font 'Sans,12'\n";
        script << "set grid\n";
        script << "set key top center\n";
        script << "plot '" << pot_file << "' with lines lw 2 lc rgb '#0066cc' title 'V(x) original', \\\n";
        script << "     '" << cond_file << "' with lines lw 2 lc rgb '#cc6600' title 'F_n^R[x] condicionado'\n";

        runGnuplot(script.str());
    }

    /**
     * Graficar todo en una figura combinada
     */
    template<typename Config>
    void plotCombined(
        const std::vector<double>& x_values,  // Valores x del histograma
        const std::vector<std::pair<Config, double>>& probability,
        const std::vector<std::pair<Config, double>>& cost,
        const std::vector<Config>& class_centers,
        const std::vector<std::vector<Config>>& sample_trajectories,
        const std::string& title = "Resultados Simulacion UE") const
    {
        // Guardar todos los datos
        // Histograma
        std::string hist_file = output_dir_ + "/hist_combined.dat";
        {
            std::ofstream out(hist_file);
            int num_bins = 50;
            double min_val = -3.0, max_val = 3.0;
            double bin_width = (max_val - min_val) / num_bins;
            std::vector<int> counts(num_bins, 0);

            for (double v : x_values) {
                int bin = static_cast<int>((v - min_val) / bin_width);
                bin = std::max(0, std::min(num_bins - 1, bin));
                counts[bin]++;
            }

            for (int i = 0; i < num_bins; ++i) {
                double x = min_val + (i + 0.5) * bin_width;
                double freq = x_values.empty() ? 0 : static_cast<double>(counts[i]) / x_values.size();
                out << x << " " << freq << "\n";
            }
        }

        // Probabilidad
        std::string prob_file = output_dir_ + "/prob_combined.dat";
        {
            std::ofstream out(prob_file);
            for (const auto& [config, prob] : probability) {
                out << config[0] << " " << prob << "\n";
            }
        }

        // Costo
        std::string cost_file = output_dir_ + "/cost_combined.dat";
        {
            std::ofstream out(cost_file);
            for (const auto& [config, c] : cost) {
                out << config[0] << " " << c << "\n";
            }
        }

        // Potencial original
        std::string pot_file = output_dir_ + "/pot_combined.dat";
        {
            std::ofstream out(pot_file);
            for (double x = -3.0; x <= 3.0; x += 0.05) {
                double v = (x*x - 1.0) * (x*x - 1.0);  // V(x) = (x²-1)²
                out << x << " " << v << "\n";
            }
        }

        // Clases
        std::string class_file = output_dir_ + "/class_combined.dat";
        {
            std::ofstream out(class_file);
            for (const auto& center : class_centers) {
                out << center[0] << "\n";
            }
        }

        // Trayectorias
        int num_traj = std::min(5, static_cast<int>(sample_trajectories.size()));
        std::vector<std::string> traj_files;
        for (int i = 0; i < num_traj; ++i) {
            std::string file = output_dir_ + "/traj_comb_" + std::to_string(i) + ".dat";
            traj_files.push_back(file);
            std::ofstream out(file);
            for (size_t t = 0; t < sample_trajectories[i].size(); ++t) {
                out << t << " " << sample_trajectories[i][t][0] << "\n";
            }
        }

        // Script combinado
        std::ostringstream script;
        script << "set terminal qt size 1200,800 font 'Sans,11'\n";
        script << "set multiplot layout 2,2 title '" << title << "' font 'Sans,14'\n";
        script << "\n";

        // Panel 1: Histograma
        script << "set title 'Distribucion p_n(x|R)'\n";
        script << "set xlabel 'x'\n";
        script << "set ylabel 'Frecuencia'\n";
        script << "set grid\n";
        script << "set style fill solid 0.5\n";
        script << "set boxwidth 0.1\n";
        script << "plot '" << hist_file << "' with boxes lc rgb '#3366cc' notitle\n";
        script << "\n";

        // Panel 2: Potencial y Paisaje
        script << "set title 'Potencial V(x) vs Paisaje F_n^R'\n";
        script << "set xlabel 'x'\n";
        script << "set ylabel 'Energia'\n";
        script << "set style fill transparent\n";
        script << "plot '" << pot_file << "' with lines lw 2 lc rgb '#666666' title 'V(x)', \\\n";
        script << "     '" << cost_file << "' with lines lw 2 lc rgb '#cc6600' title 'F_n^R'\n";
        script << "\n";

        // Panel 3: Probabilidad con clases
        script << "set title 'Probabilidad y Clases Metaestables'\n";
        script << "set xlabel 'x'\n";
        script << "set ylabel 'p_n(x|R)'\n";
        script << "set style fill solid 0.3\n";
        script << "plot '" << prob_file << "' with filledcurves y1=0 lc rgb '#3366cc' title 'p_n(x|R)'";
        if (!class_centers.empty()) {
            script << ", \\\n     '" << class_file << "' using 1:(0) with impulses lw 3 lc rgb '#cc0000' title 'Clases'";
        }
        script << "\n\n";

        // Panel 4: Trayectorias
        script << "set title 'Trayectorias de Muestra'\n";
        script << "set xlabel 'Tiempo t'\n";
        script << "set ylabel 'x(t)'\n";
        script << "set style fill transparent\n";
        script << "set arrow from graph 0, first 1 to graph 1, first 1 nohead lc rgb '#aaaaaa' dt 2\n";
        script << "set arrow from graph 0, first -1 to graph 1, first -1 nohead lc rgb '#aaaaaa' dt 2\n";

        if (!traj_files.empty()) {
            script << "plot ";
            for (size_t i = 0; i < traj_files.size(); ++i) {
                if (i > 0) script << ", \\\n     ";
                script << "'" << traj_files[i] << "' with lines lw 1.2 notitle";
            }
            script << "\n";
        }

        script << "\nunset multiplot\n";

        runGnuplot(script.str());
    }

    /**
     * Verificar si gnuplot está disponible
     */
    static bool isGnuplotAvailable() {
        return std::system("which gnuplot > /dev/null 2>&1") == 0;
    }

    /**
     * Gráfico ASCII - Alternativa cuando no hay gnuplot
     */
    void plotASCII(
        const std::vector<double>& x_values,
        const std::vector<double>& y_values,
        const std::string& title,
        int width = 70,
        int height = 20) const
    {
        if (x_values.empty() || y_values.empty()) return;

        // Encontrar rangos
        double x_min = *std::min_element(x_values.begin(), x_values.end());
        double x_max = *std::max_element(x_values.begin(), x_values.end());
        double y_min = *std::min_element(y_values.begin(), y_values.end());
        double y_max = *std::max_element(y_values.begin(), y_values.end());

        if (y_max - y_min < 1e-10) { y_max = y_min + 1.0; }

        // Crear canvas
        std::vector<std::string> canvas(height, std::string(width, ' '));

        // Dibujar puntos
        for (size_t i = 0; i < x_values.size() && i < y_values.size(); ++i) {
            int col = static_cast<int>((x_values[i] - x_min) / (x_max - x_min) * (width - 1));
            int row = static_cast<int>((y_max - y_values[i]) / (y_max - y_min) * (height - 1));
            col = std::max(0, std::min(width - 1, col));
            row = std::max(0, std::min(height - 1, row));
            canvas[row][col] = '*';
        }

        // Imprimir
        std::cout << "\n+" << std::string(width, '-') << "+\n";
        std::cout << "| " << title << std::string(std::max(0, width - 2 - static_cast<int>(title.length())), ' ') << "|\n";
        std::cout << "+" << std::string(width, '-') << "+\n";

        for (int row = 0; row < height; ++row) {
            double y_val = y_max - (row * (y_max - y_min) / (height - 1));
            std::cout << "|" << canvas[row] << "| " << std::fixed << std::setprecision(2) << y_val << "\n";
        }

        std::cout << "+" << std::string(width, '-') << "+\n";
        std::cout << " " << std::fixed << std::setprecision(2) << x_min;
        std::cout << std::string(width / 2 - 6, ' ') << "x";
        std::cout << std::string(width / 2 - 6, ' ') << x_max << "\n";
    }

    /**
     * Histograma ASCII
     */
    void plotHistogramASCII(
        const std::vector<double>& values,
        int num_bins = 40,
        const std::string& title = "Histograma",
        int height = 15) const
    {
        if (values.empty()) return;

        double min_val = *std::min_element(values.begin(), values.end());
        double max_val = *std::max_element(values.begin(), values.end());
        double range = max_val - min_val;
        if (range < 1e-10) range = 1.0;
        double bin_width = range / num_bins;

        // Calcular conteos
        std::vector<int> counts(num_bins, 0);
        for (double v : values) {
            int bin = static_cast<int>((v - min_val) / bin_width);
            bin = std::max(0, std::min(num_bins - 1, bin));
            counts[bin]++;
        }

        int max_count = *std::max_element(counts.begin(), counts.end());
        if (max_count == 0) max_count = 1;

        // Imprimir título
        std::cout << "\n" << title << "\n";
        std::cout << std::string(num_bins + 2, '=') << "\n";

        // Imprimir histograma
        for (int row = height; row >= 1; --row) {
            double threshold = static_cast<double>(row) / height;
            std::cout << "|";
            for (int bin = 0; bin < num_bins; ++bin) {
                double normalized = static_cast<double>(counts[bin]) / max_count;
                if (normalized >= threshold) {
                    std::cout << "#";
                } else if (normalized >= threshold - 0.5/height) {
                    std::cout << ".";
                } else {
                    std::cout << " ";
                }
            }
            std::cout << "|\n";
        }

        // Eje X
        std::cout << "+" << std::string(num_bins, '-') << "+\n";
        std::cout << std::fixed << std::setprecision(1);
        std::cout << min_val;
        std::cout << std::string(num_bins / 2 - 4, ' ') << "x";
        std::cout << std::string(num_bins / 2 - 4, ' ') << max_val << "\n";
    }

    /**
     * Gráfico combinado ASCII
     */
    template<typename Config>
    void plotCombinedASCII(
        const std::vector<double>& x_values,
        const std::vector<std::pair<Config, double>>& probability,
        const std::vector<std::pair<Config, double>>& cost,
        const std::vector<Config>& class_centers) const
    {
        std::cout << "\n";
        std::cout << "+================================================================+\n";
        std::cout << "|           VISUALIZACION UE (Modo ASCII)                        |\n";
        std::cout << "+================================================================+\n";

        // 1. Histograma de configuraciones
        if (!x_values.empty()) {
            plotHistogramASCII(x_values, 50, "DISTRIBUCION p_n(x|R)", 12);
        }

        // 2. Paisaje de costo
        if (!cost.empty()) {
            std::vector<double> x_cost, y_cost;
            for (const auto& [config, c] : cost) {
                x_cost.push_back(config[0]);
                y_cost.push_back(c);
            }
            plotASCII(x_cost, y_cost, "PAISAJE F_n^R[x]", 60, 15);
        }

        // 3. Distribución de probabilidad
        if (!probability.empty()) {
            std::vector<double> x_prob, y_prob;
            for (const auto& [config, p] : probability) {
                x_prob.push_back(config[0]);
                y_prob.push_back(p);
            }
            plotASCII(x_prob, y_prob, "PROBABILIDAD p_n(x|R)", 60, 12);
        }

        // 4. Mostrar centros de clases
        if (!class_centers.empty()) {
            std::cout << "\nCLASES METAESTABLES:\n";
            std::cout << "-------------------\n";
            for (size_t i = 0; i < class_centers.size(); ++i) {
                std::cout << "  Clase " << i << ": x = " << std::fixed
                          << std::setprecision(3) << class_centers[i][0] << "\n";
            }
        }

        std::cout << "\n";
    }
};

} // namespace UE
