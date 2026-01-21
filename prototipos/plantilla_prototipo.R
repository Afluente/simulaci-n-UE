# =============================================================================
# Plantilla para prototipos UE en R
# =============================================================================
# Usar este archivo como base para probar nuevos sistemas/observables
# Una vez validado, traducir a C++ en src/
# =============================================================================

library(tidyverse)

# --- PARAMETROS DEL SISTEMA ---
params <- list(
  # Espacio de estados
  x_min = -2,
  x_max = 2,

  # Dinamica
  dt = 0.01,
  noise_sigma = 0.3,

  # Condicion ERS (R)
  theta_threshold = 0.5,
  tau_persistence = 50,

  # Simulacion
  n_trajectories = 1000,
  trajectory_length = 500,
  seed = 42
)

# --- POTENCIAL V(x) ---
# Modificar segun el sistema a probar
potential <- function(x) {
  # Ejemplo: doble pozo
  (x^2 - 1)^2
}

# Gradiente del potencial
grad_potential <- function(x) {
  4 * x * (x^2 - 1)
}

# --- OBSERVABLE theta(x) ---
# Modificar segun el observable a probar
observable <- function(x) {
  # Ejemplo: distancia al origen
  abs(x)
}

# --- DINAMICA (Langevin) ---
step <- function(x, dt, sigma) {
  x - grad_potential(x) * dt + sigma * sqrt(dt) * rnorm(1)
}

# --- SIMULACION DE UNA TRAYECTORIA ---
simulate_trajectory <- function(x0, length, dt, sigma) {
  x <- numeric(length)
  x[1] <- x0
  for (t in 2:length) {
    x[t] <- step(x[t-1], dt, sigma)
  }
  x
}

# --- VERIFICAR CONDICION R ---
check_ERS <- function(traj, theta_thresh, tau_p, above = TRUE) {
  theta_vals <- observable(traj)
  if (above) {
    condition <- theta_vals >= theta_thresh
  } else {
    condition <- theta_vals <= theta_thresh
  }

  # Buscar secuencia continua de tau_p pasos
  rle_result <- rle(condition)
  any(rle_result$lengths[rle_result$values] >= tau_p)
}

# --- EJECUTAR SIMULACION ---
set.seed(params$seed)

cat("Simulando", params$n_trajectories, "trayectorias...\n")

trajectories <- list()
conditioned_indices <- c()

for (i in 1:params$n_trajectories) {
  x0 <- runif(1, params$x_min, params$x_max)
  traj <- simulate_trajectory(x0, params$trajectory_length, params$dt, params$noise_sigma)
  trajectories[[i]] <- traj

  if (check_ERS(traj, params$theta_threshold, params$tau_persistence)) {
    conditioned_indices <- c(conditioned_indices, i)
  }
}

cat("Trayectorias bajo R:", length(conditioned_indices),
    "(", round(100 * length(conditioned_indices) / params$n_trajectories, 1), "%)\n")

# --- ANALISIS BASICO ---
if (length(conditioned_indices) > 0) {
  # Extraer valores x de trayectorias condicionadas
  x_conditioned <- unlist(trajectories[conditioned_indices])

  # Histograma
  cat("\nDistribucion bajo R:\n")
  print(summary(x_conditioned))

  # Visualizacion rapida
  hist(x_conditioned, breaks = 50, main = "Distribucion bajo R",
       xlab = "x", col = "steelblue", border = "white")

  # Paisaje de costo (aproximado)
  h <- hist(x_conditioned, breaks = 100, plot = FALSE)
  freq <- h$counts / sum(h$counts)
  freq[freq == 0] <- min(freq[freq > 0]) / 10  # evitar log(0)
  cost <- -log(freq)

  plot(h$mids, cost, type = "l", lwd = 2,
       main = "Paisaje de costo F_R(x)", xlab = "x", ylab = "F_R")
}

cat("\n--- Prototipo completado ---\n")
cat("Si los resultados son correctos, traducir a C++\n")
