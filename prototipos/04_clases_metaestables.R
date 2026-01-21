# =============================================================================
# CONCEPTO 4: Clases Metaestables (A_k)
# =============================================================================
#
# Definicion en dos capas:
#
# CAPA 1 - Geometrica (Cuencas):
#   A_k = { ψ : flujo_gradiente(ψ) → minimo_k de F_R }
#   Cada estado pertenece a la cuenca del minimo al que "cae"
#
# CAPA 2 - Dinamica (Metaestabilidad):
#   τ_relax(A_k) << E[τ_exit(A_k)]
#   - τ_relax: tiempo de relajacion INTERNA (equilibrarse dentro de la clase)
#   - τ_exit: tiempo de ESCAPE (salir de la clase)
#   Una clase es metaestable si se equilibra rapido internamente
#   pero tarda mucho en escapar.
#
# =============================================================================

cat("=== PRUEBA: Clases Metaestables ===\n\n")

# --- PARAMETROS ---
set.seed(42)
n_traj <- 2000
traj_length <- 1000
dt <- 0.01
noise <- 0.35

# --- POTENCIAL TRIPLE POZO ---
V <- function(x) {
  # Tres minimos: x ≈ -1.5, 0, 1.5
  0.5 * (x^2 - 2)^2 - 0.5 * exp(-3 * x^2)
}

grad_V <- function(x) {
  2 * x * (x^2 - 2) + 3 * x * exp(-3 * x^2)
}

# --- DINAMICA ---
step <- function(x, dt, sigma) {
  x - grad_V(x) * dt + sigma * sqrt(dt) * rnorm(1)
}

simulate <- function(x0, n, dt, sigma) {
  x <- numeric(n)
  x[1] <- x0
  for (t in 2:n) {
    x[t] <- step(x[t-1], dt, sigma)
    x[t] <- max(-3, min(3, x[t]))
  }
  x
}

# --- CONDICION R: energia moderada ---
observable <- function(x) V(x)
theta_threshold <- 0.3
tau_persistence <- 20

cumple_ERS <- function(traj, theta_star, tau_star) {
  theta_t <- observable(traj)
  condicion <- theta_t >= theta_star
  if (!any(condicion)) return(FALSE)
  rachas <- rle(condicion)
  duraciones <- rachas$lengths[rachas$values == TRUE]
  if (length(duraciones) == 0) return(FALSE)
  max(duraciones) >= tau_star
}

# --- SIMULAR ---
cat("Simulando", n_traj, "trayectorias en triple pozo...\n")

trayectorias <- vector("list", n_traj)
cumple_R <- logical(n_traj)

for (i in 1:n_traj) {
  x0 <- runif(1, -2.5, 2.5)
  trayectorias[[i]] <- simulate(x0, traj_length, dt, noise)
  cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_threshold, tau_persistence)
}

n_cond <- sum(cumple_R)
cat("Condicionadas:", n_cond, "\n\n")

# --- CONSTRUIR PAISAJE F_R ---
x_cond <- unlist(trayectorias[cumple_R])

n_bins <- 100
breaks <- seq(-3, 3, length.out = n_bins + 1)
mids <- (breaks[-1] + breaks[-length(breaks)]) / 2

h <- hist(x_cond, breaks = breaks, plot = FALSE)
freq <- h$counts / sum(h$counts)
freq[freq == 0] <- min(freq[freq > 0]) / 100

epsilon <- 1 / log(length(x_cond))
F_R <- -epsilon * log(freq)
F_R <- F_R - min(F_R)

# --- IDENTIFICAR MINIMOS (centros de clase) ---
find_local_minima <- function(y, threshold = 0.01) {
  n <- length(y)
  minima <- c()
  for (i in 3:(n-2)) {
    if (y[i] < y[i-1] && y[i] < y[i+1] &&
        y[i] < y[i-2] && y[i] < y[i+2]) {
      minima <- c(minima, i)
    }
  }
  # Filtrar minimos muy cercanos
  if (length(minima) > 1) {
    filtered <- minima[1]
    for (m in minima[-1]) {
      if (abs(mids[m] - mids[filtered[length(filtered)]]) > 0.5) {
        filtered <- c(filtered, m)
      }
    }
    minima <- filtered
  }
  minima
}

minima_idx <- find_local_minima(F_R)
class_centers <- mids[minima_idx]
n_classes <- length(class_centers)

cat("Clases detectadas:", n_classes, "\n")
cat("Centros:", round(class_centers, 2), "\n\n")

# --- ASIGNAR ESTADOS A CLASES (descenso de gradiente) ---
assign_class <- function(x, centers) {
  which.min(abs(x - centers))
}

class_assignment <- sapply(x_cond, assign_class, centers = class_centers)

# --- CALCULAR TIEMPOS CARACTERISTICOS ---
cat("--- ANALISIS DE METAESTABILIDAD ---\n\n")

# Para cada trayectoria condicionada, calcular:
# - Tiempo en cada clase
# - Transiciones entre clases

analyze_trajectory <- function(traj, centers) {
  classes <- sapply(traj, assign_class, centers = centers)

  # Tiempo en cada clase
  time_in_class <- table(factor(classes, levels = 1:length(centers)))

  # Contar transiciones
  transitions <- sum(diff(classes) != 0)

  # Tiempo medio de permanencia (proxy de τ_exit)
  if (transitions > 0) {
    mean_stay <- length(traj) / transitions
  } else {
    mean_stay <- length(traj)
  }

  list(time_in_class = as.numeric(time_in_class),
       transitions = transitions,
       mean_stay = mean_stay)
}

# Analizar trayectorias condicionadas
traj_cond <- trayectorias[cumple_R]
analyses <- lapply(traj_cond, analyze_trajectory, centers = class_centers)

total_transitions <- sum(sapply(analyses, function(a) a$transitions))
mean_stay_times <- sapply(analyses, function(a) a$mean_stay)

cat("Transiciones totales entre clases:", total_transitions, "\n")
cat("Tiempo medio de permanencia (τ_exit proxy):", round(mean(mean_stay_times), 1), "pasos\n")

# Tiempo de relajacion aproximado (autocorrelacion)
calc_relax_time <- function(x) {
  n <- length(x)
  if (n < 50) return(NA)
  acf_vals <- acf(x, lag.max = min(100, n/2), plot = FALSE)$acf
  # Tiempo donde cae a 1/e
  idx <- which(acf_vals < exp(-1))[1]
  if (is.na(idx)) return(length(acf_vals))
  idx
}

relax_times <- sapply(traj_cond[1:min(100, length(traj_cond))], calc_relax_time)
mean_relax <- mean(relax_times, na.rm = TRUE)

cat("Tiempo de relajacion (τ_relax proxy):", round(mean_relax, 1), "pasos\n")
cat("Ratio τ_exit/τ_relax:", round(mean(mean_stay_times)/mean_relax, 1), "\n")

if (mean(mean_stay_times) > 5 * mean_relax) {
  cat("\n✓ Metaestabilidad verificada: τ_relax << τ_exit\n")
} else {
  cat("\n⚠ Metaestabilidad debil: τ_relax no es mucho menor que τ_exit\n")
}

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Potencial
x_grid <- seq(-3, 3, length.out = 200)
plot(x_grid, V(x_grid), type = "l", lwd = 2, col = "darkred",
     main = "Potencial V(x) - Triple pozo",
     xlab = "x", ylab = "V(x)")

# 2. Paisaje F_R con clases
plot(mids, F_R, type = "l", lwd = 2, col = "darkgreen",
     main = "Paisaje F_R y clases A_k",
     xlab = "x", ylab = "F_R(x)")
points(class_centers, F_R[minima_idx], col = "red", pch = 19, cex = 2)
for (i in 1:n_classes) {
  text(class_centers[i], F_R[minima_idx[i]] - 0.015,
       paste0("A", i), col = "red", font = 2)
}

# 3. Histograma coloreado por clase
colors <- c("steelblue", "coral", "mediumseagreen", "purple")[1:n_classes]
hist(x_cond, breaks = 60, col = "gray80", border = "white",
     main = "Distribucion por clases",
     xlab = "x", xlim = c(-3, 3))
for (i in 1:n_classes) {
  x_class <- x_cond[class_assignment == i]
  hist(x_class, breaks = 60, col = adjustcolor(colors[i], 0.6),
       border = "white", add = TRUE)
}
legend("topright", paste0("A", 1:n_classes), fill = colors)

# 4. Trayectoria ejemplo mostrando transiciones
traj_ejemplo <- traj_cond[[which.max(sapply(analyses, function(a) a$transitions))]]
plot(traj_ejemplo, type = "l", col = "gray50",
     main = "Trayectoria con transiciones",
     xlab = "tiempo", ylab = "x(t)")
abline(h = class_centers, col = colors, lty = 2)

par(mfrow = c(1, 1))

# --- PESOS DE CLASE ---
cat("\n--- PESOS DE CLASE P(A_k|R) ---\n")
class_weights <- table(class_assignment) / length(class_assignment)
for (i in 1:n_classes) {
  cat("  P(A", i, "|R) =", round(class_weights[i], 4), "\n", sep = "")
}

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("Las clases metaestables A_k son:\n")
cat("  1. Cuencas geometricas del paisaje F_R (donde 'cae' cada estado)\n")
cat("  2. Regiones donde el sistema se equilibra rapido (τ_relax bajo)\n")
cat("  3. Pero escapa lento (τ_exit alto)\n")
cat("  4. El ratio τ_exit/τ_relax mide la 'calidad' de la metaestabilidad\n")
