# =============================================================================
# CONCEPTO 3: Paisaje de Costo F_R (Funcion de Costo Condicionada)
# =============================================================================
#
# Definicion:
#   F_n^R[ψ] = -ε_n · log p_n(ψ|R)
#
# donde ε_n → 0 es parametro de escala (tipicamente 1/n o 1/log(n))
#
# Interpretacion:
#   - Estados con alta probabilidad bajo R → bajo costo
#   - Estados con baja probabilidad bajo R → alto costo
#   - Los MINIMOS de F_R son los estados TIPICOS bajo R
#   - F_R NO es energia fisica, es COSTO INFORMACIONAL
#
# Relacion con Grandes Desviaciones:
#   p_n(ψ|R) ≈ exp(-F_n^R[ψ] / ε_n)
# =============================================================================

cat("=== PRUEBA: Paisaje de Costo F_R ===\n\n")

# --- PARAMETROS ---
set.seed(42)
n_traj <- 3000
traj_length <- 500
dt <- 0.01
noise <- 0.4

# --- POTENCIAL ASIMETRICO ---
# Pozo izquierdo mas profundo que el derecho
V <- function(x) (x^2 - 1)^2 - 0.3 * x
grad_V <- function(x) 4 * x * (x^2 - 1) - 0.3

# --- DINAMICA ---
step <- function(x, dt, sigma) {
  x - grad_V(x) * dt + sigma * sqrt(dt) * rnorm(1)
}

simulate <- function(x0, n, dt, sigma) {
  x <- numeric(n)
  x[1] <- x0
  for (t in 2:n) {
    x[t] <- step(x[t-1], dt, sigma)
    # Limitar a rango razonable
    x[t] <- max(-3, min(3, x[t]))
  }
  x
}

# --- CONDICION R: energia alta sostenida ---
observable <- function(x) V(x)
theta_threshold <- 0.5  # energia por encima de 0.5
tau_persistence <- 25

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
cat("Simulando", n_traj, "trayectorias en pozo asimetrico...\n")
cat("Condicion R: V(x) >=", theta_threshold, "durante >=", tau_persistence, "pasos\n\n")

trayectorias <- vector("list", n_traj)
cumple_R <- logical(n_traj)

for (i in 1:n_traj) {
  x0 <- runif(1, -2, 2)
  trayectorias[[i]] <- simulate(x0, traj_length, dt, noise)
  cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_threshold, tau_persistence)
}

n_cond <- sum(cumple_R)
cat("Trayectorias condicionadas:", n_cond, "(", round(100*n_cond/n_traj, 1), "%)\n\n")

# --- CONSTRUIR PAISAJE DE COSTO ---
x_cond <- unlist(trayectorias[cumple_R])
n_total <- length(x_cond)

# Discretizar espacio
n_bins <- 100
x_min <- -2.5
x_max <- 2.5
breaks <- seq(x_min, x_max, length.out = n_bins + 1)
mids <- (breaks[-1] + breaks[-length(breaks)]) / 2

# Histograma = frecuencias
h <- hist(x_cond, breaks = breaks, plot = FALSE)
freq <- h$counts / sum(h$counts)

# Evitar log(0)
freq[freq == 0] <- min(freq[freq > 0]) / 100

# Paisaje de costo: F_R = -epsilon * log(p)
epsilon <- 1 / log(n_total)  # escala tipica
F_R <- -epsilon * log(freq)

# Normalizar para visualizacion
F_R <- F_R - min(F_R)

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Potencial fisico V(x)
x_grid <- seq(-2.5, 2.5, length.out = 200)
plot(x_grid, V(x_grid), type = "l", lwd = 2, col = "darkred",
     main = "Potencial fisico V(x)",
     xlab = "x", ylab = "V(x)")
abline(h = theta_threshold, col = "orange", lty = 2)
text(2, theta_threshold + 0.1, "θ*", col = "orange")

# 2. Distribucion P(x|R)
hist(x_cond, breaks = 60, prob = TRUE, col = "steelblue", border = "white",
     main = "Distribucion P(x|R)",
     xlab = "x", xlim = c(-2.5, 2.5))

# 3. Paisaje de costo F_R(x)
plot(mids, F_R, type = "l", lwd = 2, col = "darkgreen",
     main = "Paisaje de Costo F_R(x) = -ε·log P(x|R)",
     xlab = "x", ylab = "F_R(x)")
# Marcar minimos
min_idx <- which(F_R == min(F_R))
points(mids[min_idx], F_R[min_idx], col = "red", pch = 19, cex = 1.5)
text(mids[min_idx], F_R[min_idx] + 0.02, "minimo", col = "red", pos = 3)

# 4. Comparacion V(x) vs F_R(x)
# Normalizar ambos para comparar formas
V_norm <- V(mids)
V_norm <- (V_norm - min(V_norm)) / max(V_norm - min(V_norm))
F_norm <- F_R / max(F_R)

plot(mids, V_norm, type = "l", lwd = 2, col = "darkred",
     main = "Comparacion: V(x) vs F_R(x)",
     xlab = "x", ylab = "valor normalizado", ylim = c(0, 1.1))
lines(mids, F_norm, lwd = 2, col = "darkgreen")
legend("topright", c("V(x) potencial", "F_R(x) costo"),
       col = c("darkred", "darkgreen"), lwd = 2)

par(mfrow = c(1, 1))

# --- ENCONTRAR MINIMOS DE F_R ---
cat("\n--- MINIMOS DEL PAISAJE F_R ---\n")

# Buscar minimos locales
find_local_minima <- function(y) {
  n <- length(y)
  minima <- c()
  for (i in 2:(n-1)) {
    if (y[i] < y[i-1] && y[i] < y[i+1]) {
      minima <- c(minima, i)
    }
  }
  minima
}

minima_idx <- find_local_minima(F_R)
cat("Minimos locales encontrados:", length(minima_idx), "\n")
for (i in minima_idx) {
  cat("  x =", round(mids[i], 3), ", F_R =", round(F_R[i], 4), "\n")
}

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("El paisaje de costo F_R:\n")
cat("  - Se construye a partir de la distribucion condicionada P(x|R)\n")
cat("  - Sus MINIMOS indican estados TIPICOS bajo R\n")
cat("  - NO es igual al potencial fisico V(x)\n")
cat("  - La condicion R puede 'invertir' o 'deformar' el paisaje\n")
cat("  - Estados de alta energia fisica pueden tener bajo costo bajo R\n")
