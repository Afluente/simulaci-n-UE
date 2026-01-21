# =============================================================================
# CONCEPTO 2: Probabilidad Condicionada P(ψ|R)
# =============================================================================
#
# Definicion:
#   p_n(ψ|R) = P_n(ψ ∩ R) / P_n(R)
#
# Idea clave:
#   La distribucion de estados CAMBIA cuando condicionamos a R.
#   Estados que eran raros pueden volverse tipicos bajo R.
#
# Ejemplo:
#   Sin condicionar: el sistema pasa igual tiempo en ambos pozos
#   Condicionando a "energia alta": prefiere la barrera central
# =============================================================================

cat("=== PRUEBA: Probabilidad Condicionada P(ψ|R) ===\n\n")

# --- PARAMETROS ---
set.seed(42)
n_traj <- 2000
traj_length <- 500
dt <- 0.01
noise <- 0.5

# --- POTENCIAL DOBLE POZO ---
V <- function(x) (x^2 - 1)^2
grad_V <- function(x) 4 * x * (x^2 - 1)

# --- DINAMICA LANGEVIN ---
step <- function(x, dt, sigma) {
  x - grad_V(x) * dt + sigma * sqrt(dt) * rnorm(1)
}

simulate <- function(x0, n, dt, sigma) {
  x <- numeric(n)
  x[1] <- x0
  for (t in 2:n) {
    x[t] <- step(x[t-1], dt, sigma)
  }
  x
}

# --- OBSERVABLE Y CONDICION R ---
# Observable: cercania al origen (barrera)
observable <- function(x) -abs(x)  # mas alto = mas cerca del origen

theta_threshold <- -0.3  # cerca del origen
tau_persistence <- 30

cumple_ERS <- function(traj, theta_star, tau_star) {
  theta_t <- observable(traj)
  condicion <- theta_t >= theta_star  # |x| <= 0.3
  if (!any(condicion)) return(FALSE)
  rachas <- rle(condicion)
  duraciones <- rachas$lengths[rachas$values == TRUE]
  if (length(duraciones) == 0) return(FALSE)
  max(duraciones) >= tau_star
}

# --- SIMULAR ---
cat("Simulando", n_traj, "trayectorias en doble pozo...\n")
cat("Condicion R: |x| <= 0.3 durante >=", tau_persistence, "pasos\n")
cat("(permanecer cerca de la barrera central)\n\n")

trayectorias <- vector("list", n_traj)
cumple_R <- logical(n_traj)

for (i in 1:n_traj) {
  x0 <- runif(1, -2, 2)
  trayectorias[[i]] <- simulate(x0, traj_length, dt, noise)
  cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_threshold, tau_persistence)
}

n_cond <- sum(cumple_R)
cat("Trayectorias que cumplen R:", n_cond, "(", round(100*n_cond/n_traj, 1), "%)\n\n")

# --- CALCULAR DISTRIBUCIONES ---
# P(x) - distribucion sin condicionar
x_todas <- unlist(trayectorias)

# P(x|R) - distribucion condicionada
x_condicionadas <- unlist(trayectorias[cumple_R])

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Potencial
x_grid <- seq(-2, 2, length.out = 200)
plot(x_grid, V(x_grid), type = "l", lwd = 2,
     main = "Potencial V(x) = (x²-1)²",
     xlab = "x", ylab = "V(x)")
abline(v = c(-1, 1), col = "blue", lty = 2)
text(c(-1, 1), c(0.1, 0.1), c("pozo", "pozo"), col = "blue")

# 2. P(x) sin condicionar
hist(x_todas, breaks = 80, prob = TRUE, col = "gray80", border = "white",
     main = "P(x) - Sin condicionar",
     xlab = "x", xlim = c(-2, 2))
abline(v = c(-1, 1), col = "blue", lty = 2)
text(0, par("usr")[4]*0.9, "barrera", col = "red", cex = 0.8)

# 3. P(x|R) condicionada
if (n_cond > 0) {
  hist(x_condicionadas, breaks = 80, prob = TRUE, col = "steelblue", border = "white",
       main = "P(x|R) - Condicionada a R",
       xlab = "x", xlim = c(-2, 2))
  abline(v = c(-1, 1), col = "blue", lty = 2)
  abline(v = c(-0.3, 0.3), col = "red", lty = 2)
  text(0, par("usr")[4]*0.9, "zona R", col = "red", cex = 0.8)
}

# 4. Comparacion superpuesta
hist(x_todas, breaks = 80, prob = TRUE, col = rgb(0.5, 0.5, 0.5, 0.5),
     border = "white", main = "Comparacion P(x) vs P(x|R)",
     xlab = "x", xlim = c(-2, 2))
if (n_cond > 0) {
  hist(x_condicionadas, breaks = 80, prob = TRUE,
       col = rgb(0.2, 0.4, 0.8, 0.5), border = "white", add = TRUE)
}
legend("topright", c("P(x)", "P(x|R)"),
       fill = c(rgb(0.5,0.5,0.5,0.5), rgb(0.2,0.4,0.8,0.5)))

par(mfrow = c(1, 1))

# --- ESTADISTICAS ---
cat("\n--- ESTADISTICAS ---\n")
cat("Media P(x):", round(mean(x_todas), 4), "\n")
cat("Media P(x|R):", round(mean(x_condicionadas), 4), "\n")
cat("Varianza P(x):", round(var(x_todas), 4), "\n")
cat("Varianza P(x|R):", round(var(x_condicionadas), 4), "\n")

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("Al condicionar a R (permanecer cerca de la barrera):\n")
cat("  - P(x) tiene picos en los pozos (x = ±1)\n")
cat("  - P(x|R) concentra masa cerca del origen (x ≈ 0)\n")
cat("  - La barrera, antes improbable, se vuelve TIPICA bajo R\n")
cat("  - Esto es el nucleo del framework UE: cambiar la pregunta\n")
