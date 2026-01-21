# =============================================================================
# CONCEPTO 5: Metricas de Dominancia
# =============================================================================
#
# Una vez identificadas las clases A_k, medimos si alguna DOMINA:
#
# Pesos de clase:
#   p_k = P(A_k | R) = ∫_{A_k} p(ψ|R) dψ
#
# Metricas:
#   - Ratio de dominancia: D_H = p_1 / p_2 (clase mas probable / segunda)
#   - Gap: Δ = p_1 - p_2
#   - Entropia efectiva: S_eff = -Σ p_k log(p_k)
#
# Umbrales:
#   - D_H > 10 o Δ > 0.5 → Dominancia fuerte
#   - D_H ∈ [3, 10] → Dominancia moderada
#   - D_H ≈ 1 → Equilibrio (sin dominancia)
#
# =============================================================================

cat("=== PRUEBA: Metricas de Dominancia ===\n\n")

# --- PARAMETROS ---
set.seed(42)
n_traj <- 2500
traj_length <- 500
dt <- 0.01
noise <- 0.4

# --- POTENCIAL ASIMETRICO (favorece un pozo) ---
# Parametro de asimetria
asimetria <- 0.4  # Cambiar para ver efecto en dominancia

V <- function(x) (x^2 - 1)^2 - asimetria * x
grad_V <- function(x) 4 * x * (x^2 - 1) - asimetria

cat("Potencial: V(x) = (x²-1)² -", asimetria, "·x\n")
cat("(asimetria =", asimetria, "favorece pozo derecho)\n\n")

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

# --- CONDICION R ---
observable <- function(x) abs(x)  # distancia al origen
theta_threshold <- 0.5
tau_persistence <- 30

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
cat("Simulando", n_traj, "trayectorias...\n")

trayectorias <- vector("list", n_traj)
cumple_R <- logical(n_traj)

for (i in 1:n_traj) {
  x0 <- runif(1, -2, 2)
  trayectorias[[i]] <- simulate(x0, traj_length, dt, noise)
  cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_threshold, tau_persistence)
}

n_cond <- sum(cumple_R)
cat("Condicionadas:", n_cond, "\n\n")

# --- DEFINIR CLASES SIMPLES ---
# A_1: pozo izquierdo (x < 0)
# A_2: pozo derecho (x > 0)

x_cond <- unlist(trayectorias[cumple_R])

n_izq <- sum(x_cond < 0)
n_der <- sum(x_cond > 0)
n_total <- length(x_cond)

# --- CALCULAR PESOS ---
p_1 <- n_izq / n_total  # P(A_1 | R)
p_2 <- n_der / n_total  # P(A_2 | R)

# Ordenar de mayor a menor
if (p_2 > p_1) {
  temp <- p_1; p_1 <- p_2; p_2 <- temp
  clase_dominante <- "derecha (A_2)"
} else {
  clase_dominante <- "izquierda (A_1)"
}

# --- METRICAS DE DOMINANCIA ---
cat("=== PESOS DE CLASE ===\n")
cat("  P(A_1|R) - izquierda:", round(n_izq/n_total, 4), "\n")
cat("  P(A_2|R) - derecha:  ", round(n_der/n_total, 4), "\n\n")

# Ratio D_H
D_H <- p_1 / p_2
cat("=== METRICAS ===\n")
cat("  Ratio D_H = p_max / p_second =", round(D_H, 3), "\n")

# Gap
Delta <- p_1 - p_2
cat("  Gap Δ = p_max - p_second =", round(Delta, 4), "\n")

# Entropia efectiva
p_vec <- c(p_1, p_2)
S_eff <- -sum(p_vec * log(p_vec))
S_max <- log(2)  # entropia maxima con 2 clases
cat("  Entropia S_eff =", round(S_eff, 4), "(max =", round(S_max, 4), ")\n")
cat("  Entropia normalizada =", round(S_eff/S_max, 4), "\n\n")

# --- INTERPRETACION ---
cat("=== INTERPRETACION ===\n")

if (D_H > 10) {
  cat("DOMINANCIA FUERTE (D_H > 10)\n")
  cat("  La clase", clase_dominante, "domina claramente\n")
} else if (D_H > 3) {
  cat("DOMINANCIA MODERADA (3 < D_H < 10)\n")
  cat("  La clase", clase_dominante, "es preferida pero no absoluta\n")
} else if (D_H > 1.5) {
  cat("DOMINANCIA DEBIL (1.5 < D_H < 3)\n")
  cat("  Ligera preferencia por clase", clase_dominante, "\n")
} else {
  cat("EQUILIBRIO (D_H ≈ 1)\n")
  cat("  Las clases tienen probabilidad similar bajo R\n")
}

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Potencial
x_grid <- seq(-2.5, 2.5, length.out = 200)
plot(x_grid, V(x_grid), type = "l", lwd = 2, col = "darkred",
     main = paste0("Potencial (asimetria = ", asimetria, ")"),
     xlab = "x", ylab = "V(x)")
abline(v = 0, col = "gray", lty = 2)

# 2. Distribucion con clases
hist(x_cond, breaks = 60, col = "gray80", border = "white",
     main = "Distribucion P(x|R)",
     xlab = "x", xlim = c(-2.5, 2.5))
hist(x_cond[x_cond < 0], breaks = 30, col = adjustcolor("steelblue", 0.6),
     border = "white", add = TRUE)
hist(x_cond[x_cond > 0], breaks = 30, col = adjustcolor("coral", 0.6),
     border = "white", add = TRUE)
abline(v = 0, col = "black", lwd = 2)
legend("topright", c("A_1 (izq)", "A_2 (der)"),
       fill = c(adjustcolor("steelblue", 0.6), adjustcolor("coral", 0.6)))

# 3. Barras de pesos
barplot(c(n_izq/n_total, n_der/n_total),
        names.arg = c("A_1 (izq)", "A_2 (der)"),
        col = c("steelblue", "coral"),
        main = "Pesos P(A_k|R)",
        ylab = "probabilidad", ylim = c(0, 1))
abline(h = 0.5, col = "gray", lty = 2)

# 4. Metricas visuales
plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 4),
     xlab = "", ylab = "", axes = FALSE,
     main = "Resumen de Metricas")

text(5, 3.5, paste("D_H =", round(D_H, 2)), cex = 1.5, font = 2)
text(5, 2.5, paste("Δ =", round(Delta, 3)), cex = 1.2)
text(5, 1.5, paste("S_eff =", round(S_eff, 3)), cex = 1.2)

# Barra visual de dominancia
rect(1, 0.3, 1 + 8 * (D_H / max(D_H, 10)), 0.7,
     col = ifelse(D_H > 3, "coral", "steelblue"), border = NA)
rect(1, 0.3, 9, 0.7, border = "black")
text(5, 0.1, "D_H: 1 -------- 10+", cex = 0.8)

par(mfrow = c(1, 1))

# --- EXPERIMENTO: VARIAR ASIMETRIA ---
cat("\n=== EXPERIMENTO: Efecto de la asimetria ===\n\n")

asimetrias <- c(0, 0.1, 0.2, 0.4, 0.6)
resultados <- data.frame(asimetria = numeric(),
                         D_H = numeric(),
                         Delta = numeric(),
                         S_eff = numeric())

for (a in asimetrias) {
  V_a <- function(x) (x^2 - 1)^2 - a * x
  grad_V_a <- function(x) 4 * x * (x^2 - 1) - a

  step_a <- function(x, dt, sigma) {
    x - grad_V_a(x) * dt + sigma * sqrt(dt) * rnorm(1)
  }

  simulate_a <- function(x0, n, dt, sigma) {
    x <- numeric(n)
    x[1] <- x0
    for (t in 2:n) {
      x[t] <- step_a(x[t-1], dt, sigma)
      x[t] <- max(-3, min(3, x[t]))
    }
    x
  }

  # Simulacion rapida
  set.seed(42)
  x_all <- c()
  for (i in 1:500) {
    traj <- simulate_a(runif(1, -2, 2), 300, dt, noise)
    if (cumple_ERS(traj, theta_threshold, tau_persistence)) {
      x_all <- c(x_all, traj)
    }
  }

  if (length(x_all) > 100) {
    p1 <- max(mean(x_all < 0), mean(x_all > 0))
    p2 <- min(mean(x_all < 0), mean(x_all > 0))
    D_H_a <- p1 / p2
    Delta_a <- p1 - p2
    S_a <- -p1*log(p1) - p2*log(p2)
  } else {
    D_H_a <- NA
    Delta_a <- NA
    S_a <- NA
  }

  resultados <- rbind(resultados,
                      data.frame(asimetria = a, D_H = D_H_a,
                                 Delta = Delta_a, S_eff = S_a))
}

print(resultados)

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("Las metricas de dominancia cuantifican:\n")
cat("  - D_H: cuanto mas probable es la clase dominante vs la segunda\n")
cat("  - Δ: diferencia absoluta entre las dos principales\n")
cat("  - S_eff: 'uniformidad' de la distribucion entre clases\n")
cat("\nA mayor asimetria del potencial, mayor dominancia bajo R.\n")
cat("Pero la relacion no es trivial: depende de como R interactua con V(x).\n")
