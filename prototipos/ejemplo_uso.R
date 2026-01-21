# =============================================================================
# Ejemplo de uso del Framework UE modular
# =============================================================================
# Muestra como usar las funciones de ue_funciones.R
# =============================================================================

# Cargar funciones
source("ue_funciones.R")

# =============================================================================
# EJEMPLO 1: Pipeline completo en una linea
# =============================================================================

cat("\n")
cat("=============================================\n")
cat("  EJEMPLO 1: Pipeline completo\n")
cat("=============================================\n")

# Definir sistema
grad_V <- function(x) 4 * x * (x^2 - 1) - 0.3  # doble pozo asimetrico
V <- function(x) (x^2 - 1)^2 - 0.3 * x
observable <- function(x) abs(x)

# Ejecutar todo el pipeline
resultado <- ue_pipeline(
  grad_V = grad_V,
  observable = observable,
  theta_star = 0.5,
  tau_star = 25,
  n_traj = 1500,
  n_pasos = 400,
  sigma = 0.4,
  seed = 42,
  validar = FALSE
)

# Visualizar
ue_plot(resultado, potencial = V)

# =============================================================================
# EJEMPLO 2: Paso a paso (mas control)
# =============================================================================

cat("\n")
cat("=============================================\n")
cat("  EJEMPLO 2: Paso a paso\n")
cat("=============================================\n")

# 1. Crear dinamica
cat("\n[Paso 1] Crear dinamica\n")
step_fn <- ue_crear_dinamica(grad_V, dt = 0.01, sigma = 0.4)

# 2. Simular
cat("[Paso 2] Simular trayectorias\n")
trayectorias <- ue_simular_ensemble(
  n_traj = 1000,
  n_pasos = 300,
  step_fn = step_fn,
  seed = 123
)
cat("  Generadas:", length(trayectorias), "trayectorias\n")

# 3. Filtrar por ERS
cat("[Paso 3] Filtrar por condicion R\n")
filtrado <- ue_filtrar_ERS(
  trayectorias,
  observable = observable,
  theta_star = 0.5,
  tau_star = 20
)
cat("  Condicionadas:", filtrado$n_condicionadas, "/", filtrado$n_total, "\n")

# 4. Distribucion condicionada
cat("[Paso 4] Calcular distribucion P(x|R)\n")
distribucion <- ue_distribucion_condicionada(filtrado$trayectorias)
cat("  Media:", round(distribucion$media, 3), "\n")
cat("  Varianza:", round(distribucion$varianza, 3), "\n")

# 5. Paisaje de costo
cat("[Paso 5] Construir paisaje F_R\n")
paisaje <- ue_paisaje_costo(distribucion)
cat("  Epsilon:", round(paisaje$epsilon, 4), "\n")

# 6. Clases metaestables
cat("[Paso 6] Detectar clases\n")
clases <- ue_detectar_clases(paisaje, filtrado$trayectorias)
cat("  Clases:", clases$n_clases, "\n")
cat("  Centros:", round(clases$centros, 2), "\n")
cat("  Pesos:", round(clases$pesos, 3), "\n")

# 7. Metricas
cat("[Paso 7] Calcular metricas\n")
metricas <- ue_metricas_dominancia(clases)
cat("  D_H:", round(metricas$D_H, 2), "\n")
cat("  Delta:", round(metricas$Delta, 3), "\n")
cat("  Interpretacion:", metricas$interpretacion, "\n")

# =============================================================================
# EJEMPLO 3: Comparar sistemas
# =============================================================================

cat("\n")
cat("=============================================\n")
cat("  EJEMPLO 3: Comparar diferentes sistemas\n")
cat("=============================================\n")

# Sistema 1: Doble pozo simetrico
grad_V1 <- function(x) 4 * x * (x^2 - 1)

# Sistema 2: Doble pozo asimetrico
grad_V2 <- function(x) 4 * x * (x^2 - 1) - 0.5

# Sistema 3: Triple pozo
grad_V3 <- function(x) 2 * x * (x^2 - 2) + 3 * x * exp(-3 * x^2)

sistemas <- list(
  list(nombre = "Simetrico", grad = grad_V1),
  list(nombre = "Asimetrico", grad = grad_V2),
  list(nombre = "Triple pozo", grad = grad_V3)
)

cat("\nComparando sistemas con la misma condicion R:\n")
cat("  θ* = 0.5, τ* = 25\n\n")

for (s in sistemas) {
  res <- ue_pipeline(
    grad_V = s$grad,
    observable = observable,
    theta_star = 0.5,
    tau_star = 25,
    n_traj = 1000,
    n_pasos = 300,
    verbose = FALSE
  )

  cat(sprintf("%-15s: %d clases, D_H = %.2f (%s)\n",
              s$nombre,
              res$clases$n_clases,
              res$metricas$D_H,
              res$metricas$interpretacion))
}

# =============================================================================
# EJEMPLO 4: Con validacion
# =============================================================================

cat("\n")
cat("=============================================\n")
cat("  EJEMPLO 4: Pipeline con validacion\n")
cat("=============================================\n")

resultado_validado <- ue_pipeline(
  grad_V = grad_V,
  observable = observable,
  theta_star = 0.5,
  tau_star = 25,
  n_traj = 1500,
  n_pasos = 400,
  validar = TRUE,
  verbose = TRUE
)

# =============================================================================
# EJEMPLO 5: Variar parametros ERS
# =============================================================================

cat("\n")
cat("=============================================\n")
cat("  EJEMPLO 5: Efecto de variar theta*\n")
cat("=============================================\n")

thetas <- c(0.3, 0.5, 0.7, 0.9)
cat("\nFijando tau* = 25, variando theta*:\n\n")

for (th in thetas) {
  res <- ue_pipeline(
    grad_V = grad_V,
    observable = observable,
    theta_star = th,
    tau_star = 25,
    n_traj = 1000,
    verbose = FALSE
  )

  if (!is.null(res$metricas)) {
    cat(sprintf("θ* = %.1f: %d cond (%.1f%%), D_H = %.2f\n",
                th,
                res$filtrado$n_condicionadas,
                res$filtrado$tasa * 100,
                res$metricas$D_H))
  } else {
    cat(sprintf("θ* = %.1f: muy pocas trayectorias\n", th))
  }
}

cat("\n=== Ejemplos completados ===\n")
