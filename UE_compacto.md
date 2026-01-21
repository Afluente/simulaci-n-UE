# Universo Emergente — Guía Compacta para Simulación

## 1. Conceptos Fundamentales

### Selección Condicionada (SC)
- **Pregunta clave**: "¿Qué es típico dado que R se mantiene?" (no "¿qué es típico en general?")
- Cambia la medida de referencia de P a P(·|R)

### Evento Raro Sostenido (ERS / R)
Definición ex-ante con dos componentes:
```
R = { θ(t) ≥ θ* durante τ ≥ τ* }
```
- `θ(t)`: observable medible
- `θ*`: umbral (threshold)
- `τ*`: tiempo mínimo de persistencia

## 2. Núcleo Formal

### Probabilidad Condicionada
```
p_n(ψ|R) = P_n(ψ ∩ R) / P_n(R)
```

### Paisaje Condicionado (Función de Costo)
```
F_n^R[ψ] = -ε_n · log p_n(ψ|R)
```
donde `ε_n → 0` es parámetro de escala (típicamente `1/n` o `1/log n`)

### Principio de Grandes Desviaciones Condicionado
```
p_n(ψ|R) ≈ exp(-F_n^R[ψ] / ε_n)
```

## 3. Definición de Clases (A_k)

### Capa 1: Geométrica (Cuencas)
```
A_k = { ψ : flujo_gradiente(ψ) → mínimo_k de F^R }
```

### Capa 2: Dinámica (Metaestabilidad)
Criterio obligatorio:
```
τ_relax(A_k) << E[τ_exit(A_k)]
```
- `τ_relax`: tiempo de relajación interna
- `τ_exit`: tiempo de escape de la cuenca

## 4. Métricas de Dominancia

### Pesos de Clase
```
p_k = P(A_k | R) = ∫_{A_k} p_n(ψ|R) dψ
```

### Indicadores de Dominancia
- **Ratio**: `D_H = p_1 / p_2`
- **Gap**: `Δ = p_1 - p_2`
- **Entropía efectiva**: `S_eff = -Σ p_k log p_k`

### Umbrales Sugeridos
- Dominancia fuerte: `D_H > 10` o `Δ > 0.5`
- Dominancia moderada: `D_H ∈ [3, 10]`

## 5. Protocolo de Simulación

### Paso 1: Pre-registro
Fijar ANTES de simular:
- Observable `θ(t)` y su fórmula
- Umbral `θ*` y persistencia `τ*`
- Espacio de configuraciones Ω
- Dinámica/reglas de evolución

### Paso 2: Generación de Trayectorias
```python
# Pseudocódigo
trayectorias = []
while len(trayectorias) < N_objetivo:
    traj = simular_dinamica()
    if cumple_R(traj, theta_star, tau_star):
        trayectorias.append(traj)
```

### Paso 3: Construcción del Paisaje
Para cada configuración ψ visitada:
```python
frecuencia[ψ] = cuenta[ψ] / total_pasos
F_R[ψ] = -epsilon * log(frecuencia[ψ])
```

### Paso 4: Identificación de Clases
1. Encontrar mínimos locales de F_R
2. Asignar cada ψ a cuenca por descenso de gradiente
3. Verificar metaestabilidad de cada cuenca

### Paso 5: Cálculo de Dominancia
```python
p_k = sum(frecuencia[ψ] for ψ in A_k)
D_H = p_1 / p_2  # ordenar p_k descendente
```

## 6. Validación Obligatoria

### Test de Robustez
Variar parámetros ±20%:
- `θ* → θ* ± δθ`
- `τ* → τ* ± δτ`
- Verificar: clases dominantes persisten

### Test de Ablación
Romper mecanismos preservando marginales:
- Permutar temporalmente señales
- Aleatorizar conexiones
- Verificar: dominancia desaparece o cambia

### Out-of-Sample
- Dividir datos/trayectorias en entrenamiento y test
- Verificar: predicciones se mantienen

## 7. Multi-Dominio (Extensión)

### Información Condicional
```
I(X;Y|R) = H(X|R) + H(Y|R) - H(X,Y|R)
```

### Dependencia Inter-dominio
```
F_conjunta^R[ψ_1, ψ_2] ≠ F_1^R[ψ_1] + F_2^R[ψ_2]
```
La diferencia mide acoplamiento estructural.

## 8. Parámetros Típicos para Simulación

| Parámetro | Símbolo | Rango típico |
|-----------|---------|--------------|
| Trayectorias condicionadas | N | 10³ - 10⁵ |
| Pasos por trayectoria | T | 10² - 10⁴ |
| Escala de resolución | ε_n | 1/N o 1/log(N) |
| Umbral dominancia | D_H | > 3 (moderada), > 10 (fuerte) |

## 9. Checklist Pre-Simulación

- [ ] Observable θ(t) definido y medible
- [ ] Umbral θ* justificado (no ajustado post-hoc)
- [ ] Persistencia τ* especificada
- [ ] Espacio de estados Ω acotado
- [ ] Dinámica completamente especificada
- [ ] Criterio de metaestabilidad definido
- [ ] Tests de validación planificados

## 10. Uso del Simulador

### Compilación
```bash
cd build
cmake ..
make -j4
```

### Ejecución
```bash
./simulacion_ue                    # Completa con gráficos
./simulacion_ue --no-plot          # Solo texto
./simulacion_ue --validate         # Con protocolo de validación UE
./simulacion_ue --quick            # Modo rápido (menos trayectorias)
./simulacion_ue --no-plot --validate  # Validación sin gráficos
```

### Salida Esperada
- **Trayectorias bajo R**: ~10-15% para el doble pozo
- **Clases detectadas**: 2 (pozos en x=±1)
- **D_H**: ~1.0 (equilibrio simétrico)
- **Validación**: 3 tests (Robustez, Ablación, Out-of-Sample)

## 11. Errores Comunes a Evitar

1. **Posdiction**: Ajustar θ* después de ver resultados
2. **Confundir F con energía física**: F es costo informacional
3. **Ignorar metaestabilidad**: Cuencas sin τ_relax << τ_exit no son clases válidas
4. **Overfitting**: No validar out-of-sample
5. **Ablación superficial**: Debe preservar marginales mientras rompe mecanismo
6. **Ruido excesivo**: Si `noise_sigma` es muy alto, las transiciones son frecuentes y τ_exit bajo
