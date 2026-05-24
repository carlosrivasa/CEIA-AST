# Solución del primer trabajo práctico de la materia de Análisis de Series de Tiempo I, CEIA - LSE - UBA, 2026
### Carlos Rivas - a2227

---

## Solución punto 1:

Sabemos que:
$$Y_t = R \cos(2\pi(ft + \Phi))$$
Donde:
- $R$ y $\Phi$ son variables aleatorias independientes y $f$ es una frecuencia fija.
- La fase $\Phi$ se distribuye de manera uniforme sobre el intervalo $(0,1)$
- La amplitud $R$ tiene una distribución de Rayleigh con densidad:
$$f_R(r) = r e^{-r^2/2}, \quad \text{para } r > 0$$

---

#### Densidad conjunta:
Como $R$ y $\Phi$ son variables aleatorias independientes, su función de densidad de probabilidad conjunta $f_{R,\Phi}(r, \phi)$ es el producto de sus densidades marginales:
$$f_\Phi(\phi) = 1, \quad \text{ } 0 < \phi < 1$$
$$f_R(r) = r e^{-r^2/2}, \quad \text{ } r > 0$$

De ahí que:
$$f_{R,\Phi}(r, \phi) = r e^{-r^2/2}$$

#### Transformaciones:
Se usa la pista:
$$X_t = R \operatorname{sen}(2\pi(ft + \Phi))$$
$$Y_t = R \cos(2\pi(ft + \Phi))$$

y se expresa $(R, \Phi)$ en términos de $(X, Y)$:
- **Para $R$**:
   $$X^2 + Y^2 = R^2 \left(\operatorname{sen}^2(2\pi(ft + \Phi)) + \cos^2(2\pi(ft + \Phi))\right) = R^2 \implies R = \sqrt{X^2 + Y^2}$$
- **Para $\Phi$**: Si se divide $X$ sobre $Y$ se cancela R:
   $$\frac{X}{Y} = \tan(2\pi(ft + \Phi)) \implies \Phi = \frac{1}{2\pi}\arctan\left(\frac{X}{Y}\right) - ft$$

#### Jacobiano:
Se calcula el determinante del Jacobiano, con las derivadas parciales de $(X, Y)$ respecto a $(R, \Phi)$:
$$J = \begin{vmatrix} \frac{\partial X}{\partial r} & \frac{\partial X}{\partial \phi} \\ \frac{\partial Y}{\partial r} & \frac{\partial Y}{\partial \phi} \end{vmatrix} = \begin{vmatrix} \operatorname{sen}(\theta) & 2\pi r \cos(\theta) \\ \cos(\theta) & -2\pi r \operatorname{sen}(\theta) \end{vmatrix}$$
Con $\theta = 2\pi(ft + \phi)$. Evaluando el determinante:
$$|J| = -2\pi r \operatorname{sen}^2(\theta) - 2\pi r \cos^2(\theta) = -2\pi r$$

El valor absoluto del determinante del Jacobiano es $||J|| = 2\pi r$.

#### Densidad conjunta de $(X, Y)$
Se aplica el cambio de variables por la transformación:
$$f_{X,Y}(x,y) = f_{R,\Phi}(r, \phi) \cdot \frac{1}{||J||} = \frac{r e^{-r^2/2}}{2\pi r} = \frac{1}{2\pi} e^{-\frac{x^2+y^2}{2}}$$

Al separar los términos:
$$f_{X,Y}(x,y) = \left( \frac{1}{\sqrt{2\pi}} e^{-x^2/2} \right) \cdot \left( \frac{1}{\sqrt{2\pi}} e^{-y^2/2} \right)$$

Y se observa que cada término corresponde a una función de densidad de una distribución normal estándar $N(0,1)$.


---

## Solución punto 2:

Se tiene el proceso:
$$Y_t = \mu + \varepsilon_t + \varepsilon_{t-1}$$

Donde $\varepsilon_t$ es un ruido blanco con media $0$. 

La media muestral es $\bar{Y} = \frac{1}{n} \sum_{t=1}^{n} Y_t$. 
Se pide encontrar $\operatorname{Var}(\bar{Y})$.

---

La varianza de la media muestral en un proceso de medias móviles de orden 1 MA(1).
Para calcular la varianza de una suma, necesitamos saber cómo covarían los términos entre sí. 

Se asume que $\operatorname{Var}(\varepsilon_t) = \sigma_e^2$ y se calcula $\gamma_k = \operatorname{Cov}(Y_t, Y_{t-k})$ para los diferentes rezagos $k$:
* **Para $k=0$ (Varianza):**
  $$\gamma_0 = \operatorname{Var}(Y_t) = \operatorname{Var}(\mu + \varepsilon_t + \varepsilon_{t-1})$$
  $$\gamma_0 = \operatorname{Var}(\varepsilon_t) + \operatorname{Var}(\varepsilon_{t-1}) = \sigma_e^2 + \sigma_e^2 = 2\sigma_e^2$$

* **Para $k=1$ (Autocovarianza de primer orden):**
  $$\gamma_1 = \operatorname{Cov}(Y_t, Y_{t-1}) = \operatorname{Cov}(\mu + \varepsilon_t + \varepsilon_{t-1}, \mu + \varepsilon_{t-1} + \varepsilon_{t-2})$$
  $$\gamma_1 = \operatorname{Cov}(\varepsilon_t + \varepsilon_{t-1}, \varepsilon_{t-1} + \varepsilon_{t-2})$$
  $$\gamma_1 = \operatorname{Cov}(\varepsilon_t, \varepsilon_{t-1}) + \operatorname{Cov}(\varepsilon_t, \varepsilon_{t-2}) + \operatorname{Cov}(\varepsilon_{t-1}, \varepsilon_{t-1}) + \operatorname{Cov}(\varepsilon_{t-1}, \varepsilon_{t-2})$$
  $$\gamma_1 = 0 + 0 + \sigma_e^2 + 0 = \sigma_e^2$$

* **Para $k \ge 2$:**
  $$\gamma_k = \operatorname{Cov}(Y_t, Y_{t-k}) = \operatorname{Cov}(\varepsilon_t + \varepsilon_{t-1}, \varepsilon_{t-k} + \varepsilon_{t-k-1}) = 0$$

Varianza de la suma:

$$\operatorname{Var}(\bar{Y}) = \operatorname{Var}\left( \frac{1}{n} \sum_{t=1}^{n} Y_t \right) = \frac{1}{n^2} \operatorname{Var}\left( \sum_{t=1}^{n} Y_t \right)$$

Por las propiedades de la varianza aplicadas a una suma de variables dependientes:
$$\operatorname{Var}\left( \sum_{t=1}^{n} Y_t \right) = \sum_{t=1}^{n} \operatorname{Var}(Y_t) + \sum_{t \neq s} \operatorname{Cov}(Y_t, Y_s)$$

Como $\operatorname{Cov}(Y_t, Y_s) = \operatorname{Cov}(Y_s, Y_t)$:
$$\operatorname{Var}\left( \sum_{t=1}^{n} Y_t \right) = \sum_{t=1}^{n} \operatorname{Var}(Y_t) + 2 \sum_{t=1}^{n-1} \sum_{s=t+1}^{n} \operatorname{Cov}(Y_t, Y_s)$$

Al tratarse de un proceso débilmente estacionario (media y varianza constantes, varianza finita y la autocovarianza depende de k y no de t):
1. **La varianza es constante:** $\operatorname{Var}(Y_t) = \gamma_0$.
2. **Las covarianzas dependen del rezago:** Se define el rezago como $k = |t - s|$. En una muestra de tamaño $n$, existen exactamente $(n-k)$ pares de observaciones que están separados por una distancia $k$. Como la covarianza para una distancia $k$ es $\gamma_k$, la suma de todas las covarianzas para ese rezago específico es $(n-k)\gamma_k$.

Al sustituir el término de la varianza pura y el término de las covarianzas agrupadas por rezago:
$$\operatorname{Var}(\bar{Y}) = \frac{1}{n^2} \left[ n\gamma_0 + 2\sum_{k=1}^{n-1} (n-k)\gamma_k \right]$$



Como $\gamma_k = 0$ para todo $k \ge 2$, la expresión se reduce a:
$$\operatorname{Var}(\bar{Y}) = \frac{1}{n^2} [n\gamma_0 + 2(n-1)\gamma_1]$$

Al sustituir los valores de $\gamma_0 = 2\sigma_e^2$ y $\gamma_1 = \sigma_e^2$:
$$\operatorname{Var}(\bar{Y}) = \frac{1}{n^2} [n(2\sigma_e^2) + 2(n-1)\sigma_e^2]$$
$$\operatorname{Var}(\bar{Y}) = \frac{\sigma_e^2}{n^2} [2n + 2n - 2] = \frac{\sigma_e^2(4n - 2)}{n^2}$$

---

La varianza de la media muestral en un proceso de medias móviles de orden 1 MA(1).
Para calcular la varianza de una suma, necesitamos saber cómo covarían los términos entre sí:

Se calcula $\gamma_k = \operatorname{Cov}(Y_t, Y_{t-k})$ para los diferentes rezagos $k$:
* **Para $k=0$ (Varianza):**
  $$\gamma_0 = \operatorname{Var}(Y_t) = \operatorname{Var}(\mu + \varepsilon_t + \varepsilon_{t-1}) = \operatorname{Var}(\varepsilon_t) + \operatorname{Var}(\varepsilon_{t-1}) = 2\sigma^2$$
* **Para $k=1$ (Autocovarianza de primer orden):**
  $$\gamma_1 = \operatorname{Cov}(\mu + \varepsilon_t + \varepsilon_{t-1}, \mu + \varepsilon_{t-1} + \varepsilon_{t-2}) = \operatorname{Cov}(\varepsilon_{t-1}, \varepsilon_{t-1}) = \sigma^2$$
* **Para $k \ge 2$:**
  $$\gamma_k = 0$$

#### Paso 2: Desarrollar la varianza de la suma
Por propiedades de los operadores de varianza para variables dependientes:
$$\operatorname{Var}(\bar{Y}) = \operatorname{Var}\left( \frac{1}{n} \sum_{t=1}^{n} Y_t \right) = \frac{1}{n^2} \left[ n\gamma_0 + 2\sum_{k=1}^{n-1} (n-k)\gamma_k \right]$$

Dado que $\gamma_k = 0$ para todo $k \ge 2$, la expresión se reduce a:
$$\operatorname{Var}(\bar{Y}) = \frac{1}{n^2} [n\gamma_0 + 2(n-1)\gamma_1]$$

#### Paso 3: Sustitución y simplificación
Sustituyendo los valores de $\gamma_0 = 2\sigma^2$ y $\gamma_1 = \sigma^2$:
$$\operatorname{Var}(\bar{Y}) = \frac{1}{n^2} [n(2\sigma^2) + 2(n-1)\sigma^2]$$
$$\operatorname{Var}(\bar{Y}) = \frac{\sigma^2}{n^2} [2n + 2n - 2] = \frac{\sigma^2(4n - 2)}{n^2}$$

#### Conclusión
La varianza de la media muestral es:
$$\operatorname{Var}(\bar{Y}) = \frac{2\sigma^2(2n - 1)}{n^2} = \frac{4\sigma^2}{n} - \frac{2\sigma^2}{n^2}$$

---

## Solución punto 3: 

### Enunciado
Se tiene el proceso:
$$Y_t = 4 + 6\varepsilon_t + 1\varepsilon_{t-1} + 8\varepsilon_{t-2}$$

Y se asume que $\varepsilon_t$ es un ruido blanco con media $0$ y varianza $\sigma^2$.

---

#### Autocovarianzas $\gamma_k$
El proceso es un $\text{MA}(2)$ con coeficientes $\theta_0 = 6$, $\theta_1 = 1$, y $\theta_2 = 8$.
* **Rezagos $k=0$ (Varianza):**
  $$\gamma_0 = \operatorname{Var}(Y_t) = (6^2 + 1^2 + 8^2)\sigma^2 = (36 + 1 + 64)\sigma^2 = 101\sigma^2$$
* **Rezago $k=1$:**
  $$\gamma_1 = \operatorname{Cov}(Y_t, Y_{t-1}) = (\theta_1\theta_0 + \theta_2\theta_1)\sigma^2 = (1 \cdot 6 + 8 \cdot 1)\sigma^2 = 14\sigma^2$$
* **Rezago $k=2$:**
  $$\gamma_2 = \operatorname{Cov}(Y_t, Y_{t-2}) = (\theta_2\theta_0)\sigma^2 = (8 \cdot 6)\sigma^2 = 48\sigma^2$$
* **Rezagos $k \ge 3$:**
  $$\gamma_k = 0$$

Autocorrelaciones $\rho_k = \frac{\gamma_k}{\gamma_0}$
* **$\rho_0$:**
  $$\rho_0 = \frac{101\sigma^2}{101\sigma^2} = 1$$
* **$\rho_1$:**
  $$\rho_1 = \frac{14\sigma^2}{101\sigma^2} = \frac{14}{101}$$
* **$\rho_2$:**
  $$\rho_2 = \frac{48\sigma^2}{101\sigma^2} = \frac{48}{101}$$
* **$\rho_k$ para $k \ge 3$:**
  $$\rho_k = 0$$

Función de autocorrelación:
$$\rho_k = \begin{cases} 1 & \text{si } k = 0 \\ \frac{14}{101} & \text{si } k = 1 \\ \frac{48}{101} & \text{si } k = 2 \\ 0 & \text{si } k \ge 3 \end{cases}$$