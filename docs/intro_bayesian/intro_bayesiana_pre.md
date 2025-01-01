---
title: "Introducción a la estadística Bayesiana- Completa tus notas"
layout: home
nav_order: 1
index: true
redirect: false
parent: Temario
---

Este tutorial fue creado por Rosana Zenil-Ferguson (Deciembre 2025)

# ¿Qué es la estadística Bayesiana?

La estadística bayesiana es una rama del estudio de la Estadística que se enfoca en la inferencia y pronóstico de la probabilidades de eventos inciertos. El objetivo más importante de la estadística bayesiana es siempre proporcionar una medida de la incertidumbre y no solamente un estimador. Esta idea la discutiremos una y otra vez a través de los ejercicios del taller. 

Existen dos diferencias centrales de la estadística Bayesiana cuando la comparamos con la forma de hacer estadística tradicional (ya sea frecuentista o bajo el método de verosimilitud). Estas diferencias son las siguientes:


1. Cómo se miden las probabilidades en Bayesiana es distinto a cómo se define la probabilidad en estadística tradicional.

2.  Los parámetros en estadística Bayesiana son desconocidos pero también son **variables aleatorias**. En estadística tradicional, los parámetros son desconocidos pero son variables fijas. Esta diferencia es esencial. 

## Ejemplo: El experimento de un reality show "El amor es ciego"

Las citas en línea o a ciegas son difíciles, especialmente si tu personalidad es reservada. Estas citas son especialmente difíciles para los estadísticos bayesianos porque entienden sus probabilidades muy bien. Imaginemos por un momento que vamos a ver un reality show que se llama "El amor es ciego", como estadísticos bayesianos vamos es estimar la probabilidad de que uno de los participantes se enamore y se case. El experimento consiste en poner a una persona, llamada Silvia, en una serie de citas a ciegas hasta que encuentre al amor de su vida y se case. Para lograrlo Silvia debe tener un indice alto de carisma. Como estadísticos, vamos a estimar la probabilidad de que una persona tenga gran carisma. 

### Definición de eventos

Defininamos dos eventos a los que les vamos a medir la probabilidad

1. $A$= Es el evento que representa que Silvia es muy carismática

2. $B$  = Es el evento que representa cuantas de las citas a ciegas de Silvia salieron bien.

### Definición de probabilidad a priori


### El objetivo final: $P(A|B)$

En realidad lo que nos interesa como estadísticos bayesianos es saber el resultado final para Silvia ¿acabará enamorandose?, en términos de probabilidad nos interesa lo que llamamos la probabilidad a posteriori, es decir  $P(A|B)$ = la probabilidad del carisma de Silvia dadas las citas exitosas durante el reality show. En algunos contextos la probabilidad a posteriori se le conoce como la actualización de las probabilidades porque conocemos más del carisma de Silvia cada vez que tiene una cita en el reality. 

### Elicitación de la distribución a priori

Para definir una probabilidad de manera correcta necesitamos aprender un concepto nuevo   **variables aleatorias**. Las variables aleatorias son una **función matemática**  que nos lleva del espacio de eventos al espacio de los números.


```{r, eval= FALSE}
theta <- .1 # Silvia no es carismática
# Beta tiene dos parametros a y b. La media de una distribución beta es a/(a+b)


prior_distribution<-function(theta, a, b){
beta_val<-dbeta(theta, a, b)
return(beta_val)}# Por qué es más grande que uno?

a     <- 2
b     <- 5
prior_distribution (0.1, a, b)
```

```{r,eval=FALSE}
# Tidyverse es un paquete que crea y oraniza "tibbles" tablas con un formato más sencillo de organizar

# install.packages(tidyverse) #Instala estos paquetes si no los tienes 

# ggplot2 paquete para visualizar las figuras

# install.packages(ggplot2) # Instala estos paquetes si no los tienes.
# install.packages(dplyr) # Instala estos paquetes si no los tienes.

library(tidyverse,quietly=TRUE)
library(ggplot2,quietly=TRUE)
library(dplyr,quietly=TRUE)

length <- 1e4

# Creando una tabla de valores a y b para que podamos ver las diferentes formas que crea la distribución Beta

d <-
  crossing(shape1 = c(.1, 1:4),    
           shape2 = c(.1, 1:4)) %>%
  expand(nesting(shape1, shape2), 
         x = seq(from = 0, to = 1, length.out = length)) %>% 
  mutate(a     = str_c("a =", shape1),
         b     = str_c("b =", shape2),
         group = rep(1:length, each = 25))

head(d) #you can see what is in the table


# Graficando las distribuciones beta
d %>% 
  ggplot(aes(x = x, group = group)) +
  
  geom_line(aes(y = dbeta(x, shape1 = shape1, shape2 = shape2)),
            color = "hotpink", size = 1.25) +
  scale_x_continuous(breaks = c(0, .5, 1)) +coord_cartesian(ylim = c(0,3))+
  labs(x = expression(theta),
       y = expression(paste("P(",expression(theta)))) +
  theme(panel.grid = element_blank()) +
  facet_grid(b~a)
```

Pasemos un momento interpretando las distribuciones beta. Estas son distribuciones a priori potenciales sujetas a nuestras creencias sobre el carisma de Silvia. 

Con estas distribuciones mostramos el segundo punto esencial de la estadística bayesiana:

**2. Los parametros son desconocidos pero inciertos, necesitan una distribución de probabilidad.**

Logramos esto a través de $\theta$ con la **distribución a priori** $Beta(a,b)$.

## Los datos, las citas a ciegas

En el reality show, Silvia tiene tres citas completamente a ciegas  con la misma persona (ninguno se ve solo se escuchan el uno al otro). En nuestros términos estadísticos $N=3$ es el número total de citas con la misma persona. ¿Cómo modelamos los resultados de estas tres citas, sabiendo el carisma de Silvia? Si recordamos el evento $B$ es el número de citas que fueron bien, pero necesitamos llevar el evento B al espacio de los números reales. 

```{r,eval=FALSE}
# install.packages(grid) # Installa esto si no lo tienes.
library(grid)

date_number  <- 0:3
theta<-0.1
df <- data.frame(x = date_number, y = dbinom(date_number, 3, 0.1))

p1 <- ggplot(df, aes(x = x, y = y)) + 
  geom_bar(stat = "identity", col = "hotpink", fill = "hotpink") + 
  scale_y_continuous(expand = c(0.01, 0)) + xlab("x") + ylab("Density") + 
  labs(title = "dbinom(x, 20, 0.5)") + theme_bw(16, "serif") + 
  theme(plot.title = element_text(size = rel(1.2), vjust = 1.5))+
  labs(title="Binomial distribution", x ="Number of dates that went well", y = "Probability")

p1
```

Hay una segunda persona en el show para Silvia. Por el momento pensemos que con la primera persona dos de las tres citas salieron bien.

$$P(Y=2|\theta)= {3\choose 2} \theta^2 (1-\theta)^1$$
Si $\theta=0.1$ entonces $P(Y=2|\theta=0.1)= {3\choose 2} 0.1^2 (0.9)^1=0.027$. Quiere decir que si Silvia no es carismática entonces la probabilidad de tener dos citas buenas es 2.7%. Este es un escenario **inverosímil**. Esta es la palabra correcta en este contexto, vamos a empezar a introducir el concepto de la **función de verosimilitud**. 


Lo que acabamos de hacer es observar un resultado de la distribución binomial  $y_1=2$ y cuando Silvia salga con la segunda persona vamos a observar que $y_2=1$ y así hasta que Silvia haya salido con todos los candidatos. 

### La función de verosimilitud (likelihood function)

**La función de verosimilitud** se define como la probabilidad de los datos dado el modelo. El modelo en este caso es $\theta$, el carisma de Silvia. La probabilidad de los datos es la multiplicación de los resultados de las citas con cada persona, es una multiplicación porque lo que pase con una persona es independiente de lo que pase con la otra. 

Si no tenemos información a priori sobre  $\theta$, ¿qué diriamos sobre el carisma de Silvia si ya observamos 2 de 3 con la primera persona, y 1 de 3 con la segunda?

```{r,eval=FALSE}
binomial_likelihood <- function(theta, data,n) {
  
# `theta` = Probabilidad de que Silvia sea carismática
# `data`  = Datos como resultaron las citas de Silvia
  long<-length(theta)
  z<-rep(0,long)
  for(i in 1:long){
  prob.vect<-sapply(data,FUN=dbinom, size=n, prob=theta[i])
  z[i]<- prod(prob.vect)
  }
  return(z)
}

observed_dates<- c(2,1)
binomial_likelihood(theta=c(0.1,0.3), data=observed_dates,n=3)
```

La verosimilitud es muchas veces criticada o confundida como frecuentista. Esto es incorrecto, la verosimilitud es una función extremadamente útil, que conecta los datos con el modelo y se interpreta como la evidencia que los datos tienen sobre el modelo. 

## Inferencia bayesiana: La distribución posterior

Finalmente vamos a calcular $P(A|B)$ o en el nuevo lenguaje de nuestras variables aleatorias $P(\theta|Datos)$ es la distribución posterior del carisma de Silvia

La distribución posterior es **proporcional** a la verosimilitud multiplicada por la distribución a priori

```{r,eval=FALSE}
# En este gráfico vamos a comparar la a priori, la verosimilitud, y la posterior. 
trial_data <- c(1,2)
number_dates<-3
a<-2
b<-5

d <-
  tibble(theta0 = seq(from = 0, to = 1, length.out = 100)) %>% 
  mutate(`Prior (beta)`           = dbeta(theta0, 
                                          shape1 = a, 
                                          shape2 = b),
         `Likelihood (Binomial)` = binomial_likelihood(theta = theta0,                                                      data=trial_data,n=number_dates),
         `Posterior (beta)`       = dbeta(theta0, 
                                          shape1 = 4, 
                                          shape2 = 7))

glimpse(d)

d %>% 
  gather(key, value, -theta0) %>% 
  mutate(key = factor(key, levels = c("Prior (beta)", "Likelihood (Binomial)", "Posterior (beta)"))) %>% 
  
  ggplot(aes(x = theta0)) +
  # densities
  geom_ribbon(aes(ymin = 0, ymax = value),
              fill = "hotpink") +
  labs(x = expression(theta),
       y = NULL) +
  facet_wrap(~key, scales = "free_y", ncol = 1) +
  theme(panel.grid = element_blank())

```

Perfecto, pero ¿cómo calculo la posterior?


## Inferencia Bayesiana utilizando el algoritmo MCMC

El algoritmo MCMC (Markov Chain Monte Carlo) nos permite optimizar la distribución posterior. En este tutorial vamos a hacer un tipo de MCMC llamado Metropolis-Hastings que es un algoritmo de rechazo. El algoritmo va a proponer un valor para el carisma de Silvia $\theta$ y se pregunta: el valor que propuse ¿incrementa el valor de la probabilidad posterior?- si la respuesta es positiva, nos quedamos este valor, y si es negativa lo rechazamos. 

Manos a la obra

### La propuesta

Vamos a proponer un nuevo parámetro $\theta_{new}$ utilizando una distribución uniforme entre 0 y uno

```{r,eval=FALSE}
proposalfunction <- function(nvals=1){
  unif_val<-runif(nvals,min=0, max=1)
  return(unif_val)}
# Esta es una propuesta de una distribución uniforme. Selecciona valores aleatorios entre 0 y 1
```


### El algoritmo de Metropolis-Hastings

En resúmen el algoritmo de Metropolis-Hastings sigue los siguientes pasos:

1. Empieza con un valor para $\theta$ (``startvalue``) llamado $\theta_0$ 
2. Haz $\theta_{old}=\theta_0$
3. Calcula la distribución posterior $P(\theta_{old}|Datos)$
4. Proponer una distribución aleatoria $g(\theta)$ para obtener un nuevo valor $\theta_new$
5. Calcula la distribución posterior $P(\theta_{new}|Datos)$
6. Calcula los momios $momios=\frac{P(\theta_{new}|Datos)}{P(\theta_{old}|Datos)}$
7. Calcula un valor aleatorio $u$ entre 0 y 1
8. Si $u< momios$ entonces acepta $\theta_{new}$, guárdalo, sino rechaza y no lo guardes.
9. Si lo aceptas haz $\theta_{old}=\theta_{new}$ y vuelve al paso dos hasta acabar las iteraciones. Sino continua al paso dos con el mismo $\theta_{old}$ hasta acabar las iteraciones. 


Vamos a seleccionar un valor para $\theta$ de entrada que lo llamamos  ``startvalue``  y el número de veces que vamos a buscar se llama ``iterations``. El número de iteraciones va a depender de un criterio de convergencia en los pasos 5-7.

```{r,eval=FALSE}
a<-2
b<-5

run_metropolis_MCMC <- function(startvalue, iterations){
    chain = rep (0,iterations+1)
    # Empezamos un un startvalue
    chain[1] = startvalue 
    for (i in 1:iterations){
        theta_old<-chain[i]
# Llamamos a la propuesta, por favor danos un valor de theta nuevo
        theta_new = proposalfunction(1)

#Calculamos los momios (odds)
odds = (binomial_likelihood(theta_new,observed_dates,n=3)*prior_distribution(theta_new, a,b))/
(binomial_likelihood(theta_old,observed_dates,n=3)*prior_distribution(theta_old, a,b)) 
        
        if (runif(1) < odds){  # runif(1) genera una uniforme entre 0 y 1, si este valor es más pequeño que los momios entonces vamos a aceptar el nuevo theta
            chain[i+1] = theta_new
        }else{ 
  # Si no rechazamos y nos quedamos con el theta viejo
            chain[i+1] = theta_old
        }
    }
    return(chain)
}


# Hagamos una corrida y comparemos con tus compañeros
startvalue = 0.3 # Valor inicial 
chain = run_metropolis_MCMC(startvalue, 1) # ¿Qué pasó en una iteración?
chain

startvalue = 0.3 # Ahora corramos 1000 iteraciones
iter=1000
chain = run_metropolis_MCMC(startvalue, iter) # 
head(chain)
```

### Criterio de aceptación

Para aceptar $\theta_{new}$ necesitamos calcular una estadística llamad momios (odds en inglés). Los momios es una razón entre dos probabilidades


### Visualizando el resultado del MCMC para encontrar la convergencia

Para saber si el MCMC llegó al máximo de la distribución posterior, después de varias iteraciones el resultado de las muestras que se aceptan se debe ver como un milpies en un plano (una función que sube y baja todo el tiempo pero se queda horizontal). 

```{r,eval=FALSE}
mcmc <- data.frame(iterations=seq(1,iter+1,1),chain)

# Plot
ggplot(mcmc, aes(x=iterations, y=chain)) +
  geom_line(color="hotpink")+
  labs(title="MCMC run",x="Iterations", y="Posterior distribution")

hist(chain,col="hotpink")

acceptance = 1-mean(duplicated(chain)) # La proporción de los theta que si fueron aceptados, lo hicimos bien o no?
acceptance
```



## Ejercicio extra: Cambiando la propuesta


```{r,eval=FALSE}
## Una nueva propuesta bajo una distribución beta que se ve como U
proposalfunction2 <- function(nvals=1){
  beta_val<-rbeta(nvals,shape1=0.1 ,shape2=0.1) # 
  return(beta_val)}
  
a<-2
b<-5

run_metropolis_MCMC2 <- function(startvalue, iterations){
    chain = rep (0,iterations+1)
    chain[1] = startvalue # Valor inicial del parametro
    for (i in 1:iterations){
        theta_old<-chain[i]
        theta_new = proposalfunction2(1) # Nueva propuesta de una beta
      
        odds = (binomial_likelihood(theta_new,observed_dates,n=3)*prior_distribution(theta_new, a,b))/(binomial_likelihood(theta_old,observed_dates,n=3)*prior_distribution(theta_old, a,b)) 
        
        if (runif(1) < odds){  # La parte de la aceptación 
            chain[i+1] = theta_new
        }else{ # La parte del rechazo
            chain[i+1] = theta_old
        }
    }
    return(chain)
}

startvalue = 0.3 # Valor inicial 
iter=1000
chain = run_metropolis_MCMC2(startvalue, iter) # Do 10,000 steps

mcmc <- data.frame(iterations=seq(1,iter+1,1),chain)

# Plot
ggplot(mcmc, aes(x=iterations, y=chain)) +
  geom_line(color="hotpink")+
  labs(title="MCMC run",x="Iterations", y="Posterior distribution")

hist(chain,col="hotpink")
acceptance = 1-mean(duplicated(chain)) # Proporción de aceptación. Cómo lo hicimos_
acceptance
```

## Referencias

Muchas de estas ideas y código proviene de los siguientes recursos

1. Doing Bayesian Data Analysis in brms and the tidyverse by A Solomon Kurz https://bookdown.org/ajkurz/DBDA_recoded/

2. Dating for Bayesians: Here's How To Use Statistics To Improve Your Love Life by
Andy Kiersz https://www.businessinsider.com/dating-for-bayesians-heres-how-to-use-statistics-to-improve-your-love-life-2013-11
