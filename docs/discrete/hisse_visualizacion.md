---
title: "Visualización resultados del HiSSE"
layout: home
nav_order: 9
index: true
redirect: false
parent: Temario
math: katex
---

Creado por Rosana Zenil-Ferguson (Enero 2025)

## Organizando tus librerias

Para este tutorial necesitaremos los siguientes paquetes. **Importante**: El paquete RevGadgets esta disponible directamente con ``install.packages("RevGadgets")`` pero vamos a revisar una versión especial que es la que grafica mapas estocásticos. Por esa razón vamos a bajar la versión que se indica aquí, así que. sigue las instrucciones como se indican abajo. 

```
# install.packages("devtools")
# install.packages("coda")
# install.packages("phytools")
# install.packages ("ggplot2")
# install.packages ("dplyr")
## Instalar treeio
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##Instalar ggtree
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# Importante: no instales install.packages("RevGadgets") porque esto instala desde CRAN. Queremos una version especifica. Por favor sigue la instruccion de abajo. Si ya habias instalado revgadgets por favor reinicia y reinstala asi

# devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)


library(coda)
library(phytools)
library(ggplot2)
library(dplyr)
library(RevGadgets)
```

## Convergencia del MCMC

Revisemos la convergencia básica

Baja los siguientes archivos
+ [Primera corrida del MCMC](files/hisse_pollination_run_1.log)
+ [Segunda corrida del MCMC](files/hisse_pollination_run_2.log)


```
# Agrega to directorio de trabajo
setwd("~/Dropbox/Teaching/Workshops/UNAM2025/discrete_trait/files")

# Con el paquete coda
mcmc_run1 <- readTrace(path ="hisse_pollination_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)
```

Ahora hagamos un mejor trabajo con los gráficos de convergencia

```
# Esta parte se hace con ggplot2 y con RevGadgets
mcmc_run1 <- readTrace(path ="hisse_pollination_run_1.log", burnin = 0.1)
mcmc_run1<-data.frame(mcmc_run1)
mcmc_run1<- cbind(mcmc_run1,run=rep("run 1",length(mcmc_run1$Iteration)))

mcmc_run2 <- readTrace(path ="hisse_pollination_run_2.log", burnin = 0.1)
mcmc_run2<-data.frame(mcmc_run2)
mcmc_run2<- cbind(mcmc_run2,run=rep("run 2",length(mcmc_run2$Iteration)))

mcmc_table<-rbind(mcmc_run1,mcmc_run2)

trace_plot<- ggplot(mcmc_table, aes(x=Iteration,y=Posterior,group=run))+
              geom_line(aes(color=run))+
              theme_classic()
trace_plot

#Aqui se nota que hay que cortar un poco mas de burn-in (15,000 generaciones), asegúrense siempre de cortar mas valores si hace falta
```


## Estimadores de las tasas

Queremos graficar y obtener estadísticas resumen de los parametros de interés. Recordemos que para el modelo hisse los parámetros que estimamos van a ser los siguientes:

+ Tasas de transición $$q_{01A}, q_{10A}, q_{01B}, q_{10B}$$
+ Tasas de transición entre estados escondidos $$\alpha, \beta$$ (se llaman hidden_rate en el código de RevBayes)
+ Tasas de diversificación neta $$r_{0A}=\lambda_{0A}-\mu_{0A},r_{1A}=\lambda_{1A}-\mu_{1A},r_{0B}=\lambda_{0B}-\mu_{0B}, r_{1B}=\lambda_{1B}-\mu_{1B}$$
+ Alternativo: la fracción de extinción $$\epsilon_{0A}=\mu_{0A}/\lambda_{0A}, \epsilon_{1A}=\mu_{1A}/\lambda_{1A}, \epsilon_{0B}=\mu_{0B}/\lambda_{0B},\epsilon_{1B}=\mu_{1B}/\lambda_{1B}$$
+ Frecuencias de la raíz

Para los modelos más complicados siempre prefiero graficar violines, y tener más control en los colores y la figura, por eso mi selección de ggplot2. La comparación de violines es super útil.

## Graficando tasas de transición

```
# Esta parte se hace con ggplot2

traitcols<-c("#3D348B","#7678ED","#F18701", "#F35B04")

hisse<- read.table("hisse_pollination_run_1.log", header=TRUE)
hisse<- hisse[-seq(1,15000,1),] # nunca olvidemos quitar el burn-in con esta instruccion tengo mas contro de cuantas iteraciones remuevo.

transition_rates<- data.frame(dens=c(hisse$q_01A,hisse$q_01B, hisse$q_10A, hisse$q_10B) ,rate=rep(c("q_01A","q_01B","q_10A","q_10B"),each=length(hisse$q_01A)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Tasas de transicion")+
  scale_fill_manual( values = traitcols)+
  theme_classic()
violin_transitions

## Transiciones entre estados escondidos

hidden_rates<- data.frame(dens=c(hisse$hidden_rate1,hisse$hidden_rate2) ,rate=rep(c("alpha","beta"),each=length(hisse$hidden_rate1)))

violin_hidden<- ggplot(hidden_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Tasas de transición entre estados escondidos ")+
  scale_fill_manual( values = traitcols[1:2])+
  theme_classic()
violin_hidden
```

## Graficando tasas de diversificación


```
# Esta parte se hace con ggplot2

divcols<-c("#E63946","#1D3557","#F3A5AB", "#457B9D")

# Recuerda de RevBayes 1=0A, 2=1A, 3=0B, and 4=1B
netdiversification_rates<- data.frame(dens=c(hisse$speciation.1.-hisse$extinction.1.,hisse$speciation.2.-hisse$extinction.2., hisse$speciation.3.-hisse$extinction.3.,hisse$speciation.4.-hisse$extinction.4.) ,rate=rep(c("net_div_0A","net_div_1A","net_div_0B","net_div_1B"),each=length(hisse$speciation.1.)))

netdiversification_rates$rate<- factor(netdiversification_rate$rate,levels=c("net_div_0A","net_div_1A","net_div_0B","net_div_1B"))

violin_diversification<- ggplot(netdiversification_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Tasas de diversificación neta")+
  scale_fill_manual( values = divcols)+
  theme_classic()
violin_diversification
```

**¿Qué podemos concluir de este gráfico?**

## ¿Prueba de hipótesis: La polinización esta asociada a la diversificación?

Una de las ventajas más importante de las distribuciones posteriores es que rápidamente podemos concluir si dos parámetros son iguales o distintos. Formalicemos entonces la prueba de hipotésis para saber si el caracter esta asociado a la diversificación.

Construimos dos estadísticas de resumen que son las diferencias entre la diversificación neta de 0 y 1 para el estado A y el estado B respectivamente:

+ $$T_A= (\lambda_{0A}-\mu_{0A})-(\lambda_{1A}-\mu_{1A})$$ 
+ $$T_B= (\lambda_{0B}-\mu_{0B})-(\lambda_{1B}-\mu_{1B})$$ 

Si estas diferencias son 0 con alta probabilidad esto significa que la polinización y sus estados no hacen la diferencia en el proceso de diversificación.

```
# Esta parte se hace con ggplot2

difcols<-c("#FF006E","#FFC2DC")
T_diferencias<- data.frame(dens=c((hisse$speciation.1.-hisse$extinction.1.)-(hisse$speciation.2.-hisse$extinction.2.),(hisse$speciation.3.-hisse$extinction.3.)-(hisse$speciation.4.-hisse$extinction.4.)),diferencia=rep(c("T_A","T_B"),each=length(hisse$speciation.1.)))

violin_difference<- ggplot(T_diferencias,aes(x=diferencia,y=dens, fill=diferencia))+
  geom_violin(trim=FALSE)+
  labs(title="Estadisticas resumen")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

# Observemos que 0 cruza la diferencia entre las tasas lo que quiere decir que pueden ser iguales
violin_difference
```

Formalizando la prueba de hipótesis
$$ H_0: T_A=0 \textrm{ y } T_B=0$$

```
# Esta parte se hace con dplyr

cuantiles <- T_diferencias %>% group_by(diferencia)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
cuantiles
```

Observamos que el intervalo de credibilidad del 95% de la distribución de $$T_A$$ es (-0.219, 0.425) y el intervalo de credibilidad del 95% de la distribución de $$T_B$$ es (-0.554, 0.0766). Cómo 0 pertenece a estos intervalo entonces $$P(H_0\lvert Datos)>0.05$$ lo que significa que el tempo de la diversificación para los estados 0 y 1 es igual. **Importante**: Noten que en las pruebas de hipótesis bayesiana no utilizo, p-valores, ni significancia, ni rechazo, solo probabilidad e intervalos de credibilidad. Esto sucede porque la manera en que interpretamos probabilidad en bayesiana es distinto. Tengan cuidado cuando interpreten. 

**Segunda parte**: ¿Qué pasa con los estados escondidos?

Construimos dos estadísticas de resumen que son las diferencias entre la diversificación neta de A y B para el estado 0 y el estado 1 respectivamente:

+ $$T_0= (\lambda_{0A}-\mu_{0A})-(\lambda_{0B}-\mu_{0B})$$ 
+ $$T_1= (\lambda_{1A}-\mu_{1A})-(\lambda_{1B}-\mu_{1B})$$ 

Si estas diferencias son 0 con alta probabilidad esto significa que hay cambios en el tempo de la diversificación pero se deben a algo más que nosotros no medimos en la filogenia.


```
# Esta parte se hace con ggplot2

difcols<-c("#1DB32C","#BFEEC3")
T_diferencias<- data.frame(dens=c((hisse$speciation.1.-hisse$extinction.1.)-(hisse$speciation.3.-hisse$extinction.3.),(hisse$speciation.2.-hisse$extinction.2.)-(hisse$speciation.4.-hisse$extinction.4.)),diferencia=rep(c("T_0","T_1"),each=length(hisse$speciation.1.)))

violin_difference<- ggplot(T_diferencias,aes(x=diferencia,y=dens, fill=diferencia))+
  geom_violin(trim=FALSE)+
  labs(title="Estadisticas resumen")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

# Observemos que 0 cruza la diferencia entre las tasas lo que quiere decir que pueden ser iguales
violin_difference

```

Formalizando la prueba de hipótesis
$$ H_0: T_0=0 \textrm{ y } T_1=0$$


```
# Esta parte se hace con dplyr

cuantiles <- T_diferencias %>% group_by(diferencia)%>%reframe(res=quantile(dens,probs=c(0.025,0.975)))
cuantiles
```

Observamos que el intervalo de credibilidad del 95% de la distribución de $$T_0$$ es (-0.488, -0.214) y el intervalo de credibilidad del 95% de la distribución de $$T_1$$ es (-0.989, -0.321). Cómo 0 **no** pertenece a estos intervalo entonces $$P(H_0\lvert Datos)<0.05$$ lo que significa que el tempo de la diversificación para A y B es diferente. Este resultado unido al resultado anterior indican que el tempo cambia debido a algo más que nosotros no medimos pero se manifiesta en el árbol. Estos resultados son equivalentes ajustar el modelo de caracter independiente CID y seleccionarlo, indicando que la diversificación es independiente del caracter de interés. 



## ¿Qué hubieras concluido con BiSSE?

Por que la curiosidad mató al gato chequemos que paso con BiSSE

Baja los siguientes archivos
+ [Primera corrida del MCMC](files/bisse_pollination_run_1.log)
+ [Segunda corrida del MCMC](files/bisse_pollination_run_2.log)


```
# Esta parte se hace con ggplot2
bisse<- read.table("bisse_pollination_run_2.log", header=TRUE)
bisse<- bisse[-seq(1,15000,1),] 
divcols<-c("#E63946","#1D3557")

# Recuerda de RevBayes 1=0A, 2=1A, 3=0B, and 4=1B
netdiversification_rates<- data.frame(dens=c(bisse$net_diversification.1., bisse$net_diversification.2.) ,rate=rep(c("net_div_0","net_div_1"),each=length(bisse$net_diversification.1.)))

violin_diversification<- ggplot(netdiversification_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE)+
  labs(title="Tasas de diversificación neta")+
  scale_fill_manual( values = divcols)+
  theme_classic()
violin_diversification
```

Prueba estadística

```
difcols<-c("#1DB32C","#BFEEC3")
T_diferencia<- data.frame(dens=(bisse$net_diversification.1.-bisse$net_diversification.2.),diferencia=rep("T",each=length(bisse$net_diversification.1.)))

violin_difference<- ggplot(T_diferencia,aes(x=diferencia,y=dens, fill=diferencia))+
  geom_violin(trim=FALSE)+
  labs(title="Estadìstica resumen")+
  scale_fill_manual( values = difcols)+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic()

# Observemos que 0 cruza en una probabilidad muy pequeña de la distribución la diferencia entre las tasas lo que quiere decir que pueden ser iguales
violin_difference

## 0 no esta en el intervalo de probabilidad del 95%
quantile ((bisse$net_diversification.1.-bisse$net_diversification.2.),probs=c(0.025,0.975))
```

**¿Què hubieras concluido para BiSSE? ¿Cuál es la mayor falla de este modelo?**
