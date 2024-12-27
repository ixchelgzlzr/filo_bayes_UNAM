---
title: Antes del taller
layout: home
nav_order: 2
index: true
---


## Prepara tu computadora para el taller

Para que las sesiones practicas del taller transcurran sin contratiempos por favor instala los siguientes programas **con anticipación**. No dudes en contactarnos si tienes problemas con la instalación. 

### a) RevBayes

Utilizaremos [RevBayes](https://revbayes.github.io) para todos los tutoriales. Este software utiliza un lenguaje similar a R, llamado `Rev` que está enfocado en estadística filogenética Bayesiana. A diferencia de otros programas que utilizan GUIs (interfaces gráficas, por sus siglas en inglés), RevBayes requiere que los usuarios especifiquen el **modelo evolutivo** y **las características del análisis** usando código. Hacer esta especificación manual tiene la gran ventaja de proveer versatilidad en el tipo de modelos y análisis que se pueden utilizar. Sin embargo, la curva de aprendizaje en el uso de este software puede resultar intimidante. Uno de los objetivos de este taller es que los participantes se familiaricen con el código utilizado en RevBayes utilizando algunos ejemplos de análisis evolutivos. Por favor asegúrate de tener la **version más reciente** de RevBayes (v1.2.5). 

1) Descarga el archivo ejecutable (executable) correspondente a tu sistema operativo en [esta página](https://revbayes.github.io/download). 

2) Coloca la carpeta descargada en un sitio cerca de la raíz de tu computadora, por ejemplo en 
```
/Users/ixchel/Documents/revbayes-v1.2.5
```

3) Recomendado pero no requerido: Agrega el directorio donde guardaste el ejecutable de revbayes (e.g. ```/Users/ixchel/Documents/revbayes-v1.2.5/bin/```) al PATH de tu terminal. Esto permitirá que puedas abrir revbayes desde cualquier folder en tu terminal. 

En Mac o Linux, substituye ```/my/path``` por el directorio donde está tu ejecutable. Ojo! checa cual es el nombre del archivo donde está tu PATH, por ejemplo puede ser ```~/.bash_profile```, ```~/.zshrc```, o ```~/.zsh```, utiliza el que corresponda a tu computadora.
```
echo 'export PATH=/my/path:$PATH' >> ~/.bash_profile
```

### b) Algún editor de texto 

Recomiendo [Visual Studio Code](https://code.visualstudio.com/Download), pero otras opciones son [Sublime text](https://www.sublimetext.com/3), o [NotePad++](https://notepad-plus-plus.org).


### c) Tracer
Utilizaremos este programa para visualizar las MCMCs, lo puedes descargar [aquí](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

### d) FigTree
Utilizaremos este programa para visualizar los árboles, lo puedes descargar [aquí](https://github.com/rambaut/figtree/releases)
