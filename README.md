# Algoritmo SAL (Simullated Allocation)

Planificación y Operación de Redes.

Máster Universitario en Ingeniería de Telecomunicación.

Universidad de Sevilla.

## Implementación

El algoritmo SAL se ha implementado en matlab, en la función sal. 
Todo el código está comentado. 

La función aplica el algoritmo SAL, y recibe como parámetro, entre otros, los 
nodos de acceso y de núcleo que existen y los costes fijo y variable de estos 
enlaces. Es importante este detalle, ya que la función sal diferencia entre 
nodo de acceso y nodo de núcleo, y entre enlace de acceso y enlace de núcleo. 
Un enlace de acceso sería aquel que une un nodo de acceso con un nodo de núcleo, 
y un enlace de núcleo sería aquel que une dos nodos de núcleo, y ambos tipos de 
enlace se pasan a la función en parámetros separados, así como los nodos de 
núcleo y de acceso.

La función aplica dos restricciones:
 1. El número de enlaces de un nodo de núcleo v no puede ser superior a Gv.
 2. El uso de un enlace e (tanto de acceso como de núcleo) no puede ser superior
    a Me, es decir, un enlace e no podrá cursar tráfico superior a Me. (MeA para
    enlaces de acceso, MeV para enlaces de núcleo). 

Adicionalmente, se llama con un parámetro (M) que indica la profundidad del
algoritmo: Los trayectos que se van a probar y que aparecerán en la solución van 
a contener hasta un máximo de M nodos de núcleo.

Al aplicar el algoritmo, para guardar los trayectos que se van eligiendo se guardan 
sólo los nodos de núcleo que pertenecen a estos, y por otro lado se guarda el número
de trayecto (índice en la matriz de trayectos) por el que va cada unidad de demanda.

La función genera dos ficheros:
 1. Fichero 'debug.txt' que contiene mensajes de debug acerca de cómo el algoritmo
    ha ido iterando.
 2. Fichero 'solucion.txt' que contiene la solución encontrada. Esta solución también
    se imprime por pantalla en la consola de matlab.

## Autor

**Francisco Javier  Ortiz Bonilla** - [fjavierob](https://github.com/fjavierob)
