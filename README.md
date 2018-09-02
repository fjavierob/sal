# Algoritmo SAL (Simullated Allocation)

Planificaci�n y Operaci�n de Redes.

M�ster Universitario en Ingenier�a de Telecomunicaci�n.

Universidad de Sevilla.

## Implementaci�n

El algoritmo SAL se ha implementado en matlab, en la funci�n sal. 
Todo el c�digo est� comentado. 

La funci�n aplica el algoritmo SAL, y recibe como par�metro, entre otros, los 
nodos de acceso y de n�cleo que existen y los costes fijo y variable de estos 
enlaces. Es importante este detalle, ya que la funci�n sal diferencia entre 
nodo de acceso y nodo de n�cleo, y entre enlace de acceso y enlace de n�cleo. 
Un enlace de acceso ser�a aquel que une un nodo de acceso con un nodo de n�cleo, 
y un enlace de n�cleo ser�a aquel que une dos nodos de n�cleo, y ambos tipos de 
enlace se pasan a la funci�n en par�metros separados, as� como los nodos de 
n�cleo y de acceso.

La funci�n aplica dos restricciones:
 1. El n�mero de enlaces de un nodo de n�cleo v no puede ser superior a Gv.
 2. El uso de un enlace e (tanto de acceso como de n�cleo) no puede ser superior
    a Me, es decir, un enlace e no podr� cursar tr�fico superior a Me. (MeA para
    enlaces de acceso, MeV para enlaces de n�cleo). 

Adicionalmente, se llama con un par�metro (M) que indica la profundidad del
algoritmo: Los trayectos que se van a probar y que aparecer�n en la soluci�n van 
a contener hasta un m�ximo de M nodos de n�cleo.

Al aplicar el algoritmo, para guardar los trayectos que se van eligiendo se guardan 
s�lo los nodos de n�cleo que pertenecen a estos, y por otro lado se guarda el n�mero
de trayecto (�ndice en la matriz de trayectos) por el que va cada unidad de demanda.

La funci�n genera dos ficheros:
 1. Fichero 'debug.txt' que contiene mensajes de debug acerca de c�mo el algoritmo
    ha ido iterando.
 2. Fichero 'solucion.txt' que contiene la soluci�n encontrada. Esta soluci�n tambi�n
    se imprime por pantalla en la consola de matlab.

## Autor

**Francisco Javier  Ortiz Bonilla** - [Pogorelich](https://github.com/pogorelich)