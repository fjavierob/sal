%
% Francisco Javier Ortiz Bonilla
% Planificaci�n y Operaci�n de Redes
%

%% Funci�n sal: Aplica el algoritmo SAL.
% Par�metros:
    
% A     N�mero de nodos de acceso.
% V     N�mero de candidatos a nodos de n�cleo.
    
% h     Demandas: Matriz Nx3:
%                   - Columna 1: Origen.
%                   - Columna 2: Destino.
%                   - Columna 3: Valor demanda.
    
% ceAf  Costes fijos de los enlaces entre nodos acceso y nodos n�cleo. 
%       Matriz AxV:
%                   - Dimensi�n 1: Nodos de acceso.
%                   - Dimensi�n 2: Nodos de n�cleo.

% ceAv  Costes variables de los enlaces entre nodos acceso y nodos n�cleo. 
%       Matriz AxV:
%                   - Dimensi�n 1: Nodos de acceso.
%                   - Dimensi�n 2: Nodos de n�cleo.
    
% ceVf  Costes fijos de los enlaces entre nodos de n�cleo. 
%       Matriz VxV:
%                   - Dimensi�n 1: Nodos de n�cleo.
%                   - Dimensi�n 2: Nodos de n�cleo.
%       Debe ser una matriz sim�trica, de forma que 
%       ceVf(x,y) = ceVf(y,x);

% ceVv  Costes variables de los enlaces entre nodos de n�cleo. 
%       Matriz VxV:
%                   - Dimensi�n 1: Nodos de n�cleo.
%                   - Dimensi�n 2: Nodos de n�cleo.
%       Debe ser una matriz sim�trica, de forma que 
%       ceVv(x,y) = ceVv(y,x);
    
% M     L�mite enlaces maximo usados en un trayecto. N�mero entero.
    
% Gv    L�mite de enlaces a nodos de acceso desde un nodo n�cleo.
%       Vector 1xV. 

% MeA   L�mite de tr�fico por los enlaces de acceso. 
%       Matriz AxV.

% MeV   L�mite de tr�fico por los enlaces de n�cleo.
%       Matriz VxV sim�trica.
    
% cV    Coste de apertura de los nodos de n�cleo.
%       Vector 1xV.
    
% N     N�mero m�ximo de iteraciones. N�mero entero.
    
% L     L�mite inferior para el coste. N�mero entero.
    
% q     Parametro para decidir si asociar o desasociar. 
%       0.5 < q < 1
function [Fbest, yeA_sol, yeV_sol, trayec_sol, flujos_sol] = sal(A, V, h, ceAf, ceAv, ceVf, ceVv, cV, M, Gv, MeA, MeV, q, N, L)   
    
    % Mensajes de debug en el fichero debug.txt
    debug = fopen('debug.txt', 'w'); 
    % Soluci�n en el fichero solucion.txt
    solucion = fopen('solucion.txt', 'w');
    
    % Estado / soluci�n
    x = 0;
    hSinAsociar = h(:,3)';
    hAsociados  = zeros(1, length(hSinAsociar));  
    H = sum(hSinAsociar);
    yeA       = zeros(A,V);                    % Tr�fico cursado por los enlaces de nodos de acceso a los de n�cleo.
    yeV       = zeros(V,V);                    % Tr�fico cursado por los enlaces de nodos de n�cleo a los de n�cleo.
    trayectos = zeros(1,M);                    % Se guardan los diferentes trayectos por los que van las demandas.
    flujos    = zeros(size(h,1),max(h(:,3)));  % Indica el trayecto por el que va cada unidad de demanda:
                                               %   -> �ndice en la matriz trayectos. 
    
    % Inicio algoritmo.
    n = 0;
    Fbest      = inf;
    yeA_sol    = [];
    yeV_sol    = [];
    trayec_sol = [];
    flujos_sol = [];
    
    fprintf(debug, '* * * * * * * * * * * * * * * * * *\n');
    fprintf(debug, '            Iteraci�n 1            \n');
    fprintf(debug, '* * * * * * * * * * * * * * * * * *\n');

    while (n < N && Fbest > L)
        
        p = rand();
        
        if (p > q)  % Desasociar
            if (x > 0)
                % Escoger demanda aleatoria
                r = randi([1 x]);  % N�mero aleatorio entre 1 y total demanda asignada.
                for a = 1:length(h)
                    suma = sum(hAsociados(1:a));
                    if (r <= suma)
                        break;
                    end
                end
                % Desasociar la demanda a
                fprintf(debug, 'Desasociar demanda %d\n', a);
                fprintf(debug, '------------------------------------------\n');
                origen  = h(a,1);
                destino = h(a,2);
                [yeA, yeV, flujos] = disconnect(origen, destino, a, yeA, yeV, flujos, trayectos, debug);
                hSinAsociar(a) = hSinAsociar(a) + 1;
                hAsociados(a)  = hAsociados(a)  - 1;
                x = x - 1;
            end
        else        % Asociar
            r = randi([1 H-x]);  % N�mero aleatorio entre 1 y total demanda no asignada.
            for a = 1:length(h)
                suma = sum(hSinAsociar(1:a));
                if (r <= suma)
                    break;
                end
            end
                % Asociar la demanda a
                fprintf(debug, 'Asociar demanda %d\n', a);
                fprintf(debug, '------------------------------------------\n');
                origen  = h(a,1);
                destino = h(a,2);
                [trayectos, yeA, yeV, asignado] = allocate(origen, destino, V, M, Gv, MeA, MeV, cV, ceAf, ceAv, ceVf, ceVv, yeA, yeV, trayectos, debug);
                %yeA
                if (asignado)
                    hSinAsociar(a) = hSinAsociar(a) - 1;
                    hAsociados(a)  = hAsociados(a)  + 1;
                    x = x + 1;
                    index = find(flujos(a,:)==0,1,'first');
                    % La unidad de demanda n�mero 'index' de la demanda 'a'
                    % va por el trayecto n�mero 'asignado'.
                    flujos(a,index) = asignado;
                end
        end

        if (x == H) % Todo asignado
            n = n + 1;
            f = fo(yeA, yeV, cV, ceAf, ceAv, ceVf, ceVv);
            if (f < Fbest)
                Fbest      = f;
                yeA_sol    = yeA;
                yeV_sol    = yeV;
                trayec_sol = trayectos;
                flujos_sol = flujos;
                % TODO Guardar soluci�n
            end
            fprintf(debug, '\n => Iteraci�n %d finalizada. Encontrada soluci�n con coste %d\n\n\n', n, f);
            % Resetear estado / soluci�n.
            x = 0;
            hSinAsociar = h(:,3)';
            hAsociados  = zeros(1, length(hSinAsociar));  
            H = sum(hSinAsociar);
            yeA       = zeros(A,V);               
            yeV       = zeros(V,V);               
            trayectos = zeros(1,M);               
            flujos    = zeros(size(h,1),max(h(:,3)));  
            fprintf('. ');
            if (n ~= N)
                fprintf(debug, '* * * * * * * * * * * * * * * * * *\n');
                fprintf(debug, '            Iteraci�n %d           \n', n+1);
                fprintf(debug, '* * * * * * * * * * * * * * * * * *\n');
            end
        end
    end
    
    fprintf('\n');
    
    imprimirSolucion(1, Fbest, trayec_sol, flujos_sol, yeA_sol, yeV_sol, h);
    imprimirSolucion(solucion, Fbest, trayec_sol, flujos_sol, yeA_sol, yeV_sol, h);

end

%% Funci�n para imprimir la soluci�n bien bonita.
function imprimirSolucion(out, Fbest, trayectos, flujos, yeA, yeV, h)

	fprintf(out, '################################################################\n');
	fprintf(out, '# Encontrada soluci�n con coste %d \n', Fbest);
	fprintf(out, '################################################################\n\n');
    pause(1.25);
    
	% Imprimir el reparto de flujos
	fprintf(out, '# REPARTO DE FLUJOS #\n');
	fprintf(out, '################################################################\n\n');
    for d = 1:size(flujos,1)
        origen  = h(d,1);
        destino = h(d,2);
        fprintf(out, 'Demanda %d -> Origen: Nodo acceso %d, Destino: Nodo acceso %d\n', d, origen, destino);
        flujosd = flujos(d,:);
        % Contamos cu�nto tr�fico (y) va por cada trayecto (p) para una demanda (d).
        [y,p] = hist(flujosd,[0 unique(flujosd)]);
        for a = find(p)  % Cogemos los trayectos (p) que no son 0.
            y_ = y(a);
            fprintf(out, '  %d unidades de demanda por trayecto A%d-', y_, origen);
            trayecto = trayectos(p(a),:);
            % Eliminar ceros al trayecto
            u = find(trayecto,1,'last'); % �ltima posici�n que no es cero.
            trayecto = trayecto(1:u);
            for i = 1:length(trayecto)
                if (i ~= length(trayecto))
                    fprintf(out, 'N%d-',trayecto(i));
                else
                    fprintf(out, 'N%d-A%d\n',trayecto(end), destino);
                end
            end
        end
        fprintf(out, '\n');
    end
    pause(1.25);
    
    % Imprimir el tr�fico cursado por los enlaces de accceso
    fprintf(out, '\n');
    fprintf(out, '# TR�FICO CURSADO POR LOS ENLACES DE ACCESO #\n');
	fprintf(out, '################################################################\n\n');
    fprintf(out, 'Acceso[i][j] = Enlace del nodo de acceso i al nodo de n�cleo j\n');
    fprintf(out, '----------------------------------------------------------------\n');
    for i = 1:size(yeA,1)
        for j = 1:size(yeA,2)
            fprintf(out, 'Acceso[%d][%d] = %d\n', i, j, yeA(i,j));
        end
    end
    pause(1.25);
    
    % Imprimir el tr�fico cursado por los enlaces de n�cleo
    fprintf(out, '\n\n');
    fprintf(out, '# TR�FICO CURSADO POR LOS ENLACES DE N�CLEO #\n');
	fprintf(out, '################################################################\n\n');
    k = 1:size(yeV,1);
    for i = 1:size(yeV,1)
        for j = k 
            if (j ~= i)
                fprintf(out, 'Nucleo[%d][%d] = %d\n', i, j, yeV(i,j));
            end
        end
        k(find(k==i)) = [];
    end
    pause(1.25);
    
    % Imprimir los nodos de n�cleo activos
    eActivosA = yeA > 0;
    eActivosV = yeV > 0;
    numEnlacesNodoNucleo = sum(eActivosV,1) + sum(eActivosA,1);
    vActivos = find(numEnlacesNodoNucleo);
    fprintf(out, '\n\n');
    fprintf(out, '# NODOS DE N�CLEO ACTIVOS #\n');
	fprintf(out, '################################################################\n\n');
    fprintf(out, 'Nodos de n�cleo activos: ');
    for i = 1:length(vActivos)
        if (i ~= length(vActivos))
            fprintf(out, '%d, ',vActivos(i));
        else
            fprintf(out, '%d.\n\n',vActivos(end));
        end
	end
    
end


%% Funci�n objetivo.
% Par�metros:
% - yeA:        Tr�fico cursado por los enlaces de acceso.
% - yeV:        Tr�fico cursado por los enlaces de n�cleo.
% - cV:         Coste de apertura de los nodos de n�cleo.
% - ceAf:       Coste fijo de los enlaces de los nodos de acceso.
% - ceAv:       Coste variable de los enlaces de los nodos de acceso.
% - ceVf:       Coste fijo de los enlaces de los nodos de n�cleo.
% - ceVv:       Coste variable de los enlaces de los nodos de n�cleo.
function f = fo(yeA, yeV, cV, ceAf, ceAv, ceVf, ceVv)

    eActivosA = yeA > 0;
    eActivosV = yeV > 0;
    numEnlacesNodoNucleo = sum(eActivosV,1) + sum(eActivosA,1);
    vActivos = numEnlacesNodoNucleo > 0;
    
    c_nodos = sum(vActivos.*cV);
    
    c_enlA  = sum(sum(eActivosA.*ceAf)) + sum(sum(yeA.*ceAv));
        
    c_enlV  = sum(sum(eActivosV.*ceVf)) + sum(sum(yeV.*ceVv));
    % Al calcular el coste entre enlaces de n�cleo, se suman dos veces pues
    % cada enlace se hace dos veces: x -> y  =  y <- x
    c_enlV  = 0.5*c_enlV; 
    
    f = c_nodos + c_enlA + c_enlV;
end


%% Funci�n para desasociar una unidad de demanda entre origen y destino.
% Par�metros:
% - origen:     Origen de la demanda.
% - destino:    Destino de la demanda.
% - demanda:    N�mero de demanda (�ndice de la variable flujos).
% - yeA:        Tr�fico cursado por los enlaces de acceso.
% - yeV:        Tr�fico cursado por los enlaces de n�cleo.
% - flujos:     Variable que indica el trayecto por el que va cada unidad
%               de demanda.
% - trayectos:  Variable que contiene todos los trayectos que se usan.
%
% Salidas: Variables yeA, yeV, flujos actualizadas.
function [yeA, yeV, flujos] = disconnect(origen, destino, demanda, yeA, yeV, flujos, trayectos, debug)

    indices = find(flujos(demanda,:)); % Unidades de demanda asignadas de la demanda.
    indice  = indices(randi([1 length(indices)])); % Escogemos uno aleatorio.
    indice_trayecto = flujos(demanda,indice);
    flujos(demanda,indice) = 0;
    ruta = trayectos(indice_trayecto,:);
    u = find(ruta,1,'last'); % �ltima posici�n que no es cero.
    ruta = ruta(1:u);
    
    % Eliminamos una unidad de demanda de la ruta elegida:
    
    % Actualizamos el tr�fico cursado por los enlaces entre nodos de acceso
    % y n�cleo.
    % Enlace con el nodo de acceso origen.
    
    fprintf(debug, 'Desconectar -> Origen: %d, Destino: %d, Trayecto: ', origen, destino);
    fprintf(debug, '%d ', ruta);
    fprintf(debug, '\n');
        
    yeA(origen,ruta(1))    = yeA(origen,ruta(1))    - 1;
    % Enlace con el nodo de acceso destino.
    yeA(destino,ruta(end)) = yeA(destino,ruta(end)) - 1;
    
    % Actualizamos el tr�fico cursado por los enlaces entre nodos de n�cleo.
    for i = 1:length(ruta)-1
            v1 = ruta(i);
            v2 = ruta(i+1);
            % yeV es una matriz sim�trica.
            yeV(v1,v2) = yeV(v1,v2) - 1;
            yeV(v2,v1) = yeV(v1,v2);
    end
end


%% Funci�n para asociar una unidad de demanda entre origen y destino.
% Par�metros:
% - origen:     Origen de la demanda.
% - destino:    Destino de la demanda.
% - V:          N�mero de candidatos a nodos de n�cleo.
% - M:          L�mite de profundidad.
% - Gv:         L�mite de enlaces a nodos de acceso desde un nodo n�cleo.
% - cV:         Coste de apertura de los nodos de n�cleo.
% - ceAf:       Coste fijo de los enlaces de los nodos de acceso.
% - ceAv:       Coste variable de los enlaces de los nodos de acceso.
% - ceVf:       Coste fijo de los enlaces de los nodos de n�cleo.
% - ceVv:       Coste variable de los enlaces de los nodos de n�cleo.
% - yeA:        Tr�fico cursado por los enlaces de acceso.
% - yeV:        Tr�fico cursado por los enlaces de n�cleo.
% - flujos:     Variable que indica el trayecto por el que va cada unidad
%               de demanda.
% - trayectos:  Variable que contiene todos los trayectos que se usan.
%
% Salidas: 
% - Variables trayectos, yeA, yeV actualizadas.
% - asignado:   Indica el �ndice en la matriz trayectos correspondiente al 
%               trayecto que cursa la unidad de demanda asignada.
%               Vale 0 si no se ha podido asignar la unidad de demanda
%               porque no hab�a ning�n trayecto disponible.
function [trayectos, yeA, yeV, asignado] = allocate(origen, destino, V, M, Gv, MeA, MeV, cV, ceAf, ceAv, ceVf, ceVv, yeA, yeV, trayectos, debug)

    % fprintf('Allocate\n');
    asignado = 0;
   
    % Buscamos la ruta m�s barata entre origen y destino.
    eActivosA = yeA > 0;
    eActivosV = yeV > 0;
    numEnlacesNodoNucleo = sum(eActivosV,1) + sum(eActivosA,1);
    vActivos = numEnlacesNodoNucleo > 0;
    
    trayecto = calcularTrayecto(origen, destino, 1:V, [], 0, M, inf, Gv, MeA, MeV, cV, ceAf, ceAv, ceVf, ceVv, yeA, yeV, vActivos, eActivosV, eActivosA);   
    % trayecto contiene la ruta (2:end) y el coste de esta (1).
    
    coste = trayecto(1);
    if (coste == inf)
        % Si el coste que se obtiene es inf es porque no se ha conseguido
        % encontrar una ruta entre origen y destino.
        fprintf(debug, 'No encontrado ning�n trayecto v�lido \n');
        return
    end
    ruta  = trayecto(2:end);
    fprintf(debug, 'Encontrado mejor trayecto (coste %d): ', trayecto(1));
    fprintf(debug, '%d ', ruta);
    fprintf(debug, '\n');
    
    % Actualizamos el tr�fico cursado por los enlaces entre nodos de acceso
    % y n�cleo.
    % Enlace con el nodo de acceso origen.
    
    %fprintf('Origen: %d, destino: %d, nodo1: %d, nodo2: %d', origen, destino, ruta(1), ruta(end));
    
    
    yeA(origen,ruta(1))    = yeA(origen,ruta(1))    + 1;
    % Enlace con el nodo de acceso destino.
    yeA(destino,ruta(end)) = yeA(destino,ruta(end)) + 1;
    
    % Actualizamos el tr�fico cursado por los enlaces entre nodos de n�cleo.
    for i = 1:length(ruta)-1
            v1 = ruta(i);
            v2 = ruta(i+1);
            % yeV es una matriz sim�trica.
            yeV(v1,v2) = yeV(v1,v2) + 1;
            yeV(v2,v1) = yeV(v1,v2);
    end
  
    % A�adir ceros a ruta para guardarla siempre con la misma longitud.
    r = zeros(1,M);
    r(1:length(ruta)) = ruta;
    ruta = r; % Rellenada con ceros.
    % Guardar el nuevo trayecto
    [tf, index] = ismember(ruta,trayectos,'rows');
    if (~tf)
%         fprintf('\nGuardar ruta..\n');
%         ruta
        trayectos = [trayectos; ruta];
        index = size(trayectos, 1);
    end
    
    % �ndice en la matriz trayectos del trayecto por el que se cursa la
    % unidad de demanda que se ha asignado.
    asignado = index;
end


%% Funci�n para calcular los posibles trayectos para una demanda.
% Es recursiva, llam�ndose a s� misma M*(V+V-1+V-2+..+V-M+1) veces para 
% calcular las diferentes rutas posibles con m = 1, 2,.., M nodos de n�cleo.
% En cada llamada se calcula una ruta de profundidad m y su coste, y lo
% pasa como par�metro a la siguiente llamada. La siguiente llamada vuelve a
% calcular una ruta y coste, compara con la que ha recibido como par�metro
% y pasa la m�s barata a la siguiente llamada.
% Se acaba devolviendo la ruta m�s barata.
% 
% Par�metros:
% - origen:     Origen de la demanda.
% - destino:    Destino de la demanda.
% - posiblesV:  Posibles nodos de n�cleo a los que conectarse.
% - ruta:       Trozo de ruta formada.
% - m:          Profundidad en cada momento.
% - M:          L�mite de profundidad.
% - coste:      Vector con la ruta calculada previamente (posici�n 2:end)
%               y su coste (posici�n 1).
% - Gv:         L�mite de enlaces a nodos de acceso desde un nodo n�cleo.
% - MeA:        L�mite de tr�fico por los enlaces de acceso.
% - MeV:        L�mite de tr�fico por los enlaces de n�cleo.
% - yeA:        Tr�fico cursado por los enlaces de acceso.
% - yeV:        Tr�fico cursado por los enlaces de n�cleo.
% - cV:         Coste de apertura de los nodos de n�cleo.
% - ceAf:       Coste fijo de los enlaces de los nodos de acceso.
% - ceAv:       Coste variable de los enlaces de los nodos de acceso.
% - ceVf:       Coste fijo de los enlaces de los nodos de n�cleo.
% - ceVv:       Coste variable de los enlaces de los nodos de n�cleo.
% - vActivos:   Matriz que indica los nodos   de n�cleo activos (1) o no (0).
% - eactivosV:  Matriz que indica los enlaces de n�cleo activos (1) o no (0).
% - eActivosA:  Matriz que indica los enlaces de acceso activos (1) o no (0).
% 
% Primera llamada desde fuera para obtener el trayecto de menor coste:
%  - posiblesV = V
%  - ruta = []
%  - m = 0
function c = calcularTrayecto(origen, destino, posiblesV, ruta, m, M, coste, Gv, MeA, MeV, cV, ceAf, ceAv, ceVf, ceVv, yeA, yeV, vActivos, eActivosV, eActivosA)
                                      
    % Un nodo de acceso puede conectarse a varios nodos de n�cleo, pero un nodo
    % de n�cleo v s�lo puede conectarse a Gv nodos de acceso.
    
    % m indica la profundidad actual, es decir, el n�mero de nodos de 
    % n�cleo en los trayectos que se van a calcular.
    m = m + 1;
        
    for i = posiblesV
               
        r = [ruta, i];
               
        % Coste total de la ruta de esta iteraci�n.
        c = inf;
        if (rutaValida(origen, destino, r, Gv, MeA, MeV, yeA, yeV, eActivosA, eActivosV))
            c = calcularCoste(origen, destino, r, cV, ceAf, ceAv, ceVf, ceVv, vActivos, eActivosV, eActivosA);
        end
        % Guardamos en un vector el coste y la ruta que se van eligiendo.
        c = [c, r];
                
        pV = posiblesV;
        % Quitamos el nodo i como posible nodo al que conectarse en las
        % siguientes llamadas a la funci�n.
        pV(find(pV==i)) = []; 
        
        % Comparo las rutas mediante el coste y me quedo con la mejor.
        if (coste(1) < c(1))
        	c = coste;
        end
        
        % Llamamos a la funci�n de forma recursiva, hasta profundidad M.
        if (m ~= M)
            coste = calcularTrayecto(origen, destino, pV, r, m, M, c, Gv, MeA, MeV, cV, ceAf, ceAv, ceVf, ceVv, yeA, yeV, vActivos, eActivosV, eActivosA);
            % Nos vamos quedando en cada iteraci�n con la ruta de menor
            % coste.
            if (coste(1) < c(1))
                c = coste;
            end
        end        
    end
end

%% Funci�n para comprobar si un posible trayecto para una demanda cumple las
% restricciones de Gv (n�mero m�ximo de enlaces para el nodo de n�cleo v) y
% Me (tr�fico m�ximo que puede cursar el enlace e). 
% Par�metros:
% - origen:     Origen de la demanda.
% - destino:    Destino de la demanda.
% - posiblesV:  Posibles nodos de n�cleo a los que conectarse.
% - r:          Trayecto a comprobar.
% - Gv:         L�mite de enlaces a nodos de acceso desde un nodo n�cleo.
% - MeA:        L�mite de tr�fico por los enlaces de acceso.
% - MeV:        L�mite de tr�fico por los enlaces de n�cleo.
% - yeA:        Tr�fico cursado por los enlaces de acceso.
% - yeV:        Tr�fico cursado por los enlaces de n�cleo.
% - eactivosV:  Matriz que indica los enlaces de n�cleo activos (1) o no (0).
% - eActivosA:  Matriz que indica los enlaces de acceso activos (1) o no (0).
function valida = rutaValida(origen, destino, r, Gv, MeA, MeV, yeA, yeV, eActivosA, eActivosV)

    valida = 1;

    % Comprobar Gv: m�ximo n�mero de enlaces del nodo de n�cleo v. 
	numEnlacesNodoNucleo = sum(eActivosV,1) + sum(eActivosA,1);
	if (length(r) == 1) % Caso particular: Ruta de un s�lo nodo n�cleo: entrada = salida.
        G = Gv(r);
        if (G ~= inf)
            enlacesNecesarios  = 2 - eActivosA(origen,r) - eActivosA(destino,r);
            enlacesDisponibles = G - sum(eActivosA(:,r));
            if (enlacesDisponibles < enlacesNecesarios)
                valida = 0;
            end
        end
    else % entrada != salida
        for n = 1:length(r)
            v = r(n);
            if (n ~= length(r))
                vnext = r(n+1);
            end
            if (n ~= 1)
                vprev = r(n-1);
            end
            G = Gv(v);
            if (G ~= inf)
                enlacesDisponibles = G - numEnlacesNodoNucleo(v);
                if (v == r(1))
                    % Nodo de entrada
                    enlacesNecesarios = 2 - eActivosA(origen,v) - eActivosV(v,vnext);
                elseif (v == r(end))
                    % Nodo de salida
                    enlacesNecesarios = 2 - eActivosV(v,vprev) - eActivosA(destino, v);
                else
                    % Nodo intermedio
                    enlacesNecesarios = 2 - eActivosV(v,vprev) - eActivosV(v,vnext);
                end
            end
            if (enlacesDisponibles < enlacesNecesarios)
                valida = 0;
                break;
            end
        end
    end
    
       
    if (~valida)
        return
    end
    
    % Comprobar Me: Capacidad del enlace e. Comprobar la capacidad para
    % cada enlace de la ruta.
    
    % Enlaces de acceso.
    if (yeA(origen,r(1)) == MeA(origen,r(1)))
        valida = 0;
        return
    end
    if (yeA(destino,r(end)) == MeA(destino,r(end)))
        valida = 0;
        return
    end
    
    % Enlaces de n�cleo.
    for i = 1:length(r)-1
            v1 = r(i);
            v2 = r(i+1);
            % yeV es una matriz sim�trica.
            if (yeV(v1,v2) == MeV(v1,v2))
                valida = 0;
                break;
            end
	end
    
end


%% Funci�n para calcular el coste de un trayecto para una demanda.
% Par�metros:
% - origen:     Origen de la demanda.
% - destino:    Destino de la demanda.
% - ruta:       Ruta para la que se calcula la demanda.
% - cV:         Coste de apertura de los nodos de n�cleo.
% - ceAf:       Coste fijo de los enlaces de los nodos de acceso.
% - ceAv:       Coste variable de los enlaces de los nodos de acceso.
% - ceVf:       Coste fijo de los enlaces de los nodos de n�cleo.
% - ceVv:       Coste variable de los enlaces de los nodos de n�cleo.
% - vActivos:   Matriz que indica los nodos   de n�cleo activos (1) o no (0).
% - eactivosV:  Matriz que indica los enlaces de n�cleo activos (1) o no (0).
% - eActivosA:  Matriz que indica los enlaces de acceso activos (1) o no (0).
%
% Salida: Variable c que contiene la ruta m�s barata (2:end) y su coste (1).
function coste = calcularCoste(origen, destino, ruta, cV, ceAf, ceAv, ceVf, ceVv, vActivos, eactivosV, eActivosA)
    
    % Comprobamos aqu� si la ruta no cumple con el par�metro G (n�mero
    % m�ximo de enlaces
    
    % Costes de los enlaces de los nodos de acceso origen y destino a los 
    % nodos de n�cleo.
    p = ruta(1);    % Primer nodo de n�cleo del trayecto.
    u = ruta(end);  % �ltimo nodo de n�cleo del trayecto.
    % Consideramos si el enlace ya est� activo a la hora de sumar el coste
    % fijo (ceAf).
    c_origen  = ceAv(origen, p) + (eActivosA(origen, p)==0)*ceAf(origen, p);
    c_destino = ceAv(destino,u) + (eActivosA(destino,u)==0)*ceAf(destino,u);
    
    % Coste apertura de los enlaces de n�cleo. Se contempla que est�n ya
    % abiertos.
    c_apertura_nucleo = 0;
    for v = ruta
        c_apertura_nucleo = c_apertura_nucleo + (vActivos(v)==0)*cV(v);
    end
    % Costes de los enlaces entre los nodos de n�cleo.
    c_enlaces_nucleo  = 0;
    for i = 1:length(ruta)-1
        v1 = ruta(i);
        v2 = ruta(i+1);
        c_enlaces_nucleo = c_enlaces_nucleo + ceVv(v1, v2) + (eactivosV(v1, v2)==0)*ceVf(v1, v2);
    end
    coste = c_origen + c_destino + c_apertura_nucleo + c_enlaces_nucleo;
end


