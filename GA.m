function [BestFit,BestPobInd]=GA(DataProject,RRHH,TA,CP,CAType,PDim,CI,Par) 
% GA function -> Genetic Algorithm
% Autor: Lic Bolivar E. Medrano Broche. Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n: Algoritmo Genetico
% Entradas
% DataProject es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Informaci??n sobre los recursos 
%  TA -> Tiempos de arribo de los proyectos  TA=[TA1,TA2,...,TAn]
%  CP -> Longitud de las Rutas Criticas
%  Par -> [COP,MOP,RC,RM,RE] parametros de AG y de disenho
%  OPC -> [1 2 3] operadores de cruce
%  OPM -> [1 2]  operadores de mutacion
%  RC  -> Raz?n de cruce
%  RM  -> Razon de mutacion RM=1-RC
%  RE  -> Raz?n de elitismos para el caso del oprador de cruce 3 
% Salidas
% BestFit->  fitness del mejor individuo
% BestInd-> mejor individuo

% INICIALIZACI??N 

% Inicializacion de los parametros

 RC=Par(1); RM=Par(2);

 OPC=Par(3); %OPM=1;%Par(4); 
 
 RE=Par(4); % Raz?n de elitismo (si )

% Inicializaci??n 

[Pob,BestPobInd,NP,CA] = PobCreator(DataProject,RRHH,TA,PDim,CP,CAType); % Crea  la poblacion inicial

for i=1:CI % Para cada una de las  generaciones O(n)           
    
  
         %SELECI?N
         [PoolMating,MutInd]=Selection(Pob,PDim,RC);
           
         %CRUCE                          
         [Offpring]=CrossOver(PoolMating,CA,OPC,RC*PDim,RE*PDim,NP); %
         
         % MUTACI?N 
         if Par(3)==3
            
             OPM=2;
         
         else

             OPM=1;
        
         end
         
         [IndMut]=Mutation(DataProject,RRHH,TA,CP,CAType,MutInd,NP,CA,RM*PDim,OPM); %

         Pob=cat(3,Offpring,IndMut);

     for j=1:PDim % Para cada una de los individuos O(n)   
         
        [Pob(j)] = SSGSGA(Pob(j),DataProject,RRHH,TA,CP,NP,CA);
      
        if Pob(j).FO < BestPobInd.FO 
            
            BestPobInd=Pob(j);
           
        end
       
     end % For -> para cada part??cula 

end % For -> cada generaci?n
BestFit=BestPobInd.FO;
end

function [Pob,BestPobInd,NP,CA] = PobCreator(DataProject,RRHH,TA,Pdim,CP,CAType)
% PobCreator function ->  craa la poblacion inicial 
% Autor: Lic Bolivar E. Medrano Broche. o. Ciencias Basicas.
% Facultad III. Universidad de las Ciencias Informaticas
% Descripci??n: Crea el enjambre a partir de las entradas 
% Entradas
% DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Informaci??n sobre los recursos 
%  TA -> Tiempos de arribo de los proyectos  TA=[TA1,TA2,...,TAn]
%  Pdim-> Cantidad de  de indivisuos 
% SALIDAS
% Pob-> Martiz que contiene al la poblaci?n
% BestPobIndPB -> Mejor individuo.
% NP -> Cantidad de proyectos
% CA -> Cantidad de Actividades

% IMPLEMENTACI??N

Pob=[];

% Generacion aleatoria de los pesos para la asignacion de recursos
 
for i=1:Pdim
    
 p=unidrnd(9);

 w1=p/10; aux=10-p;

if aux==1
    
  w2=0.05; w3=0.05;
  
else
   
    aux1=unidrnd(aux-1);  w2=aux1/10;  w3=(aux-aux1)/10;
end

   DMR=[w1,w2,w3]; 
    
   [Schedule,FV,NP,CA] = SSGS(DataProject,RRHH,TA,CP,CAType,DMR);

   Ind=struct('TI',Schedule(:,:,1),'TF',Schedule(:,:,2),'RA',Schedule(:,:,3),'Pr',Schedule(:,:,4),'FO',FV); % Estructura del individuo 
   
   Pob=cat(3,Pob,Ind); % Crea poblaci?n inicial
       
end
 

[aa,IdBestInd]=min([Pob(1:Pdim).FO]);  %  Valor de la funcion objetivo del mejor individuo ; BestPobInd indivisuo de mejor evalucion 

aux=length(IdBestInd);

if  aux~=1

 p=unidrnd(aux);
 
 BestPobInd=IdBestInd(p);
 
else
    
  BestPobInd=Pob(IdBestInd); % mejor individuo

end

end

function [Schedule,FO,NP,CA]= SSGS(DataProject,RRHH,TA,CP,CAType,DMR)
% SSGS Summary of this function [Schedule,FO,CA,NP] = SSGSpdNew(DataProject,RRHH,TA,CP,W)
% Autor: Lic Bolivar E. Medrano Broche. Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n: Heuristica para la generaci??n de soluciones iniciales factibles basada
% en Serial Scheme Generation Sechdeule (SSGS); con enfoque Multi-Proyecto, 
% la prioridad entre los proyecto cambia en el tiempo, la prioridad entre las actividades 
% est?? dada por la regla de prioridad MTS-modificada  y prioridad para la
% asignaci??n de los recursos.
% Dise??ada para  n proyectos n=2,3,... con la misma cantidad de actividades
% con tiempos de arribo Start=[T1,T2,T3,T4,...,Tn]
% los proyectos son tomados con dos actividades ficticias
% ENTRADA
%  DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Informaci??n sobre los recursos 
%  TA -> Tiempos de arribo de los i proyectos  TA=[TA1,TA2,...,TAn]
%  CP -> Rutas criticas de los i proyectos     CP=[CP1,CP2,...,CPn]
% SALIDAS
%  Schedule -> Soluci??n  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  TI -> Tiempos de Inicio
%  TF -> Tiempos de Finalizaci??n 
%  RA -> Recurso Asignado                                         
% 

% DEFINICI??N DE VARIABLES
%tic

[CA,H,NP]=size(DataProject); % Datos de los Proyectos (CA -> cantidad de actividades, NP-> cantidad de Proyectos)

Schedule=cat(3,zeros(NP,CA),zeros(NP,CA),zeros(NP,CA),rand(NP,CA)); % Estructura de la soluci??n

PS=zeros(1,NP); % Variable que guarda los proyectos que han comenzado PS(i)={0,1} 

PF=PS;      % Variable que guarda los proyectos que han finalizado PF(i)={0,1} 

APP=zeros(NP,CA);  % Variable que guarda las actividades que se han planificado APP(i,j)={0,1} 

% INICIALIZACI??N

PS(1)=1; %  Indica que el proyecto 1 comienza procesandose

APP(1,1)=1; %  Se planifica la primera actividad del primer proyecto

ARP=APP; % Conjuntos de actividades que pueden ser planificadas se inicializa con  la actividad 1 del proyecto 1  
 
Stop=0; %  Se inicializa en cero la condici??n de parada

ind=2:NP; % indice de los proyectos 

ind2=1:NP;

v=1; % Contador de los proyectos que inician

RDE=[];

while Stop < NP % Mientras la cantidad de proyectos finalizado sea menor que la actidad de proyectos 
    
if v~=NP 

    MK=max(max(Schedule(:,2:end,2)));
            
for k=nonzeros(ind)' % Para cada uno de los Proyectos
    
    if TA(k)<=MK && TA(k)~=0 % Si el mayor si el Makespan parcial MK es mayor o igual a la fecha de arribo deque la fe
          
       PS(k)=1;
       
       Schedule(k,1,1)=TA(k); % se inician los proyectos que aun no han comenzado
       
       Schedule(k,1,2)=TA(k);
              
       v=v+1;
       
       ind(k-1)=0;
       
       ARP(k,1)=1;
       
       APP(k,1)=1;
                    
     end
       
end
    
end
       
 PR=find(PS==1);% Proyectos que han comenzado 
    
        CPS=length(PR); % Cantidad de proyectos iniciados
   
     if CPS==1 % Si la cantidad de proyectos iniciados es 1 
         
         OrdenP=PR; % Orden es PR
         
    else % en caso contrario
        
         OrdenP=zeros(1,length(PR));   
         
         MKP=OrdenP;
                         
         for l=PR
           
             MKP(l)=(CP(l)-(max(Schedule(l,:,2))-TA(l)))/CP(l);
                                                                       
         end
         
         aux=MKP(PR)';
                
         ORD=sortrows([aux,PR'],1);
         
        OrdenP=ORD(:,2)';
         
     end  

    for i=OrdenP % Para proyectos Iniciado segun OrdenP
                                                                 
         [Schedule,APP,ARP,RDO] =SFMtsRa(Schedule,DataProject,RRHH,i,nonzeros(ARP(i,:))',ARP,APP,RDE,CAType,DMR);                                                          
               
         RDE=RDO;             
  
    end
        
  for r=nonzeros(ind2)'
    
    if Schedule(r,CA,2)~=0 && sum(APP(r,:))==CA
        
        PF(r)=1;
        
        PS(r)=0;   
        
        ind2(r)=0;
        
        ARP(r,:)=0;
             
    end
    
  end
  
  Stop=sum(PF);
  
end  

FO=mean(((Schedule(:,CA,2)'-(CP+TA))));  %mean((Schedule(:,CA,2)'./CP)-1)); %mean(Schedule(:,CA,2));

%toc
end

function [Schedule,APP,ARP,RDO] = SFMtsRa(Schedule,DataProject,RRHH,i,ASched,ARP,APP,RDO,CAType,DMR)
%  SchedulingFeasible Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
%  Planifica las activitades de ActReadyPlan sin violar las restricciones de
% Precedencia y de Recursos.
% ERNTRADAS
% Schedule -> soluci??n 
% DataProyect -> Datos de los Proyectos 
%ActReadyPlan-> Conjunto de actividades que est??n lista para ser planificadas 
% RRHH-> Datos sobre los recursos humanos del proyecto
% APP(i) -> Conjunto de Actividades planificadas por Proyecto.
% DataProject -> Informaci???n sobre los proyectos.
% SALIDAS
% Schedule-> soluci??n 
% APP(i) -> Actividades planificadas por Proyecto.

% Parte 1 -> INICIALIZACI??N

 ARP(i,:)=0;
 
  TR=0;
      
  k=0;
   
  SI=DataProject(:,4:6,i);
 
  ASASched=SI(ASched,:);                   % Act-> Conjunto de actividades  sucesoras inmediatas de las actividades del conjunto ASched j    
 
  Act=nonzeros(ASASched(1:numel(ASASched)))';
  
 TipoAct=(DataProject(Act,2,i))';     %TipoAct-> Tipos de las actividades del 

 WAct=Schedule(i,Act,4);                      % WAct -> Prioridades del conjunto de activiades Act
 
 OrdenAct=sortrows([Act;TipoAct;WAct]',3);   % OrdenAct-> Actividades Ordenadas segun WAct
 
 Act=OrdenAct(:,1)';

 ModeProc=DataProject(:,7:9,i);               %  ModeProc-> guarda la duraci??n de las actividades seg??n la especilizaci??n del recurso 

 aux=length(Act);                             % aux -> Cantidad de actividades de conjunto Act
 
 
 while k < aux
     
     TFS=Schedule(:,:,2);
     
     k=k+1;
     
     [AA,p]=find(DataProject(:,[4 5 6],i)==Act(k)); %#ok<NASGU> % Conjunto de actividades antesesoras de la actividad Act(k) -> antecesora de j 
          
      CAP=sum(APP(i,AA));                           % Calcula la cantidad de actividades antecesora ya planificadas de
      
      if CAP==length(AA) && APP(i,Act(k))==0        % si la actividad Act(k) est?? lista y nunca ha sido planificada
           
          ARP(i,Act(k))=Act(k);                     % Se inserta Act(k) en las actividades planificadas para se usada de indice
            
           APP(i,Act(k))=1;                         % Se archiva Act(k) como planificada
          
            TipoR=DataProject(Act(k),2,i);          % Definir el tipo  de recurso que requiere la actividad j del proyecto i 
           
            if TipoR==0 % Si la actividad j require un tipo de recurso=0 => ??ltima actividad del proyecto 
       
            Schedule(i,Act(k),1)=max(Schedule(i,AA',2));       % Tiempo de inicio de la actividad j=CA
               
            Schedule(i,Act(k),2)=Schedule(i,Act(k),1);         %  Tiempo de finalizaci??n de la actidad j=CA 
       
            else
                
                if and(k==1,isempty(RDO)==1)==1 || TR~=TipoR || isempty(RDO)==1
              
                TR=TipoR;
           
                [RDO]=RrHhAsig(CAType,Schedule,RRHH,TipoR,DMR);
               
                end
                                
                TFAmr=TFS(Schedule(:,:,3)==RDO(1));
                                                               
                if isempty(TFAmr)==1
                    
                    % Asignar Recursos a la actividad Act(k)
                    
                   Schedule(i,Act(k),3)=RDO(1);
                                                                        
                   Schedule(i,Act(k),1)=max(Schedule(i,AA,2));
                         
                   Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RDO(1),3));
                   
                   RDO=RDO(2:length(RDO));                            
                                   
                else
                     Schedule(i,Act(k),3)=RDO(1);
                    
                     Schedule(i,Act(k),1)=max(max(TFAmr),max(Schedule(i,AA',2))); % Calcula tiempo de inicio de la actividad j del proyecto i
            
                     Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RDO(1),3)); % 
                     
                     RDO=RDO(2:length(RDO));
                     
                end
                                            
            end
           
      end
     
 end

end

function [RDO]=RrHhAsig(CAType,Schedule,RRHH,TipoR,DMR)

CTypeR=RRHH(RRHH(:,2)==TipoR,1)'; % Conjunto de recursos  igual a TipoR de recurso que requiere la actividad que esta lista para procesarce

Coef=zeros(1,length(CTypeR)); MK=max(max(Schedule(:,:,2))); TF=Schedule(:,:,2);

RA=Schedule(:,:,3);
k=0;
for r=CTypeR % Para todos los recursos de TipoR
     k=k+1;
    asig=(RA==r);
    
    CAr=sum(sum(asig));
    
    if  CAr==0
        switch RRHH(r,3)
            case 1
                 P3=(5+RRHH(r,3)/10);
            case 2
                 P3=(5+RRHH(r,3)/10);
            case 3
                 P3=(5+RRHH(r,3)/10);
        end
        Coef(k)=P3;
    else
        LastR=max(max(TF(RA==r)));
       switch RRHH(r,2)
            case 1
                P1=1-CAr/CAType(1); P2=1-LastR/MK;
            case 2
                P1=1-CAr/CAType(2); P2=1-LastR/MK;
            case 3
                P1=1-CAr/CAType(3); P2=1-LastR/MK;
           case 4
                P1=1-CAr/CAType(4); P2=1-LastR/MK;
       end
       
       P3=(1-RRHH(r,3)/10);
       
       %[RRHH(r,1);RRHH(r,3);DMR(1)*P3;DMR(2)*P2;DMR(3)*P1]
       
       Coef(k)=(DMR(3)*P1+DMR(2)*P2+DMR(1)*P3);
    end
          
end

    aux=sortrows([Coef;CTypeR]',-1); % Recursos ordenados segun Coeficiente de asignacion 

    RDO=aux(:,2)';
    
end

function [PoolMating,MutInd]=Selection(Pob,PDim,RC)
% Selection: -> Operador de seleccion (Seleccion por torneo)
% Entradas
% Pob : poblacion
% Po:
% NP: Numero de proyectos
% MO:
% RC: Razon de cruce
%Salidas
% PoolMating (piscina de cruce)
% MutInd (individuos que ser??n mutados)
IdFit=zeros(PDim,2);

IdFit(:,1)=(1:PDim)';

for i=1:PDim
    IdFit(i,2)=Pob(i).FO;
 end

IdFit=sortrows(IdFit,2);

IdPoolMating=IdFit(1:PDim*RC,1);

IdMut=IdFit(PDim*RC+1:end,1);

PoolMating=Pob(IdPoolMating);

MutInd=Pob(IdMut); 

end

function [Offspring]=CrossOver(PoolMating,CA,OPC,CDim,EDim,NP)
                            % (PoolMating,CA,OPC,RC*PDim,RE*Pdim);
% Operador de Cruce
% Entradas
%
% Salidas
%
Offspring=PoolMating;
save('PoolMating','PoolMating')
switch OPC
    
    case 1 % CRUZAMIENTO 1 (un punto en 1/2 del cromosoma)
    
for i=1:2:CDim % para los impares 
    
    % Hijo 1 
    % Prioridad
    Offspring(i).Pr(:,1:CA/2)=PoolMating(i).Pr(:,1:CA/2); 
    Offspring(i).Pr(:,CA/2+1:end)=PoolMating(i+1).Pr(:,CA/2+1:end); 
      
    % Asignaci?n de recurso
    Offspring(i).RA(:,1:CA/2)=PoolMating(i).RA(:,1:CA/2); 
    Offspring(i).RA(:,CA/2+1:end)=PoolMating(i+1).RA(:,CA/2+1:end);
      
    % Hijo 2
    % Prioridad
    Offspring(i+1).Pr(:,CA/2+1:end)=PoolMating(i+1).Pr(:,CA/2+1:end);
    Offspring(i+1).Pr(:,1:CA/2)=PoolMating(i).Pr(:,1:CA/2);
    
    % Asignacion
    Offspring(i+1).RA(:,CA/2+1:end)=PoolMating(i+1).RA(:,CA/2+1:end);
    Offspring(i+1).RA(:,1:CA/2)=PoolMating(i).RA(:,1:CA/2);
    
end
    case 2 % CRUZAMIENTO 2 (un punto de corte aleatorio)
   
        pc= ceil(2 + (CA-2).*rand);
        
 for i=1:2:CDim % para los impares 
    
    % Hijo 1 
    % Prioridad
    Offspring(i).Pr(:,1:pc)=PoolMating(i).Pr(:,1:pc); 
    Offspring(i).Pr(:,pc+1:end)=PoolMating(i+1).Pr(:,pc+1:end); 
      
    % Asignaci?n de recurso
    Offspring(i).RA(:,1:pc)=PoolMating(i).RA(:,1:pc); 
    Offspring(i).RA(:,pc+1:end)=PoolMating(i+1).RA(:,pc+1:end);
      
    % Hijo 2
    % Prioridad
    Offspring(i+1).Pr(:,pc+1:end)=PoolMating(i+1).Pr(:,pc+1:end);
    Offspring(i+1).Pr(:,1:pc)=PoolMating(i).Pr(:,1:pc);
    
    % Asignacion
    Offspring(i+1).RA(:,pc+1:end)=PoolMating(i+1).RA(:,pc+1:end);
    Offspring(i+1).RA(:,1:pc)=PoolMating(i).RA(:,1:pc);
end
 
    case 3 % CRUZAMIENTO 2 (uniforme, con umbral adaptativo)
        
        IdFit=zeros(CDim,2);

        IdFit(1:CDim,1)=1:CDim;
        
        
        for i=1:CDim

        IdFit(i,2)=PoolMating(i).FO;
        
        end
         % IdFit
        IdFit=sortrows(IdFit,2);

        BestP=PoolMating(IdFit(IdFit(1:EDim,1))); 
        
        
        Cind=1:CDim;

        for i=1:(CDim-EDim)
        
            Cind=nonzeros(Cind);

            p1=Cind(unidrnd(numel(Cind)));
            
            p2=Cind(unidrnd(numel(Cind)));

            if p1==p2
                
                if p1==CDim
                
                 p2=p2-1;
                 
                else
                    
                 p2=p2+1;
                   
                end
                
               Cind(p1)=0;
               
               Cind(p2)=0;
               
            else
                
               Cind(p1)=0;
               
               Cind(p2)=0;
             
            end
                
            FOp1=PoolMating(p1).FO;
            
            FOp2=PoolMating(p2).FO;
            
            if FOp1<FOp2 % Si padre 1 es m?s apto que padre 2
               
                R=FOp1/FOp2;
                                               
                if R>=0.5 && R <= 1
                    
                    Alpha=0.5+(1-R)/5;
                    
                elseif R < 0.5
                    
                   Alpha=1-R;
                  
                end
                
              for j=1:NP
                
                  for k=2:CA-1
                     
                      Prob=rand;
                      
                      if Prob < Alpha
                          
                      Offspring(i).Pr(j,k)=PoolMating(p1).Pr(j,k); 
                      
                      Offspring(i).RA(j,k)=PoolMating(p1).RA(j,k); 
                      
                      else
                          
                      Offspring(i).Pr(j,k)=PoolMating(p2).Pr(j,k); 
                      
                      Offspring(i).RA(j,k)=PoolMating(p2).RA(j,k); 

                      end
                  end

              end

            else % Si padre 2 es m?s apto que padre 1
             
                R=FOp2/FOp1;
                
              if R>=0.5 && R <=1
                    
                    Alpha=0.5+(1-R)/5;
                    
               elseif R < 0.5
                    
                   Alpha=1-R;
                  
              end
              
              for j=1:NP
                             
                  for k=2:CA-1
                     
                      Prob=rand;
                      
                      if Prob < Alpha
                          
                      Offspring(i).Pr(j,k)=PoolMating(p2).Pr(j,k); 
                      
                      Offspring(i).RA(j,k)=PoolMating(p2).RA(j,k); 
                      
                      else
                          
                      Offspring(i).Pr(j,k)=PoolMating(p1).Pr(j,k); 
                      
                      Offspring(i).RA(j,k)=PoolMating(p1).RA(j,k); 

                      end
                  end

               end
  
            end

        end
   
     Offspring(((CDim-EDim)+1):end)=BestP;

    case 4  % CRUZAMIENTO 4 (cruce por rasgos de los padres)
      
     for i=1:2:CDim % para los impares 
    
      % Hijo 1 
      Offspring(i).Pr(:,1:CA)=PoolMating(i).Pr(:,1:CA); 
      Offspring(i).RA(:,1:CA)=PoolMating(i+1).RA(:,1:CA); 
      
      % Hijo 2
      Offspring(i+1).Pr(:,1:CA)=PoolMating(i+1).Pr(:,1:CA);
      Offspring(i+1).RA(:,1:CA)=PoolMating(i).RA(:,1:CA);    
     end
     
     
       case 5  % CRUZAMIENTO 4 (cruce por rasgos de los padres mixto)
      
     for i=1:2:CDim % para los impares 
       
       for j=1:2:NP  
      % Hijo 1 
      Offspring(i).Pr(j,1:CA)=PoolMating(i).Pr(j,1:CA); 
      Offspring(i).RA(j,1:CA)=PoolMating(i+1).RA(j,1:CA); 
      
      % Hijo 2
      Offspring(i+1).Pr(j,1:CA)=PoolMating(i+1).Pr(j,1:CA);
      Offspring(i+1).RA(j,1:CA)=PoolMating(i).RA(j,1:CA);  
      
      
      if j+1< NP
        % Hijo 1 
      Offspring(i).Pr(j+1,1:CA)=PoolMating(i).Pr(j+1,1:CA); 
      Offspring(i).RA(j+1,1:CA)=PoolMating(i+1).RA(j+1,1:CA); 
      
      % Hijo 2
      Offspring(i+1).Pr(j+1,1:CA)=PoolMating(i+1).Pr(j+1,1:CA);
      Offspring(i+1).RA(j+1,1:CA)=PoolMating(i).RA(j+1,1:CA);  
          
      end
      
       end
     end

    case 6 % CRUZAMIENTO 6 (cruce por puntos mixtos) 
     
     for i=1:2:CDim % para los impares 
       for j=1:2:NP
           
    % Hijo 1 
    % Prioridad
    Offspring(i).Pr(j,1:CA/2)=PoolMating(i).Pr(j,1:CA/2); 
    Offspring(i).Pr(j,CA/2+1:end)=PoolMating(i+1).Pr(j,CA/2+1:end); 
      % Hijo 2
    % Prioridad
    Offspring(i+1).Pr(j,CA/2+1:end)=PoolMating(i+1).Pr(j,CA/2+1:end);
    Offspring(i+1).Pr(j,1:CA/2)=PoolMating(i).Pr(j,1:CA/2);
    
       % Hijo 1 
    % Asig Recurso
    Offspring(i).Pr(j,1:CA/2)=PoolMating(i).Pr(j,1:CA/2); 
    Offspring(i).Pr(j,CA/2+1:end)=PoolMating(i+1).Pr(j,CA/2+1:end); 
      % Hijo 2
    % A Recurso
    Offspring(i+1).Pr(j,CA/2+1:end)=PoolMating(i+1).Pr(j,CA/2+1:end);
    Offspring(i+1).Pr(j,1:CA/2)=PoolMating(i).Pr(j,1:CA/2);

    if j+1<NP
       
    %Asignaci?n
    Offspring(i).Pr(j+1,1:CA/2)=PoolMating(i+1).Pr(j+1,1:CA/2); 
    Offspring(i).Pr(j+1,CA/2+1:end)=PoolMating(i).Pr(j+1,CA/2+1:end); 
    
    Offspring(i+1).Pr(j+1,CA/2+1:end)=PoolMating(i).Pr(j+1,CA/2+1:end);
    Offspring(i+1).Pr(j+1,1:CA/2)=PoolMating(i+1).Pr(j+1,1:CA/2);
        
    % Asignaci?n
    Offspring(i).RA(j+1,1:CA/2)=PoolMating(i+1).RA(j+1,1:CA/2); 
    Offspring(i).RA(j+1,CA/2+1:end)=PoolMating(i).RA(j+1,CA/2+1:end); 
    
    Offspring(i+1).RA(j+1,CA/2+1:end)=PoolMating(i).RA(j+1,CA/2+1:end);
    Offspring(i+1).RA(j+1,1:CA/2)=PoolMating(i+1).RA(j+1,1:CA/2);

    end
      
       end
    
     end
        
end

end

function [MutInd]=Mutation(DataProject,RRHH,TA,CP,CAType,MutInd,NP,CA,MDim,MutOp)
                         %(DataProject,RRHH,TA,CP,CAType,MutInd,NP,CA,RM*PDim,MutOP)
                    
% Operador de Mutaci??n 
% Entradas
%
% Salidas
AUX=zeros(NP,CA,4);

switch MutOp
    
    case 1
        
        for i=1:MDim
            n=unidrnd(5);
            AUX(:,:,1)=MutInd(i).TI(:,:);
            AUX(:,:,2)=MutInd(i).TF(:,:);
            AUX(:,:,3)=MutInd(i).RA(:,:);
            AUX(:,:,4)=MutInd(i).Pr(:,:);
            [AUX]=Neighborhood(DataProject,AUX,NP,n);
            MutInd(i).TI(:,:)=AUX(:,:,1);
            MutInd(i).TF(:,:)=AUX(:,:,2);
            MutInd(i).RA(:,:)=AUX(:,:,3);
            MutInd(i).Pr(:,:)=AUX(:,:,4);
        end



    case 2 
        
      [MutInd,x,y,z] = PobCreator(DataProject,RRHH,TA,MDim,CP,CAType);   
   
end
end

function [Ind] = SSGSGA(Ind,DataProject,RRHH,TA,CP,NP,CA)
% SSGSGA Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche.Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
%  Descripci??n: Heuristica para la factiblilidad basada 
% en: % Serial Scheme Generation Sechdeule SSGS  
% Dise?ada para  n proyectos n=2,3,... con la misma cantidad de actividades
% con tiempos de arribo Start=[T1,T2,T3,T4,...,Tn]
% los proyectos son tomados con dos actividades ficticias
% 

% DEFINICI??N DE VARIABLES
Schedule=cat(3,zeros(NP,CA),zeros(NP,CA),Ind.RA,Ind.Pr);

PS=zeros(1,NP); % Variable que guarda los proyectos que han comenzado PS(i)={0,1} 

PF=PS;      % Variable que guarda los proyectos que han finalizado PF(i)={0,1} 

APP=zeros(NP,CA);  % Variable que guarda las actividades que se han planificado APP(i,j)={0,1} 

% INICIALIZACI??N

PS(1)=1; %  Indica que el proyecto 1 comienza procesandose

APP(1,1)=1; %  Se planifica la primera actividad del primer proyecto

ARP=APP; % Conjuntos de actividades que pueden ser planificadas se inicializa con  la actividad 1 del proyecto 1  
 
Stop=0; %  Se inicializa en cero la condici??n de parada

ind=2:NP; % indice de los proyectos 

ind2=1:NP;

v=1; % Contador de los proyectos que inician

while Stop < NP % Mientras la cantidad de proyectos finalizado sea menor que la actidad de proyectos 
    
if v~=NP 

    MK=max(max(Schedule(:,2:end,2)));
            
for k=nonzeros(ind)' % Para cada uno de los Proyectos
    
    if TA(k)<=MK && TA(k)~=0 % Si el mayor si el Makespan parcial MK es mayor o igual a la fecha de arribo deque la fe
          
       PS(k)=1;
       
       Schedule(k,1,1)=TA(k); % se inician los proyectos que aun no han comenzado
       
       Schedule(k,1,2)=TA(k);
              
       v=v+1;
       
       ind(k-1)=0;
       
       ARP(k,1)=1;
       
       APP(k,1)=1;
                    
     end
       
end
    
end
         
 PR=find(PS==1);% Proyectos que han comenzado 
    
        CPS=length(PR); % Cantidad de proyectos iniciados
   
     if CPS==1 % Si la cantidad de proyectos iniciados es 1 
         
         OrdenP=PR; % Orden es PR
         
    else % en caso contrario
        
         OrdenP=zeros(1,length(PR));   
         
         MKP=OrdenP;
                         
         for l=PR
           
             MKP(l)=(CP(l)-(max(Schedule(l,:,2))-TA(l)))/CP(l);
                                                                       
         end
         
         aux=MKP(PR)';
                
         ORD=sortrows([aux,PR'],1);
         
        OrdenP=ORD(:,2)';
         
     end  

    for i=OrdenP % Para proyectos Iniciado segun OrdenP
                              
     [Schedule,APP,ARP] = SchedulingFeasible(Schedule,DataProject,RRHH,i,nonzeros(ARP(i,:))',ARP,APP);                                                         
                                              
    end
        
  for r=nonzeros(ind2)'
    
    if Schedule(r,CA,2)~=0 && sum(APP(r,:))==CA
        
        PF(r)=1;
        
        PS(r)=0;   
        
        ind2(r)=0;
        
        ARP(r,:)=0;
             
    end
    
  end
  Stop=sum(PF);
end  
Ind.FO=mean(((Schedule(:,CA,2)'-(CP+TA))));  %mean((Schedule(:,CA,2)'./CP)-1)); %mean(Schedule(:,CA,2));
end

function [Schedule,APP,ARP] = SchedulingFeasible(Schedule,DataProject,RRHH,i,j,ARP,APP)
%  SchedulingFeasible Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
%  Planifica las activitades de ActReadyPlan sin violar las restricciones de
% Precedencia y de Recursos.
% ERNTRADAS
% Schedule -> soluci??n 
% DataProyect -> Datos de los Proyectos 
%ActReadyPlan-> Conjunto de actividades que est??n lista para ser planificadas 
% RRHH-> Datos sobre los recursos humanos del proyecto
% APP(i) -> Conjunto de Actividades planificadas por Proyecto.
% DataProject -> Informaci???n sobre los proyectos.
% SALIDAS
% Schedule-> soluci??n 
% APP(i) -> Actividades planificadas por Proyecto.

% INICIALIZACI??N

 ARP(i,:)=0;
 
 AR=nonzeros(DataProject(nonzeros(j),4:6,i))'; % Conjunto de actividades  sucesoras inmediatas de la actividad j    
 
 ModeProc=DataProject(:,7:9,i);

 aux=length(AR);
 
 AxW=[Schedule(i,AR,4);AR];
 
   OrdenA=sortrows(AxW',1);
 
   Act=OrdenA(:,2)';
 
 for k=1:aux
     
      [AA,p]=find(DataProject(:,[4 5 6],i)==Act(k)); %#ok<NASGU>
          
      CAP=sum(APP(i,AA'));
      
      if CAP==length(AA) && APP(i,Act(k))==0 % si est?? lista y nunca ha sido planificada
          
          ARP(i,Act(k))=Act(k);
            
           APP(i,Act(k))=1;
          
          TipoR=DataProject(Act(k),2,i); % Definir el tipo  de recurso que requiere la actividad j del proyecto i
           
          if TipoR==0 % Si la actividad j require un tipo de recurso=0 => ??ltima actividad del proyecto 
       
            Schedule(i,Act(k),1)=max(Schedule(i,AA',2)); % Tiempo de inicio de la actividad j=CA
               
            Schedule(i,Act(k),2)=Schedule(i,Act(k),1);         %  Tiempo de finalizaci??n de la actidad j=CA 
                      
          else
               
             RAsig=Schedule(i,Act(k),3); % Asignaci??n del Recurso
            
             TFS=Schedule(:,:,2);
                                       
             TFAmr=TFS(Schedule(:,:,3)==RAsig); % Buscar si el recurso RAsig ha sido asignado anteriormente
             
             if isempty(TFAmr)==1 && RAsig~=0 % Si el recurso asignado a la actividad j nunca ha sido asignado
                 
               Schedule(i,Act(k),1)=max(Schedule(i,AA',2)); % Calcula tiempo de inicio de la actividad j del proyecto i
               
                             
               Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalizaci??n de la actividad j del proyecto i
             
             else
                if RAsig~=0
                                 
               Schedule(i,Act(k),1)=max(max(TFAmr),max(Schedule(i,AA',2))); % Calcula tiempo de inicio de la actividad j del proyecto i
             
               Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalizaci??n de la actividad j del proyecto i
                end
             end
          end          
      end
 end
end

function [Schedule]=Neighborhood(DataProject,Schedule,NP,n)
% NeighborhoodR function  
% Autor: Lic Bolivar E. Medrano Broche. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n: operador para generar soluciones vecinas basadas en
% intercambio de asignaciones de recursos o de actividades actividades cualquiera.
% ENTRADAS
%
% 
% IMPLEMENTACI??N 

switch n
    case 1
       for i=1:NP % Para cada uno de los proyectos
   
    % Seleccionar de la actividades
     
       TypeA=unidrnd(4);
       [CAType,u]=find(DataProject(:,2,i)==TypeA);
       aux=unidrnd(numel(CAType));
       A1=CAType(aux); CAType(aux)=0;
       CAType=nonzeros(CAType);
       A2=CAType(unidrnd(numel(CAType)));
      
   % Intercambio de Recursos   
   
      ExchRA1=Schedule(i,A1,3);
      ExchRA2=Schedule(i,A2,3);
      Schedule(i,A1,3)=ExchRA2;
      Schedule(i,A2,3)=ExchRA1;

  % Intercambio de las Priodades
   
      ExchPA1=Schedule(i,A1,4);
      ExchPA2=Schedule(i,A2,4);
      Schedule(i,A1,4)=ExchPA2;
      Schedule(i,A2,4)=ExchPA1;

       end 
    case 2
        Pid=1:NP; ind=1:4;
        if NP<4
            aux=unidrnd(4);
            ind(aux)=0;
            ind=nonzeros(ind);
        end
           k=0;
    for i=ind' % Para cada uno de los  proyectos
       k=k+1;
        % Selecci??n de los proyectos P1 y P2 
       P1=unidrnd(NP);
      if k>1
       Pid(P1)=0;
      end

       Pid=nonzeros(Pid);
       P2=Pid(unidrnd(numel(Pid)));

      % Seleccionar las actividades A1, A2
       [CAType1,u]=find(DataProject(:,2,P1)==i);
       [CAType2,u]=find(DataProject(:,2,P2)==i);
       aux1=unidrnd(numel(CAType1)); aux2=unidrnd(numel(CAType2));
       A1=CAType1(aux1); A2=CAType2(aux2); 

     % Intecambio de recursos 
      ExchRA1=Schedule(P1,A1,3);
      ExchRA2=Schedule(P2,A2,3);

     Schedule(P1,A1,3)=ExchRA2;
     Schedule(P2,A2,3)=ExchRA1;

    % Intecambio de Prioridad
     ExchPA1=Schedule(P1,A1,4);
     ExchPA2=Schedule(P2,A2,4);

     Schedule(P1,A1,4)=ExchPA2;
     Schedule(P2,A2,4)=ExchPA1;
   end
    case 3
  Pid=1:NP;
 P1=unidrnd(NP);
 Pid(P1)=0;
 
 Pid=nonzeros(Pid);

 P2=Pid(unidrnd(numel(Pid)));
    
 % Seleccionar de la actividades 
    TypeA=unidrnd(4);
    [CAType1,u]=find(DataProject(:,2,P1)==TypeA);
    [CAType2,u]=find(DataProject(:,2,P2)==TypeA);
    aux1=unidrnd(numel(CAType1));A1=CAType1(aux1); 
    aux2=unidrnd(numel(CAType2));A2=CAType2(aux2); 

% Intercambio de Recursos      
  ExchRA1=Schedule(P1,A1,3);
  ExchRA2=Schedule(P2,A2,3);
  Schedule(P1,A1,3)=ExchRA2;
  Schedule(P2,A2,3)=ExchRA1;

% Intercambio de las Priodades
  ExchPA1=Schedule(P1,A1,4);
  ExchPA2=Schedule(P2,A2,4);
  Schedule(P1,A1,4)=ExchPA2;
  Schedule(P2,A2,4)=ExchPA1;
    case 4

% Para cada uno de los  proyectos P1 y P2 

Type=unidrnd(4);

P=unidrnd(NP);

% Seleccionar las actividades A1, A2

[CAType,u]=find(DataProject(:,2,P)==Type);

aux=unidrnd(numel(CAType)); 

A1=CAType(aux); 
CAType(aux)=0; 
CAType=nonzeros(CAType); 
A2=CAType(unidrnd(numel(CAType))); 

% Intecambio de recursos 
ExchRA1=Schedule(P,A1,3);
ExchRA2=Schedule(P,A2,3);

Schedule(P,A1,3)=ExchRA2;
Schedule(P,A2,3)=ExchRA1;

% Intecambio de Prioridad
ExchPA1=Schedule(P,A1,4);
ExchPA2=Schedule(P,A2,4);

Schedule(P,A1,4)=ExchPA2;
Schedule(P,A2,4)=ExchPA1;
    case 5
    % Para  los proyectos P1 y P2
 Pid=1:NP;
 P1=unidrnd(NP);
 Pid(P1)=0;
 Pid=nonzeros(Pid);
 
 P2=Pid(unidrnd(numel(Pid)));
    
 % Seleccionar de la actividades 
   
 TypeA=unidrnd(4);
 
    % Selecci??n de las actividades
    
    if P1==P2 % si los proyectos son iguales
        [CAType,u]=find(DataProject(:,2,P1)==TypeA);
        aux=unidrnd(numel(CAType)); 
        A1=CAType(aux); CAType(aux)=0; CAType=nonzeros(CAType);
        aux=unidrnd(numel(nonozeros(CAType)));
        A2=CAType(aux);
    else
       [CAType1,u]=find(DataProject(:,2,P1)==TypeA);
       [CAType2,u]=find(DataProject(:,2,P2)==TypeA);
       aux1=unidrnd(numel(CAType1)); A1=CAType1(aux1); 
       aux2=unidrnd(numel(CAType2)); A2=CAType2(aux2);   
    end
   

    
% Intercambio de Recursos      
  ExchRA1=Schedule(P1,A1,3);
  ExchRA2=Schedule(P2,A2,3);
  Schedule(P1,A1,3)=ExchRA2;
  Schedule(P2,A2,3)=ExchRA1;

% Intercambio de las Priodades
  ExchPA1=Schedule(P1,A1,4);
  ExchPA2=Schedule(P2,A2,4);
  Schedule(P1,A1,4)=ExchPA2;
  Schedule(P2,A2,4)=ExchPA1;    
        
end
end

