function [GBFO,GBest]=PSO_AGb(DataProject,RRHH,TA,CP,CAType,SS,CI,SWP) 
% DPSO  function -> Discrete Particle Swarm Optimization (DPSO)
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n: Algoritmo DPSO para la soluci??n del problema SDPSP
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
% SwarmSize -> Cantidad de part??culas 
% Salidas
% Swarm-> Martiz que contiene al Enjambre.

% INICIALIZACI??N 
 
%all=zeros(SS,CI);
%dis1=all; dis2=all;
% Inicializaci??n de los par??metros

C1=SWP(1); % Factor Personal

C2=SWP(2);  % Factor Social

%FO=zeros(1,CI);

% Inicializaci??n del Enjambre

[Swarm,SwarmPB,GBest,GB,NP,CA] = SwarmCreatorNew(DataProject,RRHH,TA,SS,CP,CAType); % Crea el enjambre O(4n???????)

for i=1:CI % Para cada una de las iteraciones O(n)
    
    Alpha=((0.6-SWP(4))/CI^2)*(CI-i)^2 + SWP(4); % 0.01-0.6
    
  %  X=((0.7-SWP(3))/CI^2)*(CI-i)^2 + SWP(3); % 0.9-0.2
                        
     for j=1:SS % Para cada una de las part??culas O(n)      
         
        m=(Swarm(j).FO-GBest.FO)/(GBest.FO+Swarm(j).FO);
           
         X=0.85 + (SWP(3)-0.85)*((exp(m)-1)/(exp(m)+1));
      
          R1=C1*rand(1,CA); % Fuerzas aleatorias
    
          R2=C2*rand(1,CA); % Fuerzas aleatorias
                                                
          Swarm(j).Va(1:CA) = X*(Swarm(j).Va(1:CA)) + R1(1:CA).*(SwarmPB(j).Xa(1:CA) - Swarm(j).Xa(1:CA)) + R2(1:CA).*(GBest.Xa(1:CA) - Swarm(j).Xa(1:CA));
                
          Swarm(j).Xa(1:CA) = Swarm(j).Xa(1:CA) + Swarm(j).Va(1:CA);
                     
          [Swarm(j)]=SSGSSwarm(Swarm(j),GBest,SwarmPB(j),DataProject,RRHH,TA,NP,CA,CP,X,R1,R2,Alpha);

          FOPB=SwarmPB(j).FO; CurrPFO=Swarm(j).FO;
        
        if  FOPB > CurrPFO
 
           SwarmPB(j)=Swarm(j);
           
        end 
              
     end % For -> para cada part??cula 
    
     % Actulizaci??n de GBest
     
      GBFOcurr=GBest.FO;
    
       [GBFO,GB]=min([Swarm(1:SS).FO]); % GBFO Valor de la funci??n objetivo de la mejor part??cula; PGB particula con la mejor evalucion global

       aux=length(GB);

       if  aux~=1

        p=unidrnd(aux);
 
        GB=GB(p);
        
        GBFO=GBFO(p);
   
       end
     
     if GBFOcurr > GBFO
         
     GBest=Swarm(GB);
           
     end
    
  % FO(i)=GBest.FO;
   
   GBFO=GBest.FO;
   
end % For -> cada interaci??n 

end

function [Swarm,SwarmPB,GBest,GB,NP,CA] = SwarmCreatorNew(DataProject,RRHH,TA,SS,CP,CAType)
% SwarmCreator function -> Creador del enjambre
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
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
% SS->SwarmSize -> Cantidad de part??culas 
% SALIDAS
% Swarm-> Martiz que contiene al Enjambre.
% SwarmPB -> Mejor Posici??n que tiene el Enjambre.
% Gbest -> Mejor Part??cula del Enjambre.
% NP -> Cantidad de proyectos
% CA -> Cantidad de Actividades

% IMPLEMENTACI??N

Swarm=[];

% Generacion aleatoria de los pesos para la asignacion de recursos
 
for i=1:SS
    
    p=unidrnd(9);

w1=p/10; aux=10-p;

if aux==1
    
  w2=0.05; w3=0.05;
  
else
   
    aux1=unidrnd(aux-1);  w2=aux1/10;  w3=(aux-aux1)/10;
end

DMR=[w1,w2,w3]; % DMR=[.33 .33 .33]
    
   [Schedule,FV,NP,CA] = SSGS(DataProject,RRHH,TA,CP,CAType,DMR);
 
   Particle=struct('TI',Schedule(:,:,1),'TF',Schedule(:,:,2),'RA',Schedule(:,:,3),'Xa',Schedule(:,:,4),'Va',rand(NP,CA),'Vr',rand(NP,CA),'Xr',zeros(NP,CA),'FO',FV); % Estructura de una part??cula
   
   Swarm=cat(3,Swarm,Particle); % Crea enjambre
       
end
 
SwarmPB=Swarm; % Inicializa SwarmPB -> mejor posici??n personal inicial, como el enjambre inicial               

[GBFO,GB]=min([Swarm(1:SS).FO]); % GBFO Valor de la funci??n objetivo de la mejor part??cula; PGB particula con la mejor evalucion global

GBest=Swarm(GB);

aux=length(GB);

if  aux~=1

 p=unidrnd(aux);
 
 GB=GB(p);
   
end
aux=[0,-1];
for i=1:SS
    for k=1:NP
        for j=1:CA
            
            if i==GB
                
             Swarm(i).Xr(k,j)=(-1)^unidrnd(2); %#ok<AGROW>
                
            else
                        
            if Swarm(i).RA(k,j)==Swarm(GB).RA(k,j) && Swarm(i).RA(k,j)==SwarmPB(i).RA(k,j)
            
            Swarm(i).Xr(k,j)=(-1)^unidrnd(2); %#ok<AGROW>
            
            else 
            
            Swarm(i).Xr(k,j)=aux(unidrnd(2)); %#ok<AGROW>
                
            end

            end
            
        end
    end
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
%  W  -> Prioridad (MTS)

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

function[CurrParticle] = SSGSSwarm(CurrParticle,GBest,PBest,DataProject,RRHH,TA,NP,CA,CP,X,R1,R2,Alpha)
% SSGSwarmA Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Controla la factibilidad de cada part??cula del enjambante est?? basada en Serial Scheme Generation Sechdeule SSGS  
% 
% ENTRADA
%  CurrParticle-> contiene la part??cula que acaba de ser modificada
%  DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  RRHH -> Informaci??n sobre los recursos 
%  TA -> Tiempos de arribo de los proyectos  TA=[TA1,TA2,...,TAn]+
%  CP -> Ruta cr??tica
% SALIDAS
%  Schedule -> Soluci??n  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  CurrParticle -> Particula con factivilidad factible
%  FO -> Evaluaci??n de la particula
%  RA -> Recurso Asignado

% DEFINICI??N DE VARIABLES
%R1=C1*rand(1,CA); R2=C2*rand(1,CA);
Schedule=cat(3,zeros(NP,CA),zeros(NP,CA),CurrParticle.RA,CurrParticle.Xa,CurrParticle.Vr,CurrParticle.Xr);

Schedule(:,1,1)=TA';

Schedule(:,1,2)=TA';

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
                                       
     [Schedule,APP,ARP] = SchedulingFeasibleSwarmA(Schedule,PBest,GBest,DataProject,RRHH,i,ARP(i,:),ARP,APP,X,R1,R2,Alpha);
                                                  
    end
        
  for r=nonzeros(ind2)'
   
    if Schedule(r,CA-1,2)~=0 && sum(APP(r,:))==CA
        
        PF(r)=1;
        
        PS(r)=0;   
        
        ind2(r)=0;
        
        ARP(r,:)=0;
             
    end
    
  end

  Stop=sum(PF);
end  
CurrParticle.FO=mean(((max(Schedule(:,:,2),[],2)'-(CP+TA))));  %mean((Schedule(:,CA,2)'./CP)-1)); %mean(Schedule(:,CA,2));
CurrParticle.TI=Schedule(:,:,1);CurrParticle.TF=Schedule(:,:,2);
CurrParticle.RA=Schedule(:,:,3);CurrParticle.Vr=Schedule(:,:,5);
CurrParticle.Xr=Schedule(:,:,6);

end

function [Schedule,APP,ARP] = SchedulingFeasibleSwarmA(Schedule,PBest,GBest,DataProject,RRHH,i,j,ARP,APP,X,R1,R2,Alpha)
%  SchedulingFeasible Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
% Inserta en la soluci??n las activitades sin violar las restricciones de  Precedencia y de Recursos.
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

Acept=0;

 ARP(i,:)=0;
  
 AR=nonzeros(DataProject(nonzeros(j),4:6,i))'; % Conjunto de actividades  sucesoras inmediatas de la actividad j    
 
 ModeProc=DataProject(:,7:9,i);

 aux=length(AR);
 
 AxW=[Schedule(i,AR,4);AR];
 
   OrdenA=sortrows(AxW',1);
 
   Act=OrdenA(:,2)';

   for k=1:aux
     
      [AA,p]=find(DataProject(:,[4 5 6],i)==Act(k));
          
      CAP=sum(APP(i,AA'));
                   
      if CAP==length(AA) && APP(i,Act(k))==0 % si est?? lista y nunca ha sido planificada
          
          ARP(i,Act(k))=Act(k);
            
           APP(i,Act(k))=1;                   
           
           TipoR=DataProject(Act(k),2,i); % Definir el tipo  de recurso que requiere la actividad j del proyecto i
           
    
           
      if TipoR~=0  % Si la actividad j require un tipo de recurso=0 => ??ltima actividad del proyecto 
         
       Schedule(i,Act(k),5) = X*Schedule(i,Act(k),5) + R1(Act(k))*(-1-Schedule(i,Act(k),6)) + R2(Act(k))*(1-Schedule(i,Act(k),6));

       H=Schedule(i,Act(k),6)+Schedule(i,Act(k),5); 

       
       %fprintf(1,'H %2.4g\n',H)
       %fprintf(1,'Alpha %2.4g\n',Alpha)
       
        if H > Alpha
        %   disp('Gbest')
           Schedule(i,Act(k),6)=1;
           
           Schedule(i,Act(k),3)=GBest.RA(i,Act(k));
            
        elseif H < -Alpha 
         %   disp('Lbest')
           Schedule(i,Act(k),6)=-1;
          
           Schedule(i,Act(k),3)=PBest.RA(i,Act(k));
           
        else
          %  disp('Rand')
            Schedule(i,Act(k),6)=0;
            
             TypeR=RRHH(RRHH(:,2)==DataProject(Act(k),2,i),1);
            
            Schedule(i,Act(k),3)=TypeR(unidrnd(length(TypeR)));

        end      
                      
             RAsig=Schedule(i,Act(k),3); % Asignaci??n del Recurso
                      
             TFS=Schedule(:,:,2);
                                       
             TFAmr=TFS(Schedule(:,:,3)==RAsig); % Buscar si el recurso RAsig ha sido asignado anteriormente
             
             if isempty(TFAmr)==1 % Si el recurso asignado a la actividad j nunca ha sido asignado
                 
                 Schedule(i,Act(k),1)=max(Schedule(i,AA',2)); % Calcula tiempo de inicio de la actividad j del proyecto i

                Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalizaci??n de la actividad j del proyecto i
             
             else
                 
              Schedule(i,Act(k),1)=max(max(TFAmr),max(Schedule(i,AA',2))); % Calcula tiempo de inicio de la actividad j del proyecto i
             
               Schedule(i,Act(k),2)=Schedule(i,Act(k),1) + ModeProc(Act(k),RRHH(RAsig,3)); % Calcula tiempo de finalizaci??n de la actividad j del proyecto i
                             
             end
             
             
      end
      end
   end
 end
