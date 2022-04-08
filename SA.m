function [M,BestSchedule,FO] = SA(DataProject,RRHH,TA,CP,CAType,CI,Alpha,nrep,n,vns)
% Simulate Annealing function 
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n:
% Entradas
% DataProject
%
% CI -> Cantidad de Inteciones
% To -> 
% Salidas
% INICIALIZACION
%tic
    
    [SolInit,FVo,NP,CA]=SSGS(DataProject,RRHH,TA,CP,CAType);

    [lb,ub]=LbUb(DataProject,NP); % Estimacion de cotas 

    BestSchedule=SolInit;
    
    Best=FVo;

    FO=zeros(1,CI+10);

    FO(1)=Best;

    To=((lb-ub)/log((Alpha)));    % valor inicial de T

    T=To;

    c=0; k=1;  M=10^3;    
    ch=0;

 while (c < CI) 
     
      i=1;
     
      k=k+1;
      
      T1=T;
    
      while  (i<nrep)  %(i < nrep)
     
          i=i+1; % contador del n?mero de repeticiones
          
          if ch==25
              if n==5
                  n=1;
                  fprintf(1,'N activo %1.0f \n',n)
              else
                  n=n+1;
                  fprintf(1,'N activo %1.0f \n',n)
              end
             
          end
     
          [Schedule]=Neighborhood(DataProject,SolInit,NP,n); % generador de vecindades
               
          [Schedule,FV]=SSGSSA(Schedule,DataProject,RRHH,TA,CP,NP,CA);  % nueva soluci?n factible de vecindades
                                  
          Deltha=FV-FVo; % diferencia entre soluci?n n-1 y n   
     
     if Deltha>0 %  calcular la probabilidad de aceptaci?n
         
     Aceptance =exp(-(Deltha)/T1); % probabilidad de aceptaci?n

     else
       
        Aceptance=0;              % probabilidad de aceptaci?n
     
     end
      
     if Deltha < 0 || Aceptance > rand  % aceptar soluci?n actual

         %fprintf(1,' Acepto %4.0f \n', i)
                 
         SolInit=Schedule; % soluci?n actual
        
         Best=FV; % FO aceptado
          
         if  Best < M % actualizar la mejor soluci?n encontrada

             M=Best;  % FO de la mejor FO encontrada  
           
             BestSchedule=struct('TI',SolInit(:,:,1),'TF',SolInit(:,:,2),'RA',SolInit(:,:,3),'FO',M);
            if vns==2
             ch=0;
            end
         else
             if vns==2
             ch=ch+1;
             end
         end
          
     else % acetpaci?n
                 
     end % if aceptaci?n 

      FVo=Best;
      
      FO(c+i)=Best;
      
      end % while n?mero de repeticiones
      
      c=c+i;

      Alpha=((0.87-0.92)/100^2)*(100-k)^2 + 0.92;

      T=T*Alpha;%/(k*Alpha);%  %decrecimiento de la temperatura log(k+1);%
       
      fprintf(1,'IT %5.0f Alpha %0.10f TEMP %4.10f\n',k,Alpha,T)
 end
%toc
end

function [Schedule,FO,NP,CA]= SSGS(DataProject,RRHH,TA,CP,CAType)
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

p=unidrnd(9);

w1=p/10; aux=10-p;

if aux==1
    
  w2=0.05; w3=0.05;
  
else
   
    aux1=unidrnd(aux-1);  w2=aux1/10;  w3=(aux-aux1)/10;
end

DMR=[w1,w2,w3]; % DMR=[.33 .33 .33]


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

%DMR=[.5 .3 .2];

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

function [Schedule]=Neighborhood(DataProject,Schedule,NP,n)
% NeighborhoodR function  
% Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
% Facultad III. Universidad de las Ciencias Inform??ticas
% Descripci??n: Realiza una intercambio de recursos del mismo tipo entre dos
% actividades cualquiera.
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

function [Schedule,FO] = SSGSSA(Schedule,DataProject,RRHH,TA,CP,NP,CA)
% SSGSnew Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
%  Descripci??n: Heuristica para la generaci??n de soluciones iniciales factibles basada
% en: % Serial Scheme Generation Sechdeule SSGS  
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
%  TA -> Tiempos de arribo de los proyectos  TA=[TA1,TA2,...,TAn]
% SALIDAS
%  Schedule -> Soluci??n  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  TI -> Tiempos de Inicio
%  TF -> Tiempos de Finalizaci??n 
%  RA -> Recurso Asignado

% DEFINICI??N DE VARIABLES

Schedule(:,:,1)=0;

Schedule(:,:,2)=0;

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
                              
     [Schedule,APP,ARP] = SchedulingFeasibleSA(Schedule,DataProject,RRHH,i,nonzeros(ARP(i,:))',ARP,APP);                                                         
                                              
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
end

function [Schedule,APP,ARP] = SchedulingFeasibleSA(Schedule,DataProject,RRHH,i,j,ARP,APP)
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

function [lb,ub]=LbUb(NET,NP)
aux=zeros(1,NP);
for m=[1,3]
    for i=1:NP
        [aux(i)] = CPM(NET,i,m);
    end
    if m==1
        lb=mean(aux);
    elseif m==3
        ub=mean(aux); 
    end
end

end


function [FO] = CPM(DataProject,i,M)
% CPM  Summary of this function 
%  Autor: Lic Bolivar E. Medrano Broche. CEGEL & Dpto. Ciencias B??sicas.
%  Facultad III. Universidad de las Ciencias Inform??ticas
%  Descripci??n: Calcula la Ruta Cr??tica de un proyecto dado um modo de
%  procesamiento
% 
% ENTRADA
%  DataProject un es un arreglo multidimencional  (P0,P1,P2,P3,...,) 
%  Donde  Pi -> Matriz (CA cantidad de actividades) x (Atributos del Proyecto i -> 9 columnas )  
%  Atributos -> 
%  Id de actividad (ID)-> col1
%  Tipo de recursos (TR)-> col2	
%  No sucesores imediatos (NSI)-> col3 
%  Sucesores inmediatos (SI1 SI2 SI3)-> col4-col6	
%  M: Modos de Procesamiento (M1	M2	M3)-> col7-col9
%  
% SALIDAS
%  Schedule -> Soluci??n  (NPxCA)  Filas Proyectos - Columnas Actividades (Matriz con los tiempos de inicio de cada actividad de cada proyecto ) 
%  TI -> Tiempos de Inicio
%  TF -> Tiempos de Finalizaci??n 

% DEFINICI??N DE VARIABLES

[CA,H,NP]=size(DataProject); % Datos de los Proyectos (CA -> cantidad de actividades, NP-> cantidad de Proyectos)

NP=1;

Schedule=cat(3,zeros(NP,CA),zeros(NP,CA)); % Estructura de la soluci??n

APP=zeros(NP,CA);  % Variable que guarda las actividades que se han planificado APP(i,j)={0,1} 

% INICIALIZACI??N

APP(1,1)=1; %  Se planifica la primera actividad del primer proyecto

ARP=APP; % Conjuntos de actividades que pueden ser planificadas se inicializa con  la actividad 1 del proyecto 1  
 
Stop=0; %  Se inicializa en cero la condici??n de parada

while Stop < NP % Mientras la cantidad de proyectos finalizado sea menor que la actidad de proyectos 
    
  
        for j=nonzeros(ARP(1,:))
                               
            if j ~= CA
                                          
              SI=nonzeros(DataProject(j,4:6,i));   % Actividades  que tienen antecesora la actividad j                        
                         
             [ActReadyP]=CheckActReady(DataProject(:,[4 5 6],i),1,SI',APP); % Chequea si las actividades ANT estan listas para planificarse 
                                           
             [Schedule,APP] =SchedulingFeasibleCPM(Schedule,DataProject,ActReadyP,i,CA,APP,M); % inserta en la soluci??n el conjunto de actividades (ActReadyP) cumpliendo las restricciones de tiempo y recurso             
                                             
             ARP(1,:)=0;
       
             ARP(1,ActReadyP)=ActReadyP;
             
            else                                        
                  [Schedule,APP] =SchedulingFeasibleCPM(Schedule,DataProject,j,i,CA,APP);
                  
                  ARP(1,:)=0;                  
            end
        end    
        
    if Schedule(1,CA,2)~=0 && sum(APP(1,:))==CA
        
       Stop=1; 
             
    end
end  

FO=max(Schedule(1,:,2));
%toc
end

function [ActReadyP] =CheckActReady(SS,pro,act,APP)
                                 
% Sumary SelectActReady 
%   Selecciona las actividades que est???n listas para planificarse. 
% Entradas
%- P-> Proyecto i
%- APP-> Actividades planificadas por proyecto.
%-
% Salida
% ActReadyPlan: Lista de actividades de cada proyecto lista para planificarse
% A variable Auxiliar 
% Implementaci??n

ActReady=zeros(1,length(act));

for c=1:length(act);
    
    [AA,p]=find(SS==act(c));
    
    CAP=sum(APP(pro,AA')); 
        
    if CAP==length(AA)
        
        ActReady(c)=act(c);    
    end

end
ActReadyP=intersect(nonzeros(ActReady)',nonzeros(ActReady)');
end

function [Schedule,APP] = SchedulingFeasibleCPM(Schedule,DataProject,ActReady,i,CA,APP,M)
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
% APP(i) -> Lista de Actividades planificadas por Proyecto.

% INICIALIZACI??N

m=[7 8 9];
      
ModeProc=DataProject(:,7:9,i); % Modos de Procesamiento del as actividades del proyecto i    

for j=ActReady % Para el orden establecido por su prioridad        
    
if APP(1,j)==0 % Si la actividad i del proyecto j nunca ha sido planificada 
    
     APP(1,j)=1;  
    
    if j~=CA
        
        switch M
            case 1
        
                [ACAnt,u]=find(DataProject(:,4:6,i)==j); % Busca el conjunto de actividades antecesoras de las la actividad j
                       
                Schedule(1,j,1)=max(Schedule(1,ACAnt',2)); % Tiempo de inicio de la actividad j=CA
               
                Schedule(1,j,2)=Schedule(1,j,1) + ModeProc(j,M); %  Tiempo de finalizaci??n de la actidad j=CA 
            case 2
                
                [ACAnt,u]=find(DataProject(:,4:6,i)==j); % Busca el conjunto de actividades antecesoras de las la actividad j
                       
                Schedule(1,j,1)=max(Schedule(1,ACAnt',2)); % Tiempo de inicio de la actividad j=CA
               
                Schedule(1,j,2)=Schedule(1,j,1) + ModeProc(j,M); %  Tiempo de finalizaci??n de la actidad j=CA 
           
            case 3 
                
                [ACAnt,u]=find(DataProject(:,4:6,i)==j); % Busca el conjunto de actividades antecesoras de las la actividad j
                       
                Schedule(1,j,1)=max(Schedule(1,ACAnt',2)); % Tiempo de inicio de la actividad j=CA
               
                Schedule(1,j,2)=Schedule(1,j,1) + ModeProc(j,M); %  Tiempo de finalizaci??n de la actidad j=CA 
                
            case 4
                
                [ACAnt,u]=find(DataProject(:,4:6,i)==j); % Busca el conjunto de actividades antecesoras de las la actividad j
                       
                Schedule(1,j,1)=max(Schedule(1,ACAnt',2)); % Tiempo de inicio de la actividad j=CA
               
                Schedule(1,j,2)=Schedule(1,j,1) + mean(ModeProc(j,:)); %  Tiempo de finalizaci??n de la actidad j=CA 
        end     
    
    else
        
     Schedule(1,j,1)=max(Schedule(1,:,2));  % Calcula tiempo de inicio de la actividad j=CA del proyecto i
    
     Schedule(1,j,2)=Schedule(1,ActReady,1); % Calcula tiempo de finalizaci??n de la actividad j=CA del proyecto i
              
    end 
        
 end
end   
end

