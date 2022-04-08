tic

Inst=40; runs=2; alg=6;

Alpha=0.85;

BestFO=zeros(6,40);

BestFO(:,:)=1000;

Data=zeros(runs,alg,Inst);

for i=[34,35,37]%38:Inst % Para cada una de las instancias

    switch i
        case 1,load('3p30a_1.mat');disp('3p30a_1.mat');
        case 2,load('3p30a_2.mat');disp('3p30a_2.mat');
        case 3,load('3p30a_3.mat');disp('3p30a_3.mat');
        case 4,load('3p30a_4.mat');disp('3p30a_3.mat');
        case 5,load('3p30a_5.mat');disp('3p30a_5.mat');
        case 6,load('3p60a_1.mat');disp('3p60a_1.mat');
        case 7,load('3p60a_2.mat');disp('3p60a_2.mat');
        case 8,load('3p60a_3.mat');disp('3p60a_3.mat');
        case 9,load('3p60a_4.mat');disp('3p60a_4.mat');
        case 10,load('3p60a_5.mat');disp('3p60a_5.mat');
        case 11,load('3p90a_1.mat');disp('3p90a_1.mat');
        case 12,load('3p90a_2.mat');disp('3p90a_2.mat');
        case 13,load('3p90a_3.mat');disp('3p90a_3.mat');
        case 14,load('3p90a_4.mat');disp('3p90a_4.mat');
        case 15,load('3p90a_5.mat');disp('3p90a_5.mat');
        case 16,load('3p120a_1.mat');disp('3p120a_1.mat');
        case 17,load('3p120a_2.mat');disp('3p120a_2.mat');
        case 18,load('3p120a_3.mat');disp('3p120a_3.mat');
        case 19,load('3p120a_4.mat');disp('3p120a_4.mat');
        case 20,load('3p120a_5.mat');disp('3p120a_5.mat');
        case 21,load('5p30a_1.mat');disp('5p30a_1.mat');           
        case 22,load('5p30a_2.mat');disp('5p30a_2.mat');  
        case 23,load('5p30a_3.mat');disp('5p30a_3.mat');  
        case 24,load('5p30a_4.mat');disp('5p30a_4.mat');  
        case 25,load('5p30a_5.mat');disp('5p30a_5.mat');  
        case 26,load('5p60a_1.mat');disp('5p60a_1.mat');  
        case 27,load('5p60a_2.mat');disp('5p60a_2.mat');  
        case 28,load('5p60a_3.mat');disp('5p60a_3.mat');  
        case 29,load('5p60a_4.mat');disp('5p60a_4.mat');  
        case 30,load('5p60a_5.mat');disp('5p60a_5.mat');  
        case 31,load('5p90a_1.mat');disp('5p90a_1.mat');  
        case 32,load('5p90a_2.mat');disp('5p90a_2.mat');  
        case 33,load('5p90a_3.mat');disp('5p90a_3.mat');  
        case 34,load('5p90a_4.mat');disp('5p90a_4.mat'); 
        case 35,load('5p90a_5.mat');disp('5p90a_5.mat'); 
        case 36,load('5p120a_1.mat');disp('5p120a_1.mat'); 
        case 37,load('5p120a_2.mat');disp('5p120a_2.mat'); 
        case 38,load('5p120a_3.mat');disp('5p120a_3.mat'); 
        case 39,load('5p120a_4.mat');disp('5p120a_4.mat'); 
        case 40,load('5p120a_5.mat');disp('5p120a_5.mat');                     
    end

    if i < 16 || i==21 || i==22 || i==23 || i==24 || i==25
      CISA=2%10000;
      nrep=2%100;
    else
      CISA=2%20000;
      nrep=2%200;
    end
    
for j=1:alg % Para cada uno de los algoritmos
                  
for k=1:runs  % Para la acantidad de corridas
    
    if j==6
        vns=2; 
    else
        vns=1;
    end
                                               
     [Data(k,j,i),BestSchedule] = SA(DataProject,RRHH,TA,CP,CAType,CISA,Alpha,nrep,j,vns);
                                 
        if  Data(k,j,i) < BestFO(j,i)
                 
             BestFO(j,i)=Data(k,j,i);
             
             switch j 
                 case 1    
                   Best1(i)=BestSchedule;
                 case 2
                   Best2(i)=BestSchedule;
                 case 3
                   Best3(i)=BestSchedule;
                 case 4
                   Best4(i)=BestSchedule;
                 case 5
                   Best5(i)=BestSchedule;
                 case 6 
                   Best6(i)=BestSchedule;
             end

        end
            
           save('Data','Data','Best1','Best2','Best3','Best4','Best5','Best6','BestFO')

end

end % Para cada ualg de los algoritmos
end% Para cada ualg de las instancias
toc