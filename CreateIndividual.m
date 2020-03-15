function It=CreateIndividual(j,i,S,Class)

It=zeros(1,17);
%Individual ID and SpeciesID
It(:,1)=j;
It(:,2)=i;
        
%Individual Age
It(:,3) = round(S(It(:,2),17)*0.7);
        
%Individual Body Size according to Gompertz Equation and
%Allocation Strategy (Cresc->K ; Sobrev->I)
It(:,4) = S(It(:,2),2)*exp(-exp(-S(It(:,2),18)*(It(:,3)-S(It(:,2),19))))-(S(It(:,2),20)-S(It(:,2),11));
        
%Individual Field Metabolic Rate (kJ/dia)
%Define Organism Class (Nagy et al 1999; Nagy et al 2005)
if strcmp(Class,'Mammals')
    It(:,5) = (3.98*(It(:,4)^0.686)*24)*0.0201; %White & Seymour, 2005
%   It(:,5)= (7.94*(It(:,4)^0.646));%Herbivores(Nagy 1999)
% 	It(:,5)= (4.82*(It(:,4)^0.734));
elseif strcmp(Class,'Birds')
	It(:,5)= (10.5*(It(:,4)^0.681));
elseif strcmp(Class,'Reptiles')
	It(:,5)= (0.196*(It(:,4)^0.889));
end
        
%Actual and Maximum Reserve(kj/g) (Starts w/ 50% reserve) 
%(Lindstedt & Schaeffer 2002)
It(:,6)= (39.3*((75*(It(:,4)/1000)^1.19)+((75*(It(:,4)/1000)^1.19)*(-0.5+S(It(:,2),6)))))/2;
It(:,7)= 39.3*((75*(It(:,4)/1000)^1.19)+((75*(It(:,4)/1000)^1.19)*(-0.5+S(It(:,2),6))));
        
%Individual Geographic Position (Coord) --> Individual dispersion from
%centre of origin (Santini et al 2013; Mean Distance -> Body Size Herbivores+Omnivores)
Vec=[-1 1];
MaxDist=(1.07*(It(:,4)/1000)^0.68)*1000;
% ActualDist = 0+ (MaxDist- 0)*sum(rand(2,1),2)/1;
ActualDist = [MaxDist;MaxDist];
It(:,8)= (datasample(Vec,1)*ActualDist(1,1))+(S(It(:,2),9));
It(:,9)= (datasample(Vec,1)*ActualDist(2,1))+(S(It(:,2),9));
        
%Individual Geographic Position (Grid)
It(:,10) = S(It(:,2),7);
It(:,11) = S(It(:,2),8);
        
%FED
It(:,12) = 0;
        
%Reproduction Timer
It(:,13) = S(It(:,2),13);
if It(:,13)<0
    It(:,13)=0;
end
        
%Alive
It(:,14)=1;
        
%Dispersed
It(:,15)=0;

%Energy Ingested
It(:,16)=0;

%Home Range Size(number of cells;Hudson & White 1985)
It(:,17)=ceil((2.3*(It(:,4)/1000)^1.02)*0.01);


