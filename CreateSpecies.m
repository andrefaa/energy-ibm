function S=CreateSpecies(BSMax,BSMin,Mundo,Grid,i)

S=zeros(1,18);

%Species ID
S(:,1)=i;

%Species Maximum Body Size
if strcmp(Mundo,'Neutro')
    S(:,2)=(BSMax-BSMin)+ BSMin;
else
    S(:,2)=round((BSMax-BSMin).*rand + BSMin); 
end
    
%Especialização por Recurso
% if strcmp(Mundo,'Neutro')
    S(:,3)=0.5;
% else
%     S(:,3)=rand;
% end
S(:,4)=1-(S(:,3));
    
%Estrategia de Alocacao de Recurso (5=Crescimento,6=Sobrevivencia,(1-6+7)=Reproducao)
aloc=rand(2,1);
S(:,5)=min(aloc);
S(:,6)=abs(min(aloc)-max(aloc));

%Species Maximum Body Size w/ Energy Priority (Chages Gompertz Assymptote)
S(:,2)=S(:,2)+(S(:,2)*(-0.5+S(:,5)));
    
%Origem Geografica da Especie (Grid)
S(:,7)=ceil(rand*Grid);
S(:,8)=ceil(rand*Grid);
    
%Origem Geografica da Especie (Coord)
S(:,9)=round(rand*1000);
S(:,10)=round(rand*1000);  
    
%Neonate Body Size (Martin & MacLarnon 1985)
% S(:,11) = 0.89*(S(:,2).^0.8);

%Neonate Body Size (Ernest et al 2003;Blueweiss et al 1978)
S(:,11) = 0.097*(S(:,2).^0.92);
    
%Maximum Litter Size (Ernest et al 2003; Blueweiss et al. 1978) --> Afetado por BS+Erep
% S(:,12)= round((0.55*(S(:,2).^0.82)./S(:,11))+((0.55*(S(:,2).^0.82)./S(:,11)).*(0.5-(S(:,5)+S(:,6)))));
S(:,12)= round(5.997*(S(:,2).^-0.142)+(5.997*(S(:,2).^-0.142).*(-0.5+(1-(S(:,5)+S(:,6))))));

%Gestation Time (Ernest 2003;Blueweiss et al 1978) --> Afetado por BS+Erep
S(:,13) = round(11.659*(S(:,2).^0.249)+(11.659*(S(:,2).^0.249).*(0.5-(1-(S(:,5)+S(:,6))))));

%Period Between Reproductive Season (De Marco Tese)
% S(:,14) = round(0.577*((S(:,2)/1000)^0.171)*365);
S(:,14) = round((0.577*((S(:,2)/1000)^0.171)+(0.577*((S(:,2)/1000)^0.171).*(0.5-(1-(S(:,5)+S(:,6))))))*365);

    
%Sexual Maturity Age(AFR;Ernest 2003)
% S(:,15) = round(0.788*((S(:,2)/1000)^0.266)*((1-S(:,5)+S(:,6)))*365);
S(:,15) = round(44.470*(S(:,2).^0.276)+(44.470*(S(:,2).^0.276.*(0.5-(1-(S(:,5)+S(:,6)))))));

%Weaning Mass(Ernest et al 2003)
S(:,16) = 0.485*(S(:,2).^0.931);
    
%Longevity (Zullinger et al. 1984; De Marco Tese)
%De Marco, 1999
% S(:,16) = round((7.38*(S(:,2)/1000)^0.197+(k1*7.38*(S(:,2)/1000).^0.197).*(-0.5+S(:,6))))*365);
%Zullinger et al 1984
S(:,17) = round(795.478*(S(:,2).^0.219)+(795.478*(S(:,2).^0.219).*(-0.5+S(:,6))));
     
%Species K and I (Gompertz Growth Equation) (Zullinger etl. al 1984;De Marco Tese)
%De Marco, 1999
% S(:,17) =0.014*((S(:,2)/1000)^-0.373);
% S(:,int18) = 66.547*((S(:,2)/1000)^0.219);
%Zullinger & Rickelfs 1984
S(:,18) =0.125*(S(:,2)^-0.302)+(0.125*(S(:,2)^-0.302).*(-0.5+S(:,5)));%Crescimento afeta velocidade de crescimento
S(:,19) = 5.343*(S(:,2)^0.354)+(0.125*(S(:,2)^-0.302).*(-0.5+S(:,6)));
% S(:,19) = 5.343*(S(:,2)^0.354)+(5.343*(S(:,2).^0.354).*(-0.5+(1-(S(:,5)+S(:,6)))));

%Fator de correcao de Gompertz
S(:,20) = S(:,2).*exp(-exp(-S(:,18).*(0-S(:,19))));