function [I1,G1] = ResourceConsumption(S,A,Y,G1,Esp,Mundo)
%Funcao que controla a retirada de recurso do ambiente pelos individuos
%Adultos

I1=A;
%Ordem de consumo
if strcmp(Esp,'Chegada')
    %Ordem de especialistas (Especialistas possuem a chance de alcançar
	%a comida antes)
	ISp=I1(:,[1 2 4]);
	ISp=horzcat(ISp,S(ISp(:,2),3));
	ISp(:,3)=ISp(:,3)/max(ISp(:,3));
	%Ordem de prioridade é definida 70% pelo grau de especialização
	%(eficiencia para encontrar recursos) e 30% pelo tamanho
	%corporal(eficiencia para competir pelo recurso)
% 	ISp(:,5)=((0.75*ISp(:,4))+(0.25*ISp(:,3)))/2;-->Afetado por grau de
% 	especializacao em recurso + body size
    ISp(:,5)=ISp(:,4);
	if strcmp(Mundo,'Nicho')
%         [VALEspe,ORDEsp]=sort(ISp(:,5),'descend');
        ORDEsp=randperm(size(ISp,1));
        I1=I1(ORDEsp,:);
    else
        I1=I1(randperm(size(ISp,1)),:);
	end
           
    %inicia o LOOP de individus para consumo de recursos
	for Ind=1:size(I1,1)
       
        %IG_Max (KJ)--> Herbivores 10KJ/g (Nagy 2001; DMI g/day-->Veio do FMR)
        IGMax = (0.859*(I1(Ind,4)^0.628))*10;
        if ~isempty(Y) & ~isempty(Y(:,3)==I1(Ind,2) & Y(:,2)==I1(Ind,1))
            IGMax = IGMax*1.5;%Femeas com juvenis se alimentam 150% do normal
        end
        %DMD(Maximum Daily Movement Distance; m/dia; Garland 1983)
        DMD = 0.875*((I1(Ind,4)/1000)^0.22)*1000;
        if ~isempty(Y) & ~isempty(Y(:,3)==I1(Ind,2) & Y(:,2)==I1(Ind,1))
            DMD = DMD*1.5;%Femeas com juvenis se alimentam 150% do normal
        end
        %ICL(Incremental Cost of Locomotion; KJ/m Garland 1983
        ICL = ((10678*(I1(Ind,4)/1000)^0.7)/1000)/1000;
        %Define o Home Range do Individuo
%         HR_ext=I1(Ind,17)-1;
        %Parametros de Consumo de Recursos
        Rmax = 0.7483*(I1(Ind,4)/1000)^0.69;%(g/min ; Shipley et al 1994)
        Vmax = 52.16*(I1(Ind,4)/1000)^0.04; %(m/min ; Shipley 1996)
        hmax = 0.0118*(I1(Ind,4)/1000)^0.03;%(min/bite ; Shipley et al 1994)
        Smax = 0.0963*(I1(Ind,4)/1000)^0.71;%(g ; Shipley et al 1994)
        
        %Consome os recursos
        [G1,I1(Ind,:)] = moviment_tese(G1,I1(Ind,:),0,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax);
	end
    %Identifica que todos so individuos passaram pelo ciclo de alimentacao
    I1(:,12)=1;
elseif strcmp(Esp,'Eficiencia')
    %Nao implementado
end%Fecha tipo de especializacao por recurso
end%Fecha funcao de consumo de recursos