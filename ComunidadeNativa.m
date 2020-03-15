function [S,G1,G2,I] = ComunidadeNativa(N,E1,Grid,tempo,Esp,Class,BSMin,BSMax,imgsim,Mundo)
%FUNCAO  PARA ENTENDER A ORIGEM, O UNIVERSO E TUDO MAIS
%function GCHH(b1,b2,d1,d2,tol,S,raio,p_disturbio)
%
%NSP  - numero de especies no inicio do modelo
%E1 - energia inicial do recurso 1 em todo o sistema
%E2 - energia inicial do recurso 2 em todo o sistema
%Grid - Numero de celulas na regiao (Grid x Grid)
%tempo - tempo para estabilizacao da comunidade nativa
%T - Temperatura do sistema
%Class - Classe de organismos no sistema (Mammals/Birds/Reptiles)
%BSMin - tamamho corporal mínimo na comunidade
%BSMax - tamamho corporal máximo na comunidade
%imgsim - faz ou não a imagem
%Mundo - Neutro / Nicho

%CONDICOES INICIAIS

%Number format
format short g
%Species Matrix (S)
S=zeros(N,18);

%Conversion to Kelvin Temperature
% T=T+273.15;

%S Matrix Columns
%[Sp_ID  Species Maximum Body Size(g)  FMR(KJ/Day)(Nagy et al 1999)  RaioE1  RaioE2  
% Estrategia(Cresc) Estrategia(Sobrev) POSL POSC POSX POSY NeonatoBodySize TempoGestacao(dias) 
% TempoentrePeriodoReprodutivo Longevidade K(Gompertz) I(Gompertz) FatorCorrecao(Gompertz)]

%Filling Species Matrix(S)
for i=1:N
    %Species ID
    S(i,1)=i;
    %Species Maximum Body Size
    if strcmp(Mundo,'Neutro')
        S(i,2)=(BSMax-BSMin)+ BSMin;
    else
       S(i,2)=round((BSMax-BSMin).*rand + BSMin); 
    end
    
    %Especialização por Recurso
    if strcmp(Mundo,'Neutro')
        S(i,3)=0.5;
    else
        S(i,3)=rand;
    end
    S(i,4)=1-(S(i,3));
    
    %Estrategia de Alocacao de Recurso (5=Crescimento,6=Sobrevivencia,(1-6-7)=Reproducao)
    sort=rand(2,1);
    S(i,5)=min(sort);
    S(i,6)=abs(min(sort)-max(sort));
    
    %Origem Geografica da Especie (Grid)
    S(i,7)=round(rand*Grid);
    S(i,8)=round(rand*Grid);
    
    %Origem Geografica da Especie (Coord)
    S(i,9)=round(rand*1000);
    S(i,10)=round(rand*1000);  
    
    %Maximum Neonate Body Size (Martin & MacLarnon 1985)
    S(i,11) = 0.89*(S(i,2)^0.8)*S(i,5);
    
    %Litter Size
    S(i,12)= ceil((0.55*(S(i,2)/1000)^0.82)*(1-(S(i,5)+S(i,6))));

    %Gestation Time (De Marco Tese)
    %S(i,12) = round(0.15*((S(i,2)/1000)^0.275)*365);

    %Period Between Reproductive Season (De Marco Tese)
    S(i,13) = round(0.577*((S(i,2)/1000)^0.171)*365);
    
    %Sexual Maturity
    S(i,14) = round(0.788*((S(i,2)/1000)^0.266)*(S(i,5)+S(i,6))*365);
    
    %Longevity (De Marco Tese)
    S(i,15) = round(7.38*((S(i,2)/1000)^0.197)*S(i,6)*365);
        
    %Species K and I (Gompertz Growth Equation) (De Marco Tese)
    S(i,16) =0.014*((S(i,2)/1000)^-0.373);
    S(i,17) = 66.547*((S(i,2)/1000)^0.219);
    
    %Fator de correcao de Gompertz
    S(i,18) = (S(i,2)*exp(-exp((-S(i,16)*S(i,5))*((0-(S(i,17)*S(i,5)))))));
    
end

%Creating the World
G1=randn(Grid,Grid);
sg1=sum(sum(G1));
G1=G1/sg1;
G1=E1*G1;

% G2=randn(Grid,Grid);
% sg2=sum(sum(G2));
% G2=G2/sg2;
% G2=E2*G2;

% G=G1+G2;
% figure(1)
% image(G1)
% figure(2)
% image(G2)
% figure(3)
% image(G)
% pause

%Individual Matrix (I)

%Columns Explanation
%[IND_ID  Sp_ID  Age_IND Body Size_IND  Actual Reserve(KJ/g of Fat)  Maximum Reserve(KJ/g of Fat)
% POSX POSY PosROW PosCOL]

I=zeros((N*20),14);
k=0;
for i=1:N
    for j=1:20
        k=k+1;
        %Individual ID and SpeciesID
        I(k,1)=j;
        I(k,2)=i;
        
        %Individual Age
        I(k,3) = round(S(I(k,2),15)*0.7);
        
        %Individual Body Size according to Gompertz Equation and
        %Allocation Strategy (Cresc->K ; Sobrev->I)
        I(k,4)= (S(I(k,2),2)*exp(-exp((-S(I(k,2),16)*S(I(k,2),5))*((I(k,3)-(S(I(k,2),17)*S(I(k,2),5)))))))-S(I(k,2),18);
        
        %Individual Field Metabolic Rate
            %Define Organism Class (Nagy et al 1999)
        if strcmp(Class,'Mammals')
            I(k,5)= (4.82*(I(k,4)^0.734));
        elseif strcmp(Class,'Birds')
            I(k,5)= (10.5*(I(k,4)^0.681));
        elseif strcmp(Class,'Reptiles')
            I(k,5)= (0.196*(I(k,4)^0.889));
        end
        
        %Actual and Maximum Reserve(kj/g) (Starts w/ 50% reserve) 
        %(Lindstedt & Schaeffer 2002)
        I(k,6)= (39.3*((75*(I(k,4)/1000)^1.19)*S(I(k,2),6)))/2;
        I(k,7)= 39.3*(75*(I(k,4)/1000)^1.19)*S(I(k,2),6);
        
        %Individual Geographic Position (Coord)
        I(k,8)= (I(k,4)/100).*randn(1)+(S(I(k,2),9));
        I(k,9)= (I(k,4)/100).*randn(1)+(S(I(k,2),10));
        
        %Individual Geographic Position (Grid)
        I(k,10) = S(I(k,2),7);
        I(k,11) = S(I(k,2),8);
        
        %FED
        I(k,12) = 0;
        
        %Reproduction Timer
        I(k,13) = S(I(k,2),13);
        
        %Alive
        I(k,14)=1;
        
        %Dispersed
        I(k,15)=0;
    end
end

%Corrigir Posicao
I=CheckPOS(I,Grid);

%Matriz de Filhotes
O=[];
 
% % DINAMICA % %
%tempo
for t=1:tempo 
  f=size(O,1)+1;
% Resource Allocation
    for l=1:Grid
        for c=1:Grid
            I1=I(I(:,9)==l & I(:,10)==c,:);
            if isempty(I1)
                continue
            end
            if strcmp(Esp,'Chegada')
                %Ordem de especialistas (Especialistas possuem a chance de alcançar
                %a comida antes)
                ISp=tabulate(I1(:,2));
                ISp=ISp(ISp(:,2)~=0,1);
                ISp=horzcat(ISp,S(ISp,4));
                if strcmp(Mundo,'Nicho')
                    ISp=sortrows(ISp,2,'descend');
                else
                    ISp=ISp(randperm(size(ISp,1)),:);
                end
                for s=1:size(ISp,1)
                   I2=I1(I1(:,2)==ISp(s,1),:);
                   %Intraespécie: Ordem dos individuos comerem é na sorte
                   I2=I2(randperm(size(I2,2)),:);
                   for Ind=1:size(I2,2)
                    %Ingestion (Clauss et al 2007 - Compar. Biochem. &
                    %Phisiol.) - Comida=Carboidrato/Proteina(17.6; Sibly et al.
                    %2013)
                    IG=(((0.047*((I2(Ind,4)/1000)^0.76))*1000)*17.6);
                    
                    %Check if there is energy at the cell
                    if IG<G1(l,c)
                        IG=0;
                    end

                    %Daily Metabolic Requirements
                    IG=IG-I2(Ind,5);
                    if (IG>0)
                        I2(Ind,12)=1;
                        G1(l,c)=G1(l,c)-IG;
                    else
                        I2(Ind,8)= (I2(Ind,4)/100).*randn(1)+(S(I2(Ind,2),9));
                        I2(Ind,9)= (I2(Ind,4)/100).*randn(1)+(S(I2(Ind,2),10));
                        I2(Ind,15)=1;
                    end

                    % % Allocation Priority % %
                    P=[S(I2(Ind,2),5:6) 1-sum(S(I2(Ind,2),5:6))];
                    [out,idx]=sort(P,'descend');

                    %%%%%% DEFINIR ORDEM DE ALOCACAO DE ACORDO COM IDX %%%%%%

                    %%% Growth %%% (GMax According to Gompertz)
                    NextW = (S(I2(Ind,2),2)*exp(-exp((-S(I2(Ind,2),16)*S(I2(Ind,2),5))*(((I2(Ind,3)+1)-(S(I2(Ind,2),17)*S(I2(Ind,2),5)))))))-S(I2(Ind,2),18);
                    GMax=NextW-I2(Ind,4);
                    %Energy required for creating tissues
                    %1g of flesh = 7kj + 6kj for synthesizing (Siby et al.
                    %2013, Peters 1983, Moses et al 2008)
                    EGMax=GMax*13;
                    if((IG-EGMax)>0)
                        IG=IG-EGMax;
                    else
                        GMax=IG/13;
                        IG=0;
                    end

                    %Consolidate Growth
                    I2(Ind,4)=I2(Ind,4)+GMax;
                    
                    %%% Fecha Crescimento %%%

                    %%% Reproduction %%% (Determinate Growers)
                    %Ver se o individuo possui filhotes & alimentar
                    %filhotes
                       if (~isempty(O) && any(O(:,3)==I2(Ind,2) & O(:,2)==I2(Ind,1)))
                           %Checa o numero de filhos
                           POS=find(O(:,3)==I2(Ind,2) & O(:,2)==I2(Ind,1));
                           %Crescimento total dos filhotes 
                           NOFW=(S(I2(Ind,2),2)*exp(-exp((-S(I2(Ind,2),16)*S(I2(Ind,2),5))*(((O(POS,4)+1)-(S(I2(Ind,2),17)*S(I2(Ind,2),5)))))))-S(I2(Ind,2),18);
                           NGMax=NOFW-O(POS,5);
                           %Energy required for creating tissues
                           %1g of flesh = 7kj + 6kj for synthesizing (Siby et al.
                           %2013, Peters 1983, Moses et al 2008)
                           NEGMax=NGMax*13;
                           %Alimentou o suficiente para dar energia aos
                           %filhotes
                           if IG>=sum(NEGMax)
                               O(POS,5)=NOFW;
                               O(POS,4)=O(POS,4)+1;
                               IG=IG-sum(NEGMax);
                           %Nao alimentou o suficiente para dar energia
                           %aos filhotes
                           else
                               %Alimenta proporcional aos que conseguirem
                               POST= find(IG>cumsum(NEGMax));
                               O(POST,5)=NOFW;
                               O(POST,4)=O(POST,4)+1;
                               IG=IG-sum(NEGMax(POST,:));
                               %Retira o restante da reserva, caso haja
                               %reserva
                               if I2(Ind,6)>=sum(NEGMax(IG<cumsum(NEGMax)))
                                    O(POS(POS~=POST),5)=NOFW;
                                    O(POS(POS~=POST),4)=O(POS(POS~=POST),4)+1;
                                    I2(Ind,6)=I2(Ind,6)-sum(NEGMax(IG<cumsum(NEGMax)));
                               %Caso nao haja reserva, absorve os filhotes 
                               %e soma sua energia à reserva     
                               else
                                    O(POS(POS~=POST),6)=0;
                                    I2(Ind,6)=I2(Ind,6)+(O(POS(POS~=POST),5)*13);
                               end
                           end   
                       end
                    
                       %Producao de Filhotes (dias>maturidade e repr. timer=0)
                       if I2(Ind,3)>=S(I2(Ind,2),14) && I2(Ind,13)==0
                         %Alocacao de energia por filhote
                         W0 = (S(I2(Ind,2),2)*exp(-exp((-S(I2(Ind,2),16)*S(I2(Ind,2),5))*((1-(S(I2(Ind,2),17)*S(I2(Ind,2),5)))))))-S(I2(Ind,2),18);
                         EW=W0*14;
                         %Numero de Filhotes
                         LS=S(I2(Ind,2),12);
                         %Maximum Energy of Litter
                         LE=EW*LS;
                         %Offspring production
                         %Se houver energia para todos os filhotes
                         if(IG-LE)>0
                            IG=IG-LE;
                            %Adicionar filhos na matriz
                            for o=1:LS
                                OF=[o I2(Ind,1) I2(Ind,2) 1 W0];
                                O(f,:)=OF;
                                f=f+1;
                            end
                         %Se nao houver energia para todos os filhotes
                         else
                            %Produz filhotes proporcional ao que alimentou
                            LST=round(IG/EW);
                            LE=EW*LST;
                            IG=IG-LE;
                            for o=1:LST
                                OF=[o I2(Ind,1) I2(Ind,2) 1 W0 1];
                                O(f,:)=OF;
                                f=f+1;
                            end
                            %Retira energia das reservas para produzir
                            %filhotes
                            if I2(Ind,6)>=((LS-LST)*EW)
                                LE=((LS-LST)*EW);
                                I2(Ind,6)=I2(Ind,6)-LE;
                                for o=1:(LS-LST)
                                    OF=[o I2(Ind,1) I2(Ind,2) 1 W0 1];
                                    O(f,:)=OF;
                                    f=f+1;
                                end
                            end
                         end
                         I2(Ind,13)=S(I2(Ind,2),13); %Reseta reproduction timer
                       end
                   %%% Fecha Reproducao %%%

                   %%% Survival %%%(Reserves)
                    if (IG<=(I2(Ind,7)-I2(Ind,6)))
                        I2(Ind,6)=I2(Ind,6)+IG;
                        IG=0;
                    else
                        I2(Ind,6)=I2(Ind,7);
                        IG=IG-(I2(Ind,7)-I2(Ind,6));
                    end 
                    %%% Fecha Survival/Reservas %%%
                    
                    %Checa se toda a energia ingerida foi utilizada
                    %Se não, retorna a energia remanescente ao sistema
                    if(IG~=0)
                        G1(l,c)=G1(l,c)+IG;
                    end
                    
                    %Update Marix
                    I(I(:,1)==I2(Ind,1) & I(:,2)==I2(Ind,2),:)=I2(Ind,:);
                   end %Fecha for intraespecifico
                end %Fecha for interespecifico
            end %Fecha especializacao por chegada
            if strcmp(Esp,'Assimilacao')
                %Ordem de especialistas (A chance de alcançar é a mesma,
                %especialistas possuem maior eficiencia energetica)
            end %Fecha especializacao por assimilacao
        end %Fecha a coluna
    end %Fecha a linha
    
    %Checa se houve dispersao
    if any(I(:,15)==1)
        IDisp=I(I(:,15)==1,:)
        ResourceAllocation(IDisp);
    end
end %Fecha o ciclo diario
end %Fecha a funcao


%Consome reservas
%I2(Ind,6)= I2(Ind,6)-(I2(Ind,5)-IG);