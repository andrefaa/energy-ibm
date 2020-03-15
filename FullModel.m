function FullModel(NSP,N,E1,Grid,Tn,tempo,Esp,Class,Mundo,BSMin,BSMax,Propagule,NProp,IRep,Figures,DirO,K)

%Parametros iniciais:
%NSP  - numero de especies no inicio do modelo
%N - numero de individuos/especie
%E1 - energia inicial do recurso 1 em todo o sistema
%E2 - energia inicial do recurso 2 em todo o sistema
%Grid - Numero de celulas na regiao (Grid x Grid)
%Tn - tempo para estabilizar a comunidade nativa
%tempo - tempo após estabilizacao da comunidade nativa (3650 dias de
%estabilizacao)
%Esp - Tipo de especialistas(Chegada/Assimilacao)
%Class - Classe de organismos no sistema (Mammals/Birds/Reptiles)
%Mundo - Neutro / Nicho
%BSMin - tamamho corporal mínimo na comunidade (kg)
%BSMax - tamamho corporal máximo na comunidade (kg)
%Propagule - type of propagule pressure (Single/Multiple)
%NProp - number of propagules introduced (Single=once/Multiple=daily for a
%year)
%IRep - number of replicates for different invaders meeting the same native
%community
%Figures - display figures (0/1)
%DirO - Diretorio de outputs
%K - replicate(used for metamodelling)


%CONDICOES INICIAIS
K=num2str(K);
%Number format
format short g
%Species Matrix (S)
S=zeros(NSP,20);

%Conversion to Kelvin Temperature
% T=T+273.15;

%Result matrix
RES=zeros(tempo,NSP+1);

%S Matrix Columns
%[Sp_ID  Species Maximum Body Size(g)  FMR(KJ/Day)(Nagy et al 1999)  RaioE1  RaioE2  
% Estrategia(Cresc) Estrategia(Sobrev) POSL POSC POSX POSY NeonatoBodySize TempoGestacao(dias) 
% TempoentrePeriodoReprodutivo Longevidade K(Gompertz) I(Gompertz) FatorCorrecao(Gompertz)]

%Filling Species Matrix(S)
strcat('Creating Species')
for i=1:NSP
    S(i,:)=CreateSpecies(BSMax,BSMin,Mundo,Grid,i);
end

%Creating the World
strcat('Creating World')
G1=rand(Grid,Grid);
sg1=sum(sum(G1));
G1=G1/sg1;
G1=E1*G1;
G0=G1;

% G2=randn(Grid,Grid);
% sg2=sum(sum(G2));
% G2=G2/sg2;
% G2=E2*G2;

% G=G1+G2;
% figure(1)
% image(G1)
% % figure(2)
% % image(G2)
% % figure(3)
% % image(G)
% pause

%Adult Matrix (A)

%Columns Explanation
%[IND_ID  Sp_ID  Age_IND Body Size_IND IND_FMR  Actual Reserve(KJ/g of Fat)  Maximum Reserve(KJ/g of Fat)
% POSX POSY PosROW PosCOL]
strcat('Creating Individuals')
A=zeros((NSP*N),17);
k=0;
for i=1:NSP
    for j=1:N
        k=k+1;
        A(k,:)=CreateIndividual(j,i,S,Class);
    end
end

%Corrigir Posicao
A=CheckPOS(A,Grid);

% %Plot do sistema
% j=0;
% for u=1:Grid
%     for p=1:Grid
%         j=j+1;
%         figure(j)
%         temp=A(A(:,10)==p & A(:,11)==u,:);
%         scatter(temp(:,8), temp(:,9),5,temp(:,2))
%     end
% end
% pause

%Matriz de Offspring, Young, Juvenile
O=[];
Y=[];

%Matriz SIZE & STRAT
SIZE=[];
STRAT=[];

% % DINAMICA % %
%tempo
strcat('Starting Model')
for t=1:Tn
    t
    STRAT(t,1)=mean(S(A(:,2),5));
    STRAT(t,2)=mean(S(A(:,2),6));
    STRAT(t,3)=mean(1-(S(A(:,2),5)+S(A(:,2),6)));
    SIZE(t,1)=size(A,1);
    SIZE(t,2)=size(Y,1);
    SIZE(t,3)=size(O,1);
    
    if Figures==1
        figure(2)

        hist(A(:,2),NSP)
        title(num2str(t))

        pause(0.0001)
    end
    
    %Resource Consumption
%     strcat('Foraging...')
    [A,G1]=ResourceConsumption(S,A,Y,G1,Esp,Mundo);
    
    %Resource Allocation
%     strcat('Allocating resources...')
    [A,O,Y,G1]=ResourceAllocation_New(S,A,O,Y,G1);

    %Atualizar posicoes
%     strcat('Updating positions...')
    A=CheckPOS(A,Grid);

    %Checa se houve dispersao
%     strcat('Individuals dispersal...')
    if any(A(:,15)==1)
        IDisp=A(A(:,15)==1,:);
        A=A(A(:,15)==0,:);
        [IDisp,O,Y,G1]=ResourceAllocation_New(S,IDisp,O,Y,G1);
        A=vertcat(A,IDisp);
    end

    %Envelhecimento
%     strcat('Individuals Ageing...')
    A(:,3)=A(:,3)+1;
    if ~isempty(O)
        O(:,4)=O(:,4)+1;
    end
    if ~isempty(Y)
        Y(:,4)=Y(:,4)+1;
    end
    
    %Longevidade/Morte
%     strcat('Killing old individuals...')
    OLD = find(A(:,3) >= S(A(:,2),17));
%     if ~isempty (O) & ~isempty(OLD)
%         O_OLD=find(O(:,3)==A(OLD,2) & O(:,2)==A(OLD,1));
%         if ~isempty(O_OLD)
%             O(O_OLD,7)=0;
%         end
%     end
%     if ~isempty (Y) & ~isempty(OLD)
%         Y_OLD=find(Y(:,3)==A(OLD,2) & Y(:,2)==A(OLD,1)); 
%         if ~isempty(J_OLD)
%             J(J_OLD,9)=0;
%         end
%     end
    
    if ~isempty(OLD)
        A(OLD,14)=0;
    end


    %Checar se filhotes viraram juvenis
    if ~isempty(O)
%         strcat('Checking offspring birth...')
        [Y,O,A]=OffspringCheck(S,Y,O,A,Class);
    end
    
    %Checa se juvenis viraram adultos
    if ~isempty(Y)
%         strcat('Juvenile weaning...')
        [A,Y]=YoungCheck(S,A,Y,Class);
    end

    %Zerar a dispersao e a alimentacao
    A(:,12)=0;
    A(:,15)=0;

    %Reduzir o reproduction timer
%     strcat('Decreasing reproduction timer...')
    A(:,13)=A(:,13)-1;
    if any(A(:,13)<0)
        
        RESET=find(A(:,13)<0);
        A(RESET,13)=S(A(RESET,2),14);
    end
    
    %Update Gestating Females Reserves (125%)
    if ~isempty(O)
%         strcat('Increasing female reserves...')
        for m=1:size(O,1)
            POS_M=find(A(:,2)==O(m,3) & A(:,1)==O(m,2));
            A(POS_M,7) = (39.3*((75*(A(POS_M,4)/1000).^1.19)+((75*(A(POS_M,4)/1000).^1.19).*(-0.5+S(A(POS_M,2),6)))))*1.25;
        end
    end
    
    %Replenish resources
%     strcat('Replenishing resources...')
    G1=G0;

    %Save species size
    for i=1:NSP+1
        RES(t,i)=size(A(A(:,2)==i,:),1);
    end
    
    if Figures==1
        figure(3)
        plot(SIZE(:,1))
        hold on
        plot(SIZE(:,2))
        hold on
        plot(SIZE(:,3))
        hold off
        grid
        legend('Adults', 'Young', 'Offspring')

        %Plot Strategies
        figure(4)
        plot(STRAT(:,1))
        hold on
        plot(STRAT(:,2))
        hold on
        plot(STRAT(:,3))
        hold off
        grid
        legend('Growth', 'Survival', 'Reproduction')
    end
end %Fecha o ciclo diario pre-invasao

xlswrite(strcat(DirO,'NativeCommunity_',K,'.xlsx'),A)
xlswrite(strcat(DirO,'CommunityStrategy_',K,'.xlsx'),STRAT)

% %Plot dda comunidade nativa
% j=0;
% for u=1:Grid
%     for p=1:Grid
%         j=j+1;
%         figure(j)
%         temp=A(A(:,10)==p & A(:,11)==u,:);
%         scatter(temp(:,8), temp(:,9),5,temp(:,2))
%     end
% end
% figure(5)
% plot(RES(:,:))

%Invasao
Result=zeros(IRep,12);
strcat('Invasion Replicates')
for m=1:IRep
    m/IRep
    RES1=RES;
    I1=A;
	S1=S;
	O1=O;
    Y1=Y;
	S1(NSP+1,:)=CreateSpecies(BSMax,BSMin,Mundo,Grid,NSP+1);        
    for t=(Tn+1):tempo
        t
        %Inserir individuos da especie invasora
        if strcmp(Propagule,'Single')
            TypeProp=1;
            if t==(Tn+1)
                k=size(I1,1);
                for a=1:NProp
                    k=k+1;
                    I1(k,:)=CreateIndividual(a,NSP+1,S1,Class);
                end
            end
        end
            
        if strcmp(Propagule,'Multiple')
            TypeProp=2;
            if t>=(Tn+1) && t<=(Tn+365)
                k=size(I1,1);
                for a=1:NProp
                    k=k+1;
                    I1(k,:)=CreateIndividual(a,NSP+1,S1,Class);
                end 
            end
        end
        
        %Corrigir Posicao
        I1=CheckPOS(I1,Grid);
        
        %Resource Consumption
        [I1,G1]=ResourceConsumption(S1,I1,Y1,G1,Esp,Mundo);
    
        %Resource Allocation
        [I1,O1,Y1,G1]=ResourceAllocation_New(S1,I1,O1,Y1,G1);

        %Atualizar posicoes
        I1=CheckPOS(I1,Grid);

        %Checa se houve dispersao
%         strcat('Individuals dispersal...')
        if any(I1(:,15)==1)
            IDisp=I1(I1(:,15)==1,:);
            I1=I1(I1(:,15)==0,:);
            [IDisp,O1,Y1,G1]=ResourceAllocation_New(S1,IDisp,O1,Y1,G1);
            I1=vertcat(I1,IDisp);
        end

        %Envelhecimento
%         strcat('Individuals Ageing...')
        I1(:,3)=I1(:,3)+1;
        if ~isempty(O1)
            O1(:,4)=O1(:,4)+1;
        end
        if ~isempty(Y)
            Y1(:,4)=Y1(:,4)+1;
        end

        %Longevidade/Morte
%         strcat('Killing old individuals...')
        OLD = find(I1(:,3) >= S1(I1(:,2),17));
    %     if ~isempty (O) & ~isempty(OLD)
    %         O_OLD=find(O(:,3)==I1(OLD,2) & O(:,2)==I1(OLD,1));
    %         if ~isempty(O_OLD)
    %             O(O_OLD,7)=0;
    %         end
    %     end
    %     if ~isempty (Y) & ~isempty(OLD)
    %         Y_OLD=find(Y(:,3)==I1(OLD,2) & Y(:,2)==I1(OLD,1)); 
    %         if ~isempty(J_OLD)
    %             J(J_OLD,9)=0;
    %         end
    %     end

        if ~isempty(OLD)
            I1(OLD,14)=0;
        end


        %Checar se filhotes viraram juvenis
        if ~isempty(O1)
%             strcat('Checking offspring birth...')
            [Y1,O1,I1]=OffspringCheck(S1,Y1,O1,I1,Class);
        end

        %Checa se juvenis viraram adultos
        if ~isempty(Y1)
%             strcat('Juvenile weaning...')
            [I1,Y1]=YoungCheck(S1,I1,Y1,Class);
        end

        %Zerar a dispersao e a alimentacao
        I1(:,12)=0;
        I1(:,15)=0;

        %Reduzir o reproduction timer
%         strcat('Decreasing reproduction timer...')
        I1(:,13)=I1(:,13)-1;
        if any(I1(:,13)<0)
            RESET=find(I1(:,13)<0);
            I1(RESET,13)=S1(I1(RESET,2),14);
        end

        %Update Gestating Females Reserves (125%)
        if ~isempty(O1)
%             strcat('Increasing female reserves...')
            for n=1:size(O1,1)
                POS_M=find(I1(:,2)==O1(n,3) & I1(:,1)==O1(n,2));
                I1(POS_M,7) = (39.3*((75*(I1(POS_M,4)/1000).^1.19)+((75*(I1(POS_M,4)/1000).^1.19).*(-0.5+S1(I1(POS_M,2),6)))))*1.5;
            end
        end

        %Replenish resources
%         strcat('Replenishing resources...')
        G1=G0;

        %Save species size
        for i=1:NSP+1
            RES1(t,i)=size(I1(I1(:,2)==i,:),1);
        end
        %Check if invader is alive
        if RES1(t,NSP+1)==0
            break
        end
    end %Fecha o ciclo diario da invasao

    NCom=I1(I1(:,2)==NSP+1,:);
    if isempty(NCom)
        NCom=0;
    else
        NCom=sum(sum(crosstab(NCom(:,10),NCom(:,11))>0));
    end
    %Matriz de Resultados
    %Abundancia , Numero Comunidades Ocupadas ,Body Size , Especialization , Egrwth, Esurv, Erepr, OriginX,
    %OriginY, Tipo de Propagulo (Single/Multiple), Numero de Propagulos,
    %Tempo que  especie invasora ficou no sistema
    Result(m,:)=[RES1(t,NSP+1), NCom ,S1(NSP+1,2), S1(NSP+1,3), S1(NSP+1,5), S1(NSP+1,6),(1-(S1(NSP+1,5)+S1(NSP+1,6))), S1(NSP+1,7), S1(NSP+1,8), TypeProp, NProp, t];
    xlswrite(strcat(DirO,'AbundanciaComunidade_',K,'_Invasor_',num2str(m),'.xlsx'),RES1)
    disp(Result(1:m,:));
end %Fecha as replicas da invasao
xlswrite(strcat(DirO,'Result_',K,'.xlsx'),Result)
% xlswrite(strcat(DirO,'Abundancia_',K,'.xlsx'),RES1)
xlswrite(strcat(DirO,'SpeciesCommunity_',K,'.xlsx'),S)
