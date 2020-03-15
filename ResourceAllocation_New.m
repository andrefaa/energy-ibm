function [I,O,Y,G1] = ResourceAllocation_New(S,I,O,Y,G1)

% Resource Allocation

%Inicia LOOP Individuos
    for Ind=1:size(I,1)                
        if I(Ind,16)<0
            I(Ind,:)
            pause
        end
        
        %Check if Individual ate more than the BMR and disperse individuals
        %Individuals can supplement the energy ingested with the reserves
        if I(Ind,15)==1 && (I(Ind,16)<I(Ind,5)) &&((I(Ind,16)+I(Ind,6))<I(Ind,5))
            I(Ind,14)=0;
%             disp('Morreu')
        elseif I(Ind,15)==0 && (I(Ind,16)<I(Ind,5)) &&((I(Ind,16)+I(Ind,6))<I(Ind,5))
            [VecX,VecY]=geramovrecurso([I(Ind,10),I(Ind,11)],G1);
            MaxDist=(1.07*(I(Ind,4)/1000)^0.68)*1000;
%             ActualDist = 0+ (MaxDist- 0)*sum(rand(2,1),2)/1;
            ActualDist = [MaxDist;MaxDist];
            I(Ind,8)= (VecY*ActualDist(1,1))+(S(I(Ind,2),9));
            I(Ind,9)= (VecX*ActualDist(2,1))+(S(I(Ind,2),9));
            I(Ind,15)=1;
            disp('Dispersou!')
            continue
        elseif I(Ind,15)==1 && (I(Ind,16)<I(Ind,5)) &&((I(Ind,16)+I(Ind,6))>I(Ind,5))
            I(Ind,6)=((I(Ind,16)+I(Ind,6))-I(Ind,5));
            I(Ind,16)=0;
            disp('Ja andou, nao comeu o suficiente e tirou das reservas')
            continue
        elseif I(Ind,15)==0 && (I(Ind,16)<I(Ind,5)) && ((I(Ind,16)+I(Ind,6))>I(Ind,5))
            I(Ind,6)=((I(Ind,16)+I(Ind,6))-I(Ind,5));
            I(Ind,16)=0;
            [VecX,VecY]=geramovrecurso([I(Ind,10),I(Ind,11)],G1);
            MaxDist=(1.07*(I(Ind,4)/1000)^0.68)*1000;
%             ActualDist = 0+ (MaxDist- 0)*sum(rand(2,1),2)/1;
            ActualDist = [MaxDist;MaxDist];
            I(Ind,8)= (VecY*ActualDist(1,1))+(S(I(Ind,2),9));
            I(Ind,9)= (VecX*ActualDist(2,1))+(S(I(Ind,2),9));
%             disp('Removeu energia das reservas para completar o necessario e dispersou tentando encontrar um lugar melhor!')
            continue
        end
        if I(Ind,16)==0
            I(Ind,:)
        end
        
        %Cenarios de Alocação de Energia (Individuo comeu mais do que sua
        %taxa metabólica basal)
        %Daily Metabolic Requirements
        I(Ind,16)=I(Ind,16)-I(Ind,5);

        % % Allocation Priority % %
        P=[S(I(Ind,2),5:6) 1-sum(S(I(Ind,2),5:6))];
        [out,idx]=sort(P,'descend');

        for ordem=1:3
            switch idx(ordem)
                case 1
                    %%% Growth %%% (GMax According to Gompertz)
                    NextW = S(I(Ind,2),2)*exp(-exp(-S(I(Ind,2),18)*((I(Ind,3)+1)-S(I(Ind,2),19))))-(S(I(Ind,2),20)-S(I(Ind,2),11));
                    GMax=NextW-I(Ind,4);
                    %Energy required for creating tissues
                    %1g of flesh = 7kj + 6kj for synthesizing (Siby et al.
                    %2013, Peters 1983, Moses et al 2008)
                    EGMax=GMax*13;
                    if (I(Ind,16)-EGMax)>0
                        I(Ind,16)=I(Ind,16)-EGMax;
                    elseif (I(Ind,16)-EGMax)<0 && I(Ind,16)~=0
                        GMax=I(Ind,16)/13;
                        I(Ind,16)=0;
                    elseif (I(Ind,16)-EGMax)<0 && I(Ind,16)==0
                        GMax=0;
                    end
                    %Consolidate Growth
                    I(Ind,4)=I(Ind,4)+GMax;
                    %Update Maximum Reserves & BMR
                    I(Ind,5) = (3.98*(I(Ind,4)^0.686)*24)*0.0201;
                    I(Ind,7) = 39.3*((75*(I(Ind,4)/1000)^1.19)+(75*(I(Ind,4)/1000)^1.19*(-0.5+S(I(Ind,2),6))));
                %%% Fecha Crescimento %%%
                
                case 2
                    %%% Survival %%%(Reserves)
                    if (I(Ind,16)<=(I(Ind,7)-I(Ind,6)))
                        I(Ind,6)=I(Ind,6)+I(Ind,16);
                        I(Ind,16)=0;
                    else
                        I(Ind,6)=I(Ind,7);
                        I(Ind,16)=I(Ind,16)-(I(Ind,7)-I(Ind,6));
                    end
                    %%% Fecha Survival/Reservas %%%
                    
                case 3 
                    %%% Reproduction %%% (Determinate Growers)
                                 
                    %Ver se o individuo possui embriões(O) & alimentar
                    %embriões
                    if (~isempty(O) && any(O(:,3)==I(Ind,2) & O(:,2)==I(Ind,1)))
                        %Checa o numero de embriões
                        POS=find(O(:,3)==I(Ind,2) & O(:,2)==I(Ind,1));
                        %Crescimento total dos embriões
                        %Linear 
%                         NOFW=(S(I(Ind,2),2)*exp(-exp((-S(I(Ind,2),16)*S(I(Ind,2),5))*(((O(POS,4)+1)-(S(I(Ind,2),17)*S(I(Ind,2),5)))))))-S(I(Ind,2),18);
%                         NGMax=NOFW-O(POS,5);
                        NOFW = O(POS,5)+O(POS,6);
                        %Energy required for creating tissues
                        %1g of flesh = 7kj + 6kj for synthesizing (Siby et al.
                        %2013, Peters 1983, Moses et al 2008)
%                       NEGMax=NGMax*13;
                        NEGMax = O(POS,6)*13;
                        
                        %Alimentou o suficiente para dar energia aos
                        %filhotes
                        if I(Ind,16)>=sum(NEGMax)
                            O(POS,5)=NOFW;
                            I(Ind,16)=I(Ind,16)-sum(NEGMax);
                            %Nao alimentou o suficiente para dar energia
                            %aos filhotes
                        else
                            %Alimenta proporcional aos que conseguirem
                            POST=POS(I(Ind,16)>cumsum(NEGMax));
                            NOFWT=NOFW(I(Ind,16)>cumsum(NEGMax));
                            if ~isempty(POST)
                                O(POST,5)=NOFWT;
                                I(Ind,16)=I(Ind,16)-sum(NEGMax(I(Ind,16)>cumsum(NEGMax),:));
                            end
                            %Retira o restante da reserva, caso haja
                            %reserva
                            if I(Ind,6)>=sum(NEGMax(I(Ind,16)<cumsum(NEGMax)))
                                if ~isempty(POST)
                                    O(POST,5)=NOFWT;
                                else
                                    O(POS,5)=NOFW;
                                end
                                I(Ind,6)=I(Ind,6)-sum(NEGMax(I(Ind,16)<cumsum(NEGMax)));
                                %Caso nao haja reserva, absorve os filhotes 
                                %e soma sua energia à reserva     
                            else
                                if ~isempty(POST)
                                    O(POS,7)=0;
                                    I(Ind,6)=I(Ind,6)+sum((O(POS,5)*13));
                                else
                                    O(POST,7)=0;
                                    I(Ind,6)=I(Ind,6)+sum((O(POST,5)*13));
                                end
                            end
                        end
                    end
                    
                    %Checa se o individuo possui filhotes(Y) &
                    %alimenta os filhotes --> Amamentar
                    if (~isempty(Y) && any(Y(:,3)==I(Ind,2) & Y(:,2)==I(Ind,1)))
                        %Femea produz leite (g Leite/dia;Hamwell & Peaker
                        %1977) --> 
%                         MILK=(0.0835*(I(Ind,4)/1000)^0.765)*1000;
                        %Daily Milk Energy (KJ/dia;Hamwell & Peaker
                        %1977)
                        MILK_E = (127.2*(I(Ind,4)/1000)^0.694)*4.184;        
                        %Milk Energy/g
%                         MILK_G = MILK_E/MILK;
                        
                        %Checa se ha energia para produzir o leite e produz
                        %proporcional, retira das reservas para
                        %complementar
                        if (I(Ind,16) >= MILK_E)
                            I(Ind,16) = I(Ind,16) - MILK_E;
                        elseif (I(Ind,16) <= MILK_E) && ((I(Ind,16)+I(Ind,6))<=MILK_E)
                            MILK_E = (I(Ind,16)+I(Ind,6));
                            I(Ind,16) = 0;
                            I(Ind,6) = 0;
                        elseif (I(Ind,16) <= MILK_E) && ((I(Ind,16)+I(Ind,6))>=MILK_E)
                            I(Ind,6) = I(Ind,6)-(MILK_E-I(Ind,16));
                            I(Ind,16) = 0;
                        end
                                     
                        %Checa o numero de filhots juvenis
                        POS=find(Y(:,3)==I(Ind,2) & Y(:,2)==I(Ind,1));
                        [id,ord]=sort(Y(POS,5),'descend');
                        POS=POS(ord);
                        %Crescimento total dos juvenis
                        NOFW=S(Y(POS,3),2).*exp(-exp(-S(Y(POS,3),18).*((Y(POS,4)+1)-S(Y(POS,3),19))))-(S(Y(POS,3),20)-S(Y(POS,3),11));
                        NGMax=NOFW-Y(POS,5);
                        %Energy required for creating tissues
                        %1g of flesh = 7kj + 6kj for synthesizing (Siby et al.
                        %2013, Peters 1983, Moses et al 2008)
                        NEGMax=NGMax*13;
                        %Somar a isto o BMR dos filhotes
                        NEGMax = NEGMax + ((3.98*(Y(POS,5).^0.686)*24)*0.0201);
                                                              
                        %Alimentou o suficiente para dar
                        %energia para os filhotes suprirem
                        %seu BMR e crescerem
                        if MILK_E>=sum(NEGMax)
                            Y(POS,5)=NOFW;
                            %Update Maximum Reserves & BMR
                            Y(POS,6) = (3.98*(Y(POS,5).^0.686)*24)*0.0201;
                            Y(POS,8) = 39.3*((75*(Y(POS,5)/1000).^1.19)+(75*(Y(POS,5)/1000).^1.19.*(-0.5+S(Y(POS,3),6))));
                            
                            %Se sobrar leite preencher a reserva dos
                            %filhotes de forma proporcional(mesma
                            %quantidade para cada filhote)
                            MILK_E=MILK_E-sum(NEGMax);
                            if MILK_E>0
                                MILK_E=MILK_E/size(POS,1);
                                FAM = find(Y(POS,7)<Y(POS,8));
                                if ~isempty(FAM)
                                    Y(POS(FAM),7) = Y(POS(FAM),7)+MILK_E;
                                    if Y(POS(FAM),7)>Y(POS(FAM),8)
                                        Y(POS(FAM),7)=Y(POS(FAM),8);
                                    end
                                end
                            end
                            
                        %Nao alimentou o suficiente para dar energia
                        %aos filhotes
                        else
                            %O que fazer quando nao há energia na reserva para
                            %alimentar todos os filhotes?
                            %1-Filhotes maiores tem preferencia na ordem de alimentacao
                            %2-Filhotes podem retirar a energia de sua propria reserva
                            %3-Caso nao possuam reserva estes filhotes
                            %morrem
                            POST=POS(MILK_E>cumsum(NEGMax));
                            NOFWT=NOFW(MILK_E>cumsum(NEGMax));
                            if ~isempty(POST)
                                Y(POST,5)=NOFWT;
                                %Update Maximum Reserves & BMR
                                Y(POST,6) = (3.98*(Y(POST,5).^0.686)*24)*0.0201;
                                Y(POST,8) = 39.3*((75*(Y(POST,5)/1000).^1.19)+(75*(Y(POST,5)/1000).^1.19.*(-0.5+S(Y(POST,3),6))));
                                %Filhotes que não se alimentaram
                                POS(1:size(POST,1))=[];
                                NEGMax(1:size(POST,1))=[];
                                NOFW(1:size(NOFWT,1))=[];
                            end
                            %Checar se esses filhotes possuem reservas
                            %para suprir a falta de leite
                            if ~isempty(POS) & Y(POS,7)>=NEGMax%ERRO aqui!
                                Y(POS,7) = Y(POS,7) - NEGMax;
                                Y(POS,5)=NOFW;
                                %Update Maximum Reserves & BMR
                                Y(POS,6) = (3.98*(Y(POS,5).^0.686)*24)*0.0201;
                                Y(POS,8) = 39.3*((75*(Y(POS,5)/1000).^1.19)+(75*(Y(POS,5)/1000).^1.19.*(-0.5+S(Y(POS,3),6))));
                            %Quando nao possuem reservas estes filhotes
                            %morrem
                            elseif ~isempty(POS) & Y(POS,7) < NEGMax
                                Y(POS,9) = 0;
                            end
                        end
                    end
                    
                    % Producao de Filhotes (dias>maturidade e repr. timer=0) %
                    if (I(Ind,3)>=S(I(Ind,2),15)) && I(Ind,13)==0 && I(Ind,3)<=(0.9*S(I(Ind,2),17))
                        %Numero de Filhotes
                        LS=S(I(Ind,2),12);
                        %Adicionar filhos na matriz
                        f=size(O,1);%Contador filhotes
                        for o=1:LS
                            f=f+1;
                            OF=[o I(Ind,1) I(Ind,2) 0 0 (S(I(Ind,2),11)/S(I(Ind,2),13)) 1];
                            O(f,:)=OF;
                        end
                        I(Ind,13)=S(I(Ind,2),14); %Reseta reproduction timer
                    end
            end %%% Fecha Swith %%%
        end %%% Fecha Ordem de Alocação %%
        
        %Checa se toda a energia ingerida foi utilizada
        %Se não, retorna a energia remanescente ao sistema
        if I(Ind,16)~=0
            G1(I(Ind,11),I(Ind,10))=G1(I(Ind,11),I(Ind,10))+I(Ind,16);
            I(Ind,16)=0;
        end

        %Correcao matriz reservas
        if I(Ind,6)<0
            I(Ind,6)=0;
        end

        %Update Adult Marix
        I(I(:,1)==I(Ind,1) & I(:,2)==I(Ind,2),:)=I(Ind,:);
    end %%% Fecha Loop Individuos %%
    

       
%Matar individuos, embrioes e filhotes
I=I(I(:,14)==1,:);
if ~isempty(O)
    O=O(O(:,7)==1,:);
end
if ~isempty(Y)
    Y=Y(Y(:,9)==1,:);
end