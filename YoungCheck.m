function [A,Y]=YoungCheck(S,A,Y,Class)
VAL=[];
POSVAL=[];
    %Checagem inicial de filhotes que atingiram o tamanho de desmame minimo possivel
    if any(Y(:,5)>=min(S(:,16)))
        YVal=Y(Y(:,5)>=min(S(:,16)),:);
        POSY=find(Y(Y(:,5)>=min(S(:,16)),:));
        %Checagem se os filhotes atingiram o WeaningMass da sua especie
        for i=1:size(YVal,1)
            if YVal(i,5)>=S(YVal(i,3),16)
                VAL=vertcat(VAL,i);
                POSVAL=vertcat(POSVAL,POSY(i));
            end
        end
        %Inclusao dos que atingiram na matriz de individuos
        if ~isempty(VAL)
            Gr=[];
            for i=1:size(VAL,1)
                if isempty(A(find(A(:,2)==YVal(VAL(i),3) & A(:,1)==YVal(VAL(i),2)),:))
                    continue
                end
                Gr(i,1) = max(A(A(:,2) == YVal(VAL(i),3),1))+i;
                Gr(i,2) = YVal(VAL(i),3);
                Gr(i,3) = YVal(VAL(i),4);
                Gr(i,4) = YVal(VAL(i),5);
                if strcmp(Class,'Mammals')
                    Gr(i,5) = (3.98*(Gr(i,4)^0.686)*24)*0.0201; %White & Seymour, 2005
%                     Gr(i,5)= (4.82*(Gr(i,4)^0.734));
                elseif strcmp(Class,'Birds')
                    Gr(i,5)= (10.5*(Gr(i,4)^0.681));
                elseif strcmp(Class,'Reptiles')
                    Gr(i,5)= (0.196*(Gr(i,4)^0.889));
                end
                Gr(i,6) = YVal(VAL(i),7);
                Gr(i,7) = YVal(VAL(i),8);
                Gr(i,8) = A(find(A(:,2)==YVal(VAL(i),3) & A(:,1)==YVal(VAL(i),2)),8);
                Gr(i,9) = A(find(A(:,2)==YVal(VAL(i),3) & A(:,1)==YVal(VAL(i),2)),9);
                Gr(i,10) = A(find(A(:,2)==YVal(VAL(i),3) & A(:,1)==YVal(VAL(i),2)),10);
                Gr(i,11) = A(find(A(:,2)==YVal(VAL(i),3) & A(:,1)==YVal(VAL(i),2)),11);
                Gr(i,12) = 0;
                Gr(i,13) = S(Gr(i,2),15)-Gr(i,3);
                Gr(i,14) = 1;
                Gr(i,15) = 0;
                Gr(i,16) = 0;
                Gr(i,17) = ceil((2.3*(Gr(i,4)/1000)^1.02)*0.01);
            end
            %Junta Filhotes na matriz de individuos e remove da matriz de filhotes
            if ~isempty(Gr)
                A=vertcat(A,Gr);
                %Correcao de erro louco que nao quero corrigir agora (esta
                %criando uma linha só com zeros)
                A=A(A(:,14)==1,:);
            end
            Y(POSVAL,:)=[];
        end
    end