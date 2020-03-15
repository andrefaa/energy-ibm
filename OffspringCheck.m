function [Y,O,I]=OffspringCheck(S,Y,O,I,Class)
VAL=[];
POSVAL=[];
    %Checagem inicial de filhotes que atingiram ao minimo o minimo possivel
    if any(O(:,4)>=min(S(:,13)))
        OVal=O(O(:,4)>=min(S(:,13)),:);
        POSO=find(O(O(:,4)>=min(S(:,13)),:));
        %Checagem se os filhotes atingiram o NeonateBS da sua especie
        for i=1:size(OVal,1)
            if OVal(i,5)>=S(OVal(i,3),11)
                VAL=vertcat(VAL,i);
                POSVAL=vertcat(POSVAL,POSO(i));
            end
        end
        %Inclusao dos que atingiram na matriz de juvenis
        if ~isempty(VAL)
            Gr=[];
            for i=1:size(VAL,1)
                if isempty(I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),:))
                    continue
                end
                
                if ~isempty(Y) && any(Y(:,2) == OVal(VAL(i),3))
                    Gr(i,1) = max(Y(Y(:,2) == OVal(VAL(i),3),1))+i;
                else 
                    Gr(i,1) = 1;
                end
                Gr(i,2) = OVal(VAL(i),2);
                Gr(i,3) = OVal(VAL(i),3);
                Gr(i,4) = 0;
                Gr(i,5) = S(OVal(VAL(i),3),11);
                if strcmp(Class,'Mammals')
                    Gr(i,6) = (3.98*(Gr(i,5)^0.686)*24)*0.0201;
%                     Gr(i,6)= (4.82*(Gr(i,5)^0.734));
                elseif strcmp(Class,'Birds')
                    Gr(i,6)= (10.5*(Gr(i,5)^0.681));
                elseif strcmp(Class,'Reptiles')
                    Gr(i,6)= (0.196*(Gr(i,5)^0.889));
                end
                Gr(i,7) = 39.3*((75*(Gr(i,5)/1000)^1.19)+(75*(Gr(i,5)/1000)^1.19*(-0.5+S(Gr(i,3),6))));
                Gr(i,8) = Gr(i,7);
                Gr(i,9) = 1;
            end
            %Junta Filhotes na matriz de juvenis e remove da matriz de filhotes
            Y=vertcat(Y,Gr);
            O(POSVAL,:)=[];
        end
    end