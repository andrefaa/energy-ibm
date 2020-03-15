function [I,O]=JuvenileCheck(S,I,O,Class)
VAL=[];
    %Checagem inicial de filhotes que atingiram ao minimo o minimo possivel
    if any(O(:,4)>=min(S(:,11)))
        OVal=O(O(:,4)>=min(S(:,11)),:);
        %Checagem se os filhotes atingiram o NeonateBS da sua especie
        for i=1:size(OVal,1)
            if OVal(i,5)>=S(OVal(i,3),11)
                VAL=vertcat(VAL,i);
            end
        end
        %Inclusao dos que atingiram na matriz de individuos
        if ~isempty(VAL)
            Gr=[];
            for i=1:size(VAL,1)
                if isempty(I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),:))
                    continue
                end
                Gr(i,1) = max(I(I(:,2) == OVal(VAL(i),3),1))+i;
                Gr(i,2) = OVal(VAL(i),3);
                Gr(i,3) = OVal(VAL(i),4);
                Gr(i,4) = OVal(VAL(i),5);
                if strcmp(Class,'Mammals')
                    Gr(i,5)= (4.82*(Gr(i,4)^0.734));
                elseif strcmp(Class,'Birds')
                    Gr(i,5)= (10.5*(Gr(i,4)^0.681));
                elseif strcmp(Class,'Reptiles')
                    Gr(i,5)= (0.196*(Gr(i,4)^0.889));
                end
                Gr(i,6) = (39.3*((75*(Gr(i,4)/1000)^1.19)*S(Gr(i,2),6)))/2;
                Gr(i,7) = 39.3*(75*(Gr(i,4)/1000)^1.19)*S(Gr(i,2),6);
                Gr(i,8) = I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),8);
                Gr(i,9) = I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),9);
                Gr(i,10) = I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),10);
                Gr(i,11) = I(find(I(:,2)==OVal(VAL(i),3) & I(:,1)==OVal(VAL(i),2)),11);
                Gr(i,12) = 0;
                Gr(i,13) = S(OVal(VAL(i),3),13);
                Gr(i,14) = 1;
                Gr(i,15) = 0;
            end
            %Junta Filhotes na matriz de individuos e remove da matriz de filhotes
            I=vertcat(I,Gr);
            O(VAL,:)=[];
        end
    end