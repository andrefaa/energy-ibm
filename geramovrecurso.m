function [movx,movy]=geramovrecurso(P,R)
%gera o movimento como uma escolha que maximiza apenas o uso de rescursos
%com uma escolha aleatória entre as quatro opções com maior recurso
   [lin,cols]=size(R);d=zeros(8,1);r=zeros(8,1);
   
   m=[-1 -1;-1 0; -1 1;1 -1;1 0;1 1;0 -1;0 1];
   for i=1:8
      O(i,1)=P(1)+m(i,1);
      O(i,2)=P(2)+m(i,2);
      %Corrigir a posicao no grid
      if any(O(i,:)<1) || (O(i,1)>cols) || (O(i,2)>lin)
          if O(i,1)<1 & O(i,2)<1
              O(i,1) = cols;
              O(i,2) = lin;
          elseif O(i,1)>cols & O(i,2)>lin
              O(i,1) = 1;
              O(i,2) = 1;
          elseif O(i,1)<1 & O(i,2)>lin
              O(i,1) = cols;
              O(i,2) = 1;
          elseif O(i,1)>cols & O(i,2)<1
              O(i,1) = 1;
              O(i,2) = lin;
          elseif O(i,1)<1
              O(i,1) = cols;
          elseif O(i,2)<1
              O(i,2) = lin;
          elseif O(i,1)>cols
              O(i,1)=1;
          elseif O(i,2)>lin
              O(i,2)=1;
          end
      end
      r(i)=R(O(i,2),O(i,1));
   end
   [r,I]=sort(r,'descend');
   plim=randi(4);
   
   movx=m(I(plim),1);
   movy=m(I(plim),2);