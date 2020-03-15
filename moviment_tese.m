function [R,I] = moviment_tese(R,I,nograph,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax)
% R - é um mapa da distribuição da quantidade de recursos
% I - é o indivíduo que ira forragear
% nograph - plota(0) ou não os scatterplots
% DMD - Total Daily Movement Distance of an individual
% IGMax - Maximum energy that can be ingested in a day
% ICL - Incremental cost of locomotion (energy lost in the process of
% foraging)
% Rmax - Food intake(g) per time(min)
% Vmax - Individual foraging velocity (m/min)
% hmax - Time taken in each bite (min/bite)
% Smax - Maximum bite size (g)

if nograph==0
    figure(1)
    imagesc(R);
end
[lin,cols]=size(R);

%Definir pontos de origem e destino
Porigin = I([10,11]);
Ptarget=Porigin;

%o ponto de saida é P origin
mov=1;
vid=VideoWriter('C:\Users\decoa\OneDrive\Doutorado\Apresentacao\Teste_Video\Caminhos12.avi');
vid.FrameRate=10;vid.Quality=100;
open(vid)
if nograph==0
    hold on
    scatter(Porigin(1),Porigin(2),90,[1 0 0],'filled')
    scatter(Ptarget(1),Ptarget(2),90,[1 0 0],'s','filled')
    IM=getframe;
	writeVideo(vid,IM);
end

%Primeiro movimento
P1=[];
DMDt=DMD;
IGMaxt=IGMax;
z=1;
R1=R;%Matriz de recursos espelhada
[movx,movy]=geramovrecurso(Porigin(1,:),R1);
P1(1,1)=Porigin(1,1)+movx;
P1(1,2)=Porigin(1,2)+movy;
%Correcao da posicao de P
if any(P1(1,:)<1) | (P1(1,1)>cols) | (P1(1,2)>lin)
	if P1(z,1)<1 & P1(z,2)<1
        P1(z,1) = cols;
        P1(z,2) = lin;
    elseif P1(z,1)>cols & P1(z,2)>lin
        P1(z,1) = 1;
        P1(z,2) = 1;
	elseif P1(z,1)<1 & P1(z,2)>lin
        P1(z,1) = cols;
        P1(z,2) = 1;
	elseif P1(z,1)>cols & P1(z,2)<1
        P1(z,1) = 1;
        P1(z,2) = lin;
    elseif P1(z,1)<1
        P1(1,1) = cols;
	elseif P1(1,2)<1
        P1(1,2) = lin;
	elseif P1(1,1)>cols
        P1(1,1)=1;
	elseif P1(z,2)>lin
        P1(1,2)=1;
	end
end
%Processo de forrageio entre as duas celuals (distancia percorrida
%pelo individuo(m) e energia consumida
[R,I,DMDt,IGMaxt] = forage_cells(R1,P1,I,z,DMDt,IGMaxt,ICL,Rmax,Vmax,hmax,Smax);

if nograph==0 
	scatter([Porigin(1),P1(z,1)],[Porigin(2),P1(z,2)],[],[1 0 0],'filled')
    IM=getframe;
	writeVideo(vid,IM);
	pause(1)
end

G(z,1)=(IGMax-IGMaxt);%quantidade de recurso ingerida
G(z,2)=(DMD-DMDt);%distancia percorrida em cada tempo
G(z,3)=G(1,1);%Recursos acumulados
G(z,4)=G(1,2);%Movimento acumulado

while (DMDt~=0) | P1(z,:)~=Ptarget %para quando chega no Ptarget ou quando terminar o tempo
%     if DMDt ==0
%             strcat('Tá dando ruim...')
%     end
    z=z+1;
    if any(Porigin~=Ptarget) %dispersal moviment
    %    [movx,movy]=geramov(mov,0);
    %    [movx,movy]=geramovvolta(mov,P1(z-1,:),Ptarget,z,lin,cols);
        [movx,movy]=geramovquant(P1(z-1,:),Ptarget,z,d0+1,R);
        P1(z,1)=P1(z-1,1)+movx;
        P1(z,2)=P1(z-1,2)+movy;
    else %central place foraging
      if DMDt>1*DMD/2
        %na primeira metade do caminho ele está andando de acordo com a
        %quantidade de recursos na vizinhança
        [movx,movy]=geramovrecurso(P1(z-1,:),R1);
        P1(z,1)=P1(z-1,1)+movx;
        P1(z,2)=P1(z-1,2)+movy;
        %Correcao da posicao de P
        if any(P1(z,:)<1) | (P1(z,1)>cols) | (P1(z,2)>lin)
          if P1(z,1)<1 & P1(z,2)<1
              P1(z,1) = cols;
              P1(z,2) = lin;
          elseif P1(z,1)>cols & P1(z,2)>lin
              P1(z,1) = 1;
              P1(z,2) = 1;
          elseif P1(z,1)<1 & P1(z,2)>lin
              P1(z,1) = cols;
              P1(z,2) = 1;
          elseif P1(z,1)>cols & P1(z,2)<1
              P1(z,1) = 1;
              P1(z,2) = lin;
          elseif P1(z,1)<1
              P1(z,1) = cols;
          elseif P1(z,2)<1
              P1(z,2) = lin;
          elseif P1(z,1)>cols
              P1(z,1)=1;
          elseif P1(z,2)>lin
              P1(z,2)=1;
          end
        end
        j=0;
        while (any((P1(1:z-1,1)==P1(z,1)) & (P1(1:z-1,2)==P1(z,2)))) || (any((P1(z,1)==Ptarget(1,1)) & (P1(z,2)==Ptarget(1,2))))
            j=j+1;
            [movx,movy]=geramovrecurso(P1(z-1,:),R);
            P1(z,1)=P1(z-1,1)+movx;
            P1(z,2)=P1(z-1,2)+movy;
            %Correcao da posicao de P
            if any(P1(z,:)<1) | (P1(z,1)>cols) | (P1(z,2)>lin)
              if P1(z,1)<1 & P1(z,2)<1
                P1(z,1) = cols;
                P1(z,2) = lin;
              elseif P1(z,1)>cols & P1(z,2)>lin
                P1(z,1) = 1;
                P1(z,2) = 1;
              elseif P1(z,1)<1 & P1(z,2)>lin
                P1(z,1) = cols;
                P1(z,2) = 1;
              elseif P1(z,1)>cols & P1(z,2)<1
                P1(z,1) = 1;
                P1(z,2) = lin;
              elseif P1(z,1)<1
                  P1(z,1) = cols;
              elseif P1(z,2)<1
                  P1(z,2) = lin;
              elseif P1(z,1)>cols
                  P1(z,1)=1;
              elseif P1(z,2)>lin
                  P1(z,2)=1;
              end
            end
            if j==100
                z=2;
            end
        end
%       Processo de forrageio entre as duas celuals (distancia percorrida
%       pelo individuo(m) e energia consumida
        IG1=IGMaxt;
        DMD1=DMDt;
        [R,I,DMDt,IGMaxt] = forage_cells(R1,P1,I,z,DMDt,IGMaxt,ICL,Rmax,Vmax,hmax,Smax);
      else
        %na segunda metade do caminho ele retorna mais direcionado para casa
        [movx,movy]=geramovquant(P1(z-1,:),Ptarget,DMDt,DMD,R,lin,cols);
        P1(z,1)=P1(z-1,1)+movx;
        P1(z,2)=P1(z-1,2)+movy;
        %Correcao da posicao de P
        if any(P1(z,:)<1) | (P1(z,1)>cols) | (P1(z,2)>lin)
            if P1(z,1)<1 & P1(z,2)<1
              P1(z,1) = cols;
              P1(z,2) = lin;
          elseif P1(z,1)>cols & P1(z,2)>lin
              P1(z,1) = 1;
              P1(z,2) = 1;
          elseif P1(z,1)<1 & P1(z,2)>lin
              P1(z,1) = cols;
              P1(z,2) = 1;
          elseif P1(z,1)>cols & P1(z,2)<1
              P1(z,1) = 1;
              P1(z,2) = lin;
          elseif P1(z,1)<1
                P1(z,1) = cols;
            elseif P1(z,2)<1
                P1(z,2) = lin;
            elseif P1(z,1)>cols
                P1(z,1)=1;
            elseif P1(z,2)>lin
                P1(z,2)=1;
            end
        end
%       Processo de forrageio entre as duas celuals (distancia percorrida
%       pelo individuo(m) e energia consumida
        IG1=IGMaxt;
        DMD1=DMDt;
        [R,I,DMDt,IGMaxt] = forage_cells(R1,P1,I,z,DMDt,IGMaxt,ICL,Rmax,Vmax,hmax,Smax);
      end
    end
   
%   G(z,1)=forage;%quantidade de recurso ingerido nas células visitadas
%   G(z,2)=dist;%distancia percorrida em cada tempo
    G(z,1)=(IG1-IGMaxt);%quantidade de recurso na célula visitada
    G(z,2)=(DMD1-DMDt);%distancia percorrida em cada tempo
    G(z,3)=G(z,1)+G(z-1,3);%Recursos acumulados --> Fazer um comparativo com IGMax
    G(z,4)=G(z,2)+G(z-1,4);%distância acumulada --> Fazer um comparativo com DMD

    if DMDt<0
        DMDt = 0;
    end
    if IGMaxt<0
        IGMaxt = 0;
    end
         
    if nograph==0
        if (DMDt>(DMD/2))
            scatter(P1(z,1),P1(z,2),[],[1 0 0],'filled')
            IM=getframe;
            writeVideo(vid,IM);
        else
            scatter(P1(z,1),P1(z,2),[],[0 0 0],'filled')
            IM=getframe;
            writeVideo(vid,IM);
        end
        pause(2)
    end
end
close(vid)
if nograph==0
    plot(P1(:,1),P1(:,2),'r-')
    hold off
    figure(2)
    plot(G(:,4),G(:,3),'-ro')
    xlabel('Distância percorrida');ylabel('Recursos acumulados')
end
% [tmax,c]=size(G(:,1));
% dmax=G(tmax,3);
% rmax=G(tmax,4);
   
