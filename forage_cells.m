function [R,I,DMD,IGMax] = forage_cells(R,P1,I,t,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax)
%Function to caculate dispersal movement among cells (in m) and energy
%consumed in the process

%Gravar locais de origem
if t==1
    CellX = I(10);
    CellY = I(11);
else
    CellX = P1(t-1,1);
    CellY = P1(t-1,2);
end
%Nova celula
NewCellX = P1(t,1);
NewCellY = P1(t,2);
%Posicao dentro da celula
PosX = I(8);
PosY = I(9);
%Sortear Posição na Nova Célula
NewX = round(1000*rand());
NewY = round(1000*rand());

%Identificar a celula para qual o individuo se locomoveu e a distancia que
%percorreu na celula atual até alcançar a borda
DCell = 0;
DNew = 0;
if NewCellX>CellX & NewCellY>CellY
    DCell = sqrt((0-PosY)^2+(1000-PosX)^2);
    DNew = sqrt((1000-NewY)^2+(0-NewX)^2);
elseif NewCellX<CellX & NewCellY<CellY
    DCell = sqrt((1000-PosY)^2+(0-PosX)^2);
    DNew = sqrt((0-NewY)^2+(1000-NewX)^2);
elseif NewCellX>CellX & NewCellY<CellY
    DCell = sqrt((1000-PosY)^2+(1000-PosX)^2);
    DNew = sqrt((0-NewY)^2+(0-NewX)^2);
elseif NewCellX<CellX & NewCellY>CellY
    DCell = sqrt((0-PosY)^2+(0-PosX)^2);
    DNew = sqrt((1000-NewY)^2+(1000-NewX)^2);    
elseif NewCellX>CellX
    DCell = 1000-PosX;
    DNew = NewX;
elseif NewCellX<CellX
    DCell = PosX;
    DNew = 1000-NewX;
elseif NewCellY>CellY
    DCell = PosY;
    DNew = 1000-NewY;
elseif NewCellY<CellY
    DCell = 1000-PosY;
    DNew = NewY;
end

%Energy consumed in the foraging step
% DensCell = 1/sqrt(((R(CellY,CellX)/17)/100000000000)*10^6);
% DensNew = 1/sqrt((R(NewCellY,NewCellX)/17)/100000);
%Tempo em cada celula
TCell = DCell/Vmax;
TNewCell = DNew/Vmax;

%Energia na celula (Consumo+Gasto com movimento)
% if DensCell>((Vmax*Smax)/Rmax)
    IG_Cell = TCell*(((Vmax*sqrt((R(CellY,CellX)/17)/10000)*Smax)/(1+hmax*Vmax*sqrt((R(CellY,CellX)/17)/10000)))*10);
% else
% 	IG_Cell = TCell*((Rmax*Smax)/(Rmax*hmax+Smax))*10;
% end
Gasto_Cell=ICL*DCell;
LiquiCell = IG_Cell-Gasto_Cell;

%Energia Consumida na nova celula
% if DensNew>((Vmax*Smax)/Rmax)
    IG_New = TNewCell*(((Vmax*sqrt((R(NewCellY,NewCellX)/17)/10000)*Smax)/(1+hmax*Vmax*sqrt((R(NewCellY,NewCellX)/17)/10000)))*10);
% else
% 	IG_New = TNewCell*((Rmax*Smax)/(Rmax*hmax+Smax))*10;
% end
Gasto_New=ICL*DNew;
LiquiNew = IG_New-Gasto_New;

%Corrige excessos de consumo
if IG_Cell>R(CellY,CellX)
    IG_Cell = R(CellY,CellX);
end
if IG_New>R(NewCellY,NewCellX)
    IG_New = R(NewCellY,NewCellX);
end
%Remove a energia do sistema
R(CellY,CellX) = R(CellY,CellX)-IG_Cell;
R(NewCellY,NewCellX) = R(NewCellY,NewCellX)-IG_New;

%Atualiza energia consumida, posicao, budget de movimento e de ingestão
I(8) = NewX;
I(9) = NewY;
I(10) = NewCellX;
I(11) = NewCellY;
I(16) = I(16)+LiquiCell+LiquiNew;
DMD = DMD - (DCell+DNew);
IGMax = IGMax - (IG_Cell+IG_New);
