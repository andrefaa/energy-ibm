function [A] = CheckPOS(A,Grid)
%Argumentos
%A = Matriz de Individuos
%Grid = Tamanho da Regiao (proveniente da funcao superior)

%Check Species Position
%Move to the Right
pos=find(A(:,8)>1000);
while size(pos,1) ~=0
    A(pos,8)=A(pos,8)-1000;
    A(pos,10)=A(pos,10)+1;
    A(A(:,10)>Grid,10)=1;
    pos=find(A(:,8)>1000);
end

%Move to the Left
pos=find(A(:,8)<0);
while size(pos,1) ~=0
    A(pos,8)=A(pos,8)+1000;
    A(pos,10)=A(pos,10)-1;
    A(A(:,10)<1,10)=Grid;
    pos=find(A(:,8)<0);
end

%Move Up
pos=find(A(:,9)>1000);
while size(pos,1) ~=0
    A(pos,9)=A(pos,9)-1000;
    A(pos,11)=A(pos,11)-1;
    A(A(:,11)<1,11)=Grid;
    pos=find(A(:,9)>1000);
end
%Move Down
pos=find(A(:,9)<0);
while size(pos,1) ~=0
    A(pos,9)=A(pos,9)+1000;
    A(pos,11)=A(pos,11)+1;
    A(A(:,11)>Grid,11)=1;
    pos=find(A(:,9)<0);
end