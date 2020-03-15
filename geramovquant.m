function [movx,movy]=geramovquant(P,Ptarget,DMDt,DMD,R,lin,cols) 
%gera o movimento como uma escolha que minimiza a distância para o Ptarget
%cada vez mais de acordo com o numero de passos
   dmax=lin*cols;
   [lin,cols]=size(R);d=zeros(8,1);r=zeros(8,1);
   d0=sqrt(((P-Ptarget)*(P-Ptarget)'));
   m=[-1 -1;-1 0; -1 1;1 -1;1 0;1 1;0 -1;0 1];
   for i=1:8
      O(i,1)=P(1)+m(i,1);
      O(i,2)=P(2)+m(i,2);    
      if any(O(i,:)<1) | (O(i,1)>cols) | (O(i,2)>lin)
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
      dt(1,1)=(((O(i,1)-Ptarget(1,1))^2)+((O(i,2)-Ptarget(1,2))^2))/dmax;
      dt(1,2)=((lin-O(i,1))^2+((1-Ptarget(1,1))^2+(O(i,2)-Ptarget(1,2))^2))/dmax;
      dt(1,3)=((1-O(i,1))^2+((lin-Ptarget(1,1))^2+(O(i,2)-Ptarget(1,2))^2))/dmax;
      dt(1,4)=((cols-O(i,2))^2+((O(i,1)-Ptarget(1,1))^2+(1-Ptarget(1,2))^2))/dmax;
      dt(1,5)=((1-O(i,2))^2+((O(i,1)-Ptarget(1,1))^2+(cols-Ptarget(1,2))^2))/dmax;
      d(i)=dt(find(dt==min(dt),1));
   end
%     [d,A]=sort(d,'descend');
% 	plim=floor(1+8*((DMD-DMDt)/DMD));
% 	if plim>8
%         plim=8;
% 	end
% 	d=d(A)
   i=find(d==min(d),1);
   movx=m(i,1);
   movy=m(i,2);