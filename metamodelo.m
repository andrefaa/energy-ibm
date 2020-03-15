parpool(4)
parfor K=1:20
    FullModel(300,10,2000000,10,3650,7300,'Chegada','Mammals','Nicho',10,100000,'Single',50,50,0,...
        'C:\Users\decoa\OneDrive\Doutorado\Bloco2-IBM\Cap3-Competicao\Resultados_2020_01_26\',K)
end
poolobj = gcp('nocreate');
delete(poolobj);