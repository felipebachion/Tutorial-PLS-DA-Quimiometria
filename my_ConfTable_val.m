function [Tval Tval2] = my_ConfTable_val(yval,yprev_val,ts)
% Rotina desenvolvida para calcular a tabela de confusão/contingencia do conjunto de
% amostras da validação
% INPUT : 
%          >yval       Referência das classes das amostras de validação
%          >yprev_val  Valores previstos pelo modelo das amostras de
%          validação
%          >ts         Limiares entre as classes
% OUTPUT:
%          > Tabela de confusão 1 e 2

YVAL = vet_matrix(yval);
yprev_val_new=[];
nc=yprev_val-repmat(ts,[size(yprev_val,1),1]);
for i = 1:size(yprev_val,1)
[ww,yprev_val_new(i,1)]=max(nc(i,:));
end

class_param_val = calc_class_param(yprev_val_new,YVAL)
C1v=class_param_val.conf_mat(:,1);
C2v=class_param_val.conf_mat(:,2);
C3v=class_param_val.conf_mat(:,3);
C4v=class_param_val.conf_mat(:,4);
Sem_Classev=class_param_val.conf_mat(:,5);
figmerit={'Classe_Orignal_Oliva','Classe_Orignal_Canola','Classe_Orignal_Milho','Classe_Orignal_Soja'};
Tval=table(C1v,C2v,C3v,C4v,Sem_Classev,'RowNames',figmerit)
Tval.Properties.VariableNames = {'Previsto_Oliva','Previsto_Canola','Previsto_Milho','Previsto_Soja', 'Previsto_Sem_Classe'}
writetable(Tval,'Resultados_PLSDA.xlsx','Sheet',1,'Range','I1','WriteRowNames',true)
clear C1v C2v C3v C4v

C1v=[class_param_val.accuracy;class_param_val.sensitivity(:,1);class_param_val.specificity(:,1)];
C2v=[class_param_val.accuracy;class_param_val.sensitivity(:,2);class_param_val.specificity(:,2)];
C3v=[class_param_val.accuracy;class_param_val.sensitivity(:,3);class_param_val.specificity(:,3)];
C4v=[class_param_val.accuracy;class_param_val.sensitivity(:,4);class_param_val.specificity(:,4)];
figmerit={'Acurácia','Sensibilidade','Especificidade'};
Tval2=table(C1v,C2v,C3v,C4v,'RowNames',figmerit)
Tval2.Properties.VariableNames = {'Azeite_de_Oliva', 'Oleo_de_Canola', 'Oleo_de_Milho', 'Oleo_de_Soja'}
writetable(Tval2,'Resultados_PLSDA.xlsx','Sheet',1,'Range','I8','WriteRowNames',true)

end

