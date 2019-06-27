function [Tcal Tcal2] = my_ConfTable_cal(ycal,yprev_cal,ts)
% Rotina desenvolvida para calcular a tabela de confusão/contingencia do conjunto de
% amostras da calibração
% INPUT : 
%          >ycal       Referência das classes das amostras de calibração
%          >yprev_cal  Valores previstos pelo modelo das amostras de calibração
%          >ts         Limiares entre as classes
% OUTPUT:
%          > Tabela de confusão 1 e 2
          

YCAL = vet_matrix(ycal);
nc=yprev_cal-repmat(ts,[size(yprev_cal,1),1]);
yprev_cal_new=[];
    for i = 1:size(yprev_cal,1)
    [ww,yprev_cal_new(i,1)]=max(nc(i,:));
    end
class_param_cal = calc_class_param(yprev_cal_new,YCAL)
C1=class_param_cal.conf_mat(:,1);
C2=class_param_cal.conf_mat(:,2);
C3=class_param_cal.conf_mat(:,3);
C4=class_param_cal.conf_mat(:,4);
Sem_Classe=class_param_cal.conf_mat(:,5);
figmerit={'Classe_Orignal_Oliva','Classe_Orignal_Canola','Classe_Orignal_Milho','Classe_Orignal_Soja'};
Tcal=table(C1,C2,C3,C4,Sem_Classe,'RowNames',figmerit)
Tcal.Properties.VariableNames = {'Previsto_Oliva','Previsto_Canola','Previsto_Milho','Previsto_Soja', 'Previsto_Sem_Classe'}
writetable(Tcal,'Resultados_PLSDA.xlsx','Sheet',1,'Range','A1','WriteRowNames',true)
clear C1 C2 C3 C4

C1=[class_param_cal.accuracy;class_param_cal.sensitivity(:,1);class_param_cal.specificity(:,1)];
C2=[class_param_cal.accuracy;class_param_cal.sensitivity(:,2);class_param_cal.specificity(:,2)];
C3=[class_param_cal.accuracy;class_param_cal.sensitivity(:,3);class_param_cal.specificity(:,3)];
C4=[class_param_cal.accuracy;class_param_cal.sensitivity(:,4);class_param_cal.specificity(:,4)];
figmerit={'Acurácia','Sensibilidade','Especificidade'};
Tcal2=table(C1,C2,C3,C4,'RowNames',figmerit)
Tcal2.Properties.VariableNames = {'Azeite_de_Oliva','Oleo_de_Canola','Oleo_de_Milho','Oleo_de_Soja'}
writetable(Tcal2,'Resultados_PLSDA.xlsx','Sheet',1,'Range','A8','WriteRowNames',true)

end

