%% Rotina desenvolvida em função do Artigo publicado no periodico Quimica nova na escola
% Fazer o donwload de todas as funções necessárias em:
% 

%% 1) carregar o conjunto de dados
load amostras
%% 2) plotar o cojunto completo de espectros
figure(1)
subplot(211)
plot(num,X);
axis tight
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(212)
plot(num(:,1300:1750),X(:,1300:1750));
axis tight
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B','LineStyle','none');
[ax,h1]=suplabel('Número de onda (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('Absorbância','y'); 
set(h1,'FontSize',24,'FontName','Times New Roman') 
set(h2,'FontSize',24,'FontName','Times New Roman') 
%% 3) criar o vetor de classes para segregar os espectros de cada tipo de óleo
y=[ones(108,1);2*ones(54,1);3*ones(54,1);4*ones(54,1)];
%% 4) usar o algoritmo CALTESTDA para separar as amostras dos conjuntos de calibração, teste e preprocessar
% observe que uma região um pouco maior que a região da impressão digital é
% selecionada.
% [~,xcal,xval,ycal,yval]=caltestda(X(:,1300:1750),y,70,'k',[],{'deriv';[9,2,1]});
[~,xcal,xval,ycal,yval]=caltestda(X(:,1300:1750),y,70,'k',[],[]);
numc=num(1300:1750); % Realizando o corte da faixa espectral.
%% 5) plotar os espectros de calibração e validação após a seleção da região espectral de interesse
figure(2)
subplot(2,1,1);
plot(numc,xcal)
axis tight
set(gca,'xtick',[],'FontSize',18);...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A');
title('Amostras de Calibração','fontsize',14);
subplot(2,1,2);
plot(numc,xval)
axis tight
set(gca,'xtick',[],'FontSize',18);...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B');
title('Amostras de Validação','fontsize',14); hold on
[ax,h1]=suplabel('Número de onda (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('Absorbância','y',[.11 .11 .84 .84]); 
set(h1,'FontSize',24) 
set(h2,'FontSize',24) 
%% 6) Usar a função vet_matrix para transformar um vetor y1 em uma matriz y2
%   y1[1        y2[1 0 0
%      1           1 0 0
%      2           0 1 0
%      2           0 1 0
%      3           0 0 1
%      3]          0 0 1]

ycal = vet_matrix(ycal);
yval = vet_matrix(yval);
% digitar ycal e yval no pronpt da janela de comandos para visualizar a
% criação da matriz de classes.
%% 7) Usar a rotina PRETRAT para o pretamento dos espectros; Lembrar que não centrei na média
[xcal,xval]=pretrat(xcal,xval,{'deriv';[9,2,1]});
% [xcal,xval]=pretrat(xcal,xval,{'deriv';[9,2,1];'center'});
% Cortar o vetor de número de ondas de acordo com a seleção da faixa
%% 8) plotar as matrizes de calibração e validação preprocessadas para a
% visualização do preprocessamento espectral empregado.
figure(3)
subplot(211);
plot(numc,xcal)
axis tight
title('Amostras de calibração','fontname','Times New Roman','fontsize',14);
hold on
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(212);
plot(numc,xval)
axis tight
title('Amostras de validação','fontname','Times New Roman','fontsize',14); hold on
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.3,0.1,0.1],'String', 'B','LineStyle','none');
[ax,h1]=suplabel('Número de onda (cm^-^1)','x',[.11 .11 .84 .84]); 
[ax,h2]=suplabel('Primeira derivada','y',[.08 .08 .84 .84]); 
set(h1,'FontSize',24,'FontName','Times New Roman') 
set(h2,'FontSize',24,'FontName','Times New Roman') 
%% 9) validação Cruzada para determinar o Número de Variáveis Latentes
for u = 1:size(ycal,2)
[cvv(:,u)] = mycv(xcal,ycal(:,u),10,10);
end
figure(4)
for u = 1:size(ycal,2)
    plot(mean(cvv(u).perc_correc_class_val_cruzada,2));
    hold on
end
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
legend ('Classe 1', 'Classe 2', 'Classe 3', 'Classe 4')
xlabel('Número de variáveis latentes','FontSize',24,'FontName','Times New Roman');
ylabel('Porcentagem de amostras classificadas corretamente','FontSize',24,'FontName','Times New Roman');
legend({'Classe 1','Classe 2', 'Classe 3', 'Classe 4'},'FontSize',14,'FontName','Times New Roman')
disp(['Escolha o Número de variáveis latentes ','. ']);
A=input('');

%% 10) Calcular o modelo final PLS
[yprev_cal]=previsto_pls(xcal,ycal,xcal,0,A);
% Calculando o threshould
ts=[];
for ki=1:size(ycal,2)
    plsda_thres = plsdafindthr(yprev_cal(:,ki),ycal(:,ki));
    ts=[ts,plsda_thres.class_thr];
end
% Prevendo as amostras
for u = 1:size(ycal,2)
    yprev_calts(:,u)=yprev_cal(:,u)>=ts(:,u);
end
%% 11) Calcular o modelo de validação
[yprev_val]=previsto_pls(xcal,ycal,xval,0,A);
for u = 1:size(yval,2)
    yprev_valts(:,u)=yprev_val(:,u)>=ts(:,u);
end
%% 12) Plotar o gráfico das amostras de calibração
Nome{1,1}='Oliva' ;  Nome{1,2} = 'k';
Nome{2,1}='Canola' ; Nome{2,2} = 'm';
Nome{3,1}='Milho' ;  Nome{3,2} = 'r';
Nome{4,1}='Soja' ;   Nome{4,2} = 'g';
figure(5)
for ki=1:size(yprev_cal,2)
    %figure(ki)
    subplot(2,2,ki)
    tp1=find(ycal(:,ki));tp2=setxor(1:length(ycal),tp1);
    plot(1:length(yprev_cal),yprev_cal(:,ki),'o','MarkerFaceColor','b'), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    marc=strcat(Nome{ki,2},'o');
    plot(tp1,yprev_cal(tp1,ki),marc,'MarkerFaceColor',marc(:,1)), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
        hline(ts(ki),'r')
    title(Nome{ki,1},'FontSize',14,'FontName','Times New Woman')
    xlabel('Amostra','FontSize',18,'FontName','Times New Woman')
    ylabel(sprintf('Classe %g',ki),'FontSize',18,'FontName','Times New Woman')
    end
%% 13) Plotar o gráfico das amostras de validação
figure(6)
for ji=1:size(yprev_val,2)
    %figure(ji)
    subplot(2,2,ji)
    tv1=find(yval(:,ji));tv2=setxor(1:length(yval),tv1);
    plot(1:length(yprev_val),yprev_val(:,ji),'o','MarkerFaceColor','b'), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',14,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    marc=strcat(Nome{ji,2},'o');
    plot (tv1,yprev_val(tv1,ji),marc,'MarkerFaceColor',marc(:,1)), hold on
    set(gcf,'Color','white')
    set(gca,'xtick',[],'FontSize',12,'FontName','Times New Roman');
    set(gca,'xtickMode', 'auto');hold on;
    hline(ts(ji),'r')
    title(Nome{ji,1},'FontSize',14,'FontName','Times New Woman')
    xlabel('Amostra','FontSize',18,'FontName','Times New Woman')
    ylabel(sprintf('Classe %g',ji),'FontSize',18,'FontName','Times New Woman')
end
    
%% 14 Calculando a tabela de confusão para o conjunto de Calibração
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

C1=[class_param_cal.accuracy;class_param_cal.sensitivity(:,1);class_param_cal.specificity(:,1)]*100;
C2=[class_param_cal.accuracy;class_param_cal.sensitivity(:,2);class_param_cal.specificity(:,2)]*100;
C3=[class_param_cal.accuracy;class_param_cal.sensitivity(:,3);class_param_cal.specificity(:,3)]*100;
C4=[class_param_cal.accuracy;class_param_cal.sensitivity(:,4);class_param_cal.specificity(:,4)]*100;
figmerit={'Acurácia %','Sensibilidade %','Especificidade %'};
Tcal2=table(C1,C2,C3,C4,'RowNames',figmerit)
Tcal2.Properties.VariableNames = {'Azeite_de_Oliva','Oleo_de_Canola','Oleo_de_Milho','Oleo_de_Soja'}
writetable(Tcal2,'Resultados_PLSDA.xlsx','Sheet',1,'Range','A8','WriteRowNames',true)

%% 15 Calculando a tabela de confusão para o conjunto de Validação
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

C1v=[class_param_val.accuracy;class_param_val.sensitivity(:,1);class_param_val.specificity(:,1)]*100;
C2v=[class_param_val.accuracy;class_param_val.sensitivity(:,2);class_param_val.specificity(:,2)]*100;
C3v=[class_param_val.accuracy;class_param_val.sensitivity(:,3);class_param_val.specificity(:,3)]*100;
C4v=[class_param_val.accuracy;class_param_val.sensitivity(:,4);class_param_val.specificity(:,4)]*100;
figmerit={'Acurácia %','Sensibilidade %','Especificidade %'};
Tval2=table(C1v,C2v,C3v,C4v,'RowNames',figmerit)
Tval2.Properties.VariableNames = {'Azeite_de_Oliva', 'Oleo_de_Canola', 'Oleo_de_Milho', 'Oleo_de_Soja'}
writetable(Tval2,'Resultados_PLSDA.xlsx','Sheet',1,'Range','I8','WriteRowNames',true)

% clearvars -except class_param_cal class_param_val cvv num numc xcal xval ycal yval YCAL YVAL yprev_cal yprev_cal_new yprev_calts yprev_val yprev_val_new yprev_valts Tcal2 Tcal Tval2 Tval
