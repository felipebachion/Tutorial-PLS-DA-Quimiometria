%% Rotina desenvolvida em função do Artigo publicado no periodico Quimica nova
% Fazer o donwload de todas as funções necessárias em:
% https://github.com/felipebachion/Tutorial-PLS-DA-Quimiometria
% As funções e o banco de amostras também podem ser solicitados via e-mail:
% rjpoppi@unicamp.br

%% 1) carregar o conjunto de dados
load amostras
%% 2) plotar o cojunto completo de espectros
figure(1)
subplot(2,1,1)
plot(num,X);
axis tight
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(2,1,2)
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

%% 8) plotar as matrizes de calibração e validação preprocessadas para a
% visualização do preprocessamento espectral empregado.
figure(3)
subplot(2,1,1);
plot(numc,xcal)
axis tight
title('Amostras de calibração','fontname','Times New Roman','fontsize',14);
hold on
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
    set(gca,'xtickMode', 'auto');hold on;
annotation('textbox', [0.8,0.8,0.1,0.1],'String', 'A','LineStyle','none');
subplot(2,1,2);
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
cvvc = my_cross_validation(xcal,ycal,10,10,4,0.8);

figure(4)
plot(cvvc.porc_am_class_cor')
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
legend ('Classe 1', 'Classe 2', 'Classe 3', 'Classe 4')
xlabel('Número de variáveis latentes','FontSize',24,'FontName','Times New Roman');
ylabel('Porcentagem de amostras classificadas corretamente','FontSize',24,'FontName','Times New Roman');
disp(['Escolha o Número de variáveis latentes ','. ']);
A=input('');

%% 10) Analisando a presença de amostras anômalas para o número de variaveis latentes escolhido
[model] = my_calc_qt_limits_cal (xcal,ycal,xval,A);
% Amostras anômalas
am_anomalasQ=find(model.Qres>=model.qlim);
am_anomalasT2=find(model.Thot>=model.tlim);
a=ismember(am_anomalasQ,am_anomalasT2);
am_anomalas=am_anomalasQ(1,a);

%% 11) Eliminando estas amostras e selecionando novamente o número ideal de variáveis latentes.
xcal(am_anomalas,:)=[];
ycal(am_anomalas,:)=[];

cvvc = my_cross_validation(xcal,ycal,10,10,4,0.8);

figure
plot(cvvc.porc_am_class_cor')
set(gca,'xtick',[],'FontSize',18,'FontName','Times New Roman');...
set(gca,'xtickMode', 'auto');hold on;
legend ('Classe 1', 'Classe 2', 'Classe 3', 'Classe 4')
xlabel('Número de variáveis latentes','FontSize',24,'FontName','Times New Roman');
ylabel('Porcentagem de amostras classificadas corretamente','FontSize',24,'FontName','Times New Roman');
disp(['Escolha o Número de variáveis latentes ','. ']);
A=input('');

%% 12) Analisando novamente a presença de amostras anômalas para o novo modelo PLS-DA nos conjuntos de calibração e validação.
[model_cal] = my_calc_qt_limits_cal (xcal,ycal,xval,A);
[model_val] = my_calc_qt_limits_val (xcal,ycal,xval,A);
% Desta vez não temos amostras anômalas no conjunto de calibração e
% validação

%% 13) Calcular o modelo final PLS-DA
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
%% 14) Calcular os valores de previsão das amostras de validação
[yprev_val]=previsto_pls(xcal,ycal,xval,0,A);
for u = 1:size(yval,2)
    yprev_valts(:,u)=yprev_val(:,u)>=ts(:,u);
end

%% 15) Plotar o gráfico das amostras de calibração
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
%% 16) Plotar o gráfico das amostras de validação
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
    
%% 17 Calculando a tabela de confusão para o conjunto de Calibração

[Tcal Tcal2] = my_ConfTable_cal (ycal,yprev_cal,ts)
[Tval Tval2] = my_ConfTable_val (yval,yprev_val,ts)

clearvars -except model_cal model_val class_param_cal class_param_val num numc xcal xval ycal yval YCAL YVAL yprev_cal yprev_cal_new yprev_calts yprev_val yprev_val_new yprev_valts Tcal2 Tcal Tval2 Tval
