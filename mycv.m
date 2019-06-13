function [cvv] = mycv(X,y,comp,cv_groups);
% Rotina que realiza a validação cruzada selecionando aleatóriamente 80%
% das amostras de cada classe e com o restante das amostras realizada a
% validação cruzada.
% X é conjunto de calibracao
% y são os valores de classe da forma y = [1 0 0..
%                                          [0 1 0..]
% comp é o número máximo de variaveis que serão testadas
% cv_groups é o número de vezes que as amostras serão aleatórizadas para
% fazer a validaçaõ cruzada.

% Rotina desenvolvida por Felipe Bachion para fins não comerciais
% contato: felipebachion@gmail.com
rng('default') % Para reprodutibilidade
x=X;
assigned_class = []; 
    perc_in = 0.8;
    a1=find(y==1);
    a0=find(y==0);
    x1=x(a1,:); x0=x(a0,:);
    y1=y(a1,:); y0=y(a0,:);
    nobj1=size(y1,1);
    nobj0=size(y0,1);
    take_in1 = round(nobj1*perc_in);
    take_in0 = round(nobj0*perc_in);
    out_rand1 = zeros(nobj1,1);
    out_rand0= zeros(nobj0,1);
    class_true = [];
    ts=zeros(comp,cv_groups); % Pré-alocando a matriz
    n_correc_class_val_cruzada=[];
    h=waitbar(0,'Construindo os modelos PLS-DA');
for n=1:comp %1
       h=waitbar(n/comp);
    for i=1:cv_groups
        out1=ones(nobj1,1);
        out0=ones(nobj0,1);
        whos_in1 = randperm(nobj1);
        whos_in1 = whos_in1(1:take_in1);
        whos_in0 = randperm(nobj0);
        whos_in0 = whos_in0(1:take_in1);
        out1(whos_in1) = 0;
        out0(whos_in0) = 0;
        out_rand1(find(out1 == 1)) = out_rand1(find(out1 == 1)) + 1;
        out_rand0(find(out0 == 1)) = out_rand0(find(out0 == 1)) + 1;
        x_outt = x1(find(out1 == 1),:);
        x_outtt = x0(find(out0 == 1),:);
        x_out=[x_outt;x_outtt];
        y_outt = y1(find(out1 == 1));
        y_outtt = y0(find(out0 == 1));
        y_out=[y_outt;y_outtt];
        x_inn  = x1(whos_in1,:);
        x_innn  = x0(whos_in0,:);
        x_in=[x_inn;x_innn];
        y_inn  = y1(whos_in1,:);
        y_innn  = y0(whos_in0,:);
        y_in=[y_inn;y_innn];


        [yprev_val_cruzada]=previsto_pls(x_in,y_in,x_out,0,n);
        [yprev_cal_cruzada]=previsto_pls(x_in,y_in,x_in,0,n);
        cvv.pred_valcruzada(:,i)=yprev_val_cruzada;
        cvv.ref_valcruzada(:,i)=y_out;
        cvv.pred_calcruzada(:,i)=yprev_cal_cruzada;
        cvv.ref_calcruzada(:,i)=y_in;
         % calculo do threshould para cada vl nas amostras de validação
         % cruzada e já encontrando o número de amostras classificadas
         % corretamente
                for ki=1:size(cvv.pred_valcruzada,2)
                plsda_thres=plsdafindthr(cvv.pred_calcruzada(:,ki),cvv.ref_calcruzada(:,ki));
                ts (n,ki)=[plsda_thres.class_thr];
                     cv=yprev_val_cruzada;
                     cv=cv>=ts(n,ki);
                     cv=double(cv);
                     aw=cv==y_out;
                     n_correc_class_val_cruzada(n,ki)=sum(aw);
                end
     cvv.ts=ts;
     cvv.ts_mean=mean(ts,2)';
     cvv.perc_correc_class_val_cruzada=100*(n_correc_class_val_cruzada/size(y_out,1));
    end
   
end
    close (h)
end

