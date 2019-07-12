function [x_inn,y_inn,x_outt,y_outt] = myPerm (X,y,n_class,perc_in)
% Rotina que seleciona aleatóriamente 80% das amostras de cada classe para realizar a validação cruzada
% X é conjunto de calibracao
% y são os valores de classe da forma y = [1 0 0..
%                                          [0 1 0..]
% n_class = o numero de classes
% perc_in = 100% - a porcentagem de amostras que serão utilizadas na
% validação cruzada.

% Rotina desenvolvida por Felipe Bachion para fins não comerciais
% contato: felipebachion@gmail.com
if size(y,2)~=n_class
f = msgbox('O número de classes informado (n_class) não é igual ao número de classes do vetor y');
end
nobj1=[0];
for p=1:n_class
    classes=find(y(:,p)==1);
    nobj=size(classes,1);
    nobj1=[nobj+nobj1];
    take_in = round(nobj*perc_in);
    whos_in = randperm(nobj);
    whos_in=whos_in+(nobj1-nobj);
    whos_in1 = whos_in(1:take_in);
    x_in1=X(whos_in1,:);
    y_in1=y(whos_in1,:);
    x_out1=X(whos_in(1,take_in+1:end),:);
    y_out1=y(whos_in(1,take_in+1:end),:);
    x_in0{p}=x_in1;
    y_in0{p}=y_in1;
    x_out0{p}=x_out1;
    y_out0{p}=y_out1;
end
x_inn=[];
y_inn=[];
x_outt=[];
y_outt=[];

for i=1:n_class
x_in=[x_in0{i}];
x_inn=[x_inn;x_in];

y_in=[y_in0{i}];
y_inn=[y_inn;y_in];

x_out=[x_out0{i}];
x_outt=[x_outt;x_out];

y_out=[y_out0{i}];
y_outt=[y_outt;y_out];

end