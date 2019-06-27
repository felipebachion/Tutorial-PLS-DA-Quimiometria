function [Y] = matrix_vet(yvector)
classes=[1:size(yvector,2)];
Y=[];
for lol = 1:classes(end,end)
a=find(yvector(:,lol)==1);
Y1=[ones(size(a,1),1)*lol];
Y=[Y;Y1];
end
end

