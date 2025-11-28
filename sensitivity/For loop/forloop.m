%***** pause to click enter to see output of each iteration*********

x =[10 5 0 -5 -10];

for i=1:length(x)
    fprintf("entry %d of vector x is %d\n",i,x(i))
    %pause

end

%% ***** Copy vector into another vector*********

x=[1:8]';


y = 0;  %works but it is slower 
for i=1:length(x)
    %copy
    y(i,1)=x(i,1);
    fprintf("%d\n",y(i,1));
end
