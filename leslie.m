
L= [0.4 0.7 0.5 0.1;
    0.6 0 0 0;
    0 0.3 0 0;
    0 0 0.2 0;
   ];
L

%********** Eigen values
[Eigen_V,D]=eig(L)

Eigen_V

% initial pop =[100;40;50]

% ******* population in 2 years(X2)
X0=[100;40;50]

