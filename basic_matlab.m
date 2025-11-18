
% =================== Variables ========================

a = sqrt(2)

b = 4 + 5

%.....latest printed
%      ans
%.....clear screen
clc
%=====================clear screen======================
clearvars



%==================== Variables datatype ===============
a = sqrt(2)

b = 4 + 5

whos

character = 'win'
stringvar = "hi am lucky"

whos

%========~ To view the data/ matrix click on variable square table======


%====================== Matrices and Vectors ========================


% 1.~~~~ vector or array is a 1xN matrix 
% middle is a step 
x = 1:10
x = 1:2:10
%~~~~~ Transpose
x'

%2.~~~~~ linspace(evenly spaced data(default 100)

y = linspace(0,1,6)

%3. ~~~~~ using list
list_vector_comma = [2,3,4,5]
list_vector_nocomma = [2 3 4 5]


%=====================Matrix=========================

A = [1 2; 3 4; 1 1]

% times to match we have to take transpose
A* A'

% element wise power
A .^ 2

% sqaure of ones
One = ones (3)

% nxm ones

one_n_m = ones(3,1)

zero = zeros(2)

% identity

i = eye(2)


%=================== Matrix Slicing ================
A = [1 2 0;-1 3 4; 1 0 7]

% count number
A(5)
% row column
A(2,2)

% end 
A(end)

% second last
A(end - 1)

A(2,1:3)

%===================== delete column ==================

%........ we have to update the modified 
A = A(1:2 , 1:3)


