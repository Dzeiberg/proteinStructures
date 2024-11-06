clear
% clc

global BaseGraphs

% alphabet = 'ACDEFGHIKLMNPQRSTVWY';

alphabet = 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy';

load UnrootedGraphlets.mat
load phosY
% G = read_svml('data/P04637/residues.adj');
% L = fileread('data/P04637/wildtype.labels');
disp(size(L))
% graphlets to count
%ns = [1 2 3 4 5];
ns = [1 2 3 4 5];

tic
[v, g] = countgraphlets(G, L, ns, alphabet);
toc
