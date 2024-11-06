function [g] = grabgraphlets(G, L, n)

% This code takes a symmetric adjacency matric G and string of labels L to
% identify all unrooted labeled graphlets of size n. The output g is a cell
% array where each element contains a graphlet adjacency matrix `.G` and the
% string `.L`. These then need to be later permuted to identify automorphic
% structure of miminal "score" which is what we actually need to count.
%
% Predrag Radivojac
% Northeastern University
%
% December 13, 2023
% Boston, Massachusetts 02115

g = {};
k = 1;

% the graph must include the first node
if n == 1
    t = G(1, 1);
    g{1}.G = t;
    g{1}.L = L(1);
elseif n == 2
    no_nodes = size(G, 1);
    for i = 2 : no_nodes
        t = G([1 i], [1 i]);
        if numConnComp(graph(t)) == 1
            g{k}.G = t;
            g{k}.L = L([1 i]);
            k = k + 1;
        end
    end
elseif n == 3
    no_nodes = size(G, 1);
    for i = 2 : no_nodes - 1
        for ii = i + 1 : no_nodes
            t = G([1 i ii], [1 i ii]);
            if numConnComp(graph(t)) == 1
                g{k}.G = t;
                g{k}.L = L([1 i ii]);
                k = k + 1;
            end
        end
    end
elseif n == 4
    no_nodes = size(G, 1);
    for i = 2 : no_nodes - 2
        for ii = i + 1 : no_nodes - 1
            for iii = ii + 1 : no_nodes
                t = G([1 i ii iii], [1 i ii iii]);
                if numConnComp(graph(t)) == 1
                    g{k}.G = t;
                    g{k}.L = L([1 i ii iii]);
                    k = k + 1;
                end
            end
        end
    end
elseif n == 5
    no_nodes = size(G, 1);
    for i = 2 : no_nodes - 3
        for ii = i + 1 : no_nodes - 2
            for iii = ii + 1 : no_nodes - 1
                for iiii = iii + 1 : no_nodes
                    t = G([1 i ii iii iiii], [1 i ii iii iiii]);
                    if numConnComp(graph(t)) == 1
                        g{k}.G = t;
                        g{k}.L = L([1 i ii iii iiii]);
                        k = k + 1;
                    end
                end
            end
        end
    end
else
    error('Incorrect graphlet size.')
end

return
