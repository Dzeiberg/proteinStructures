function [countvector, graphlets] = countgraphlets (G, L, ns, alphabet)

% This function takes a graph adjacency matric G and a string of node
% labels L to then provide a count vector counts that can be used as
% feature vector. The output element graphlets contains the same thing as
% counts but it does give each graph as well as the label string. If we
% only want a feature vector, we don't need graphlets, we just need counts.
%
% Predrag Radivojac
% Northeastern University
%
% December 13, 2023
% Boston, Massachusetts 02115

global BaseGraphs

% Cell array where each element is one graphlet with its matrix (G), label
% sequence (L), index of sorted labels (S) and type of graphlet (B). 
% Graphlets observed more than once will appear here as duplicates.
graphlets = {};
% TO DO: Pre-allocate graphlets

% Sparse feature vector of counts. If a graphlet is observed more than
% once, the countvector will have an appropriate count for that graphlet
countvector = [];

for n = ns
    P = perms(1 : n);
    no_perm = size(P, 1);

    for i = 1 : size(G, 1)
        % distances from node i to nodes with index at least i
        d = distances(graph(G), i, i : size(G, 1));
        
        % find nodes where the shortest path from i is less than n
        q = find(d < n) + i - 1;

        % find graphlets using indices q, but they are not "oriented" yet
        g = grabgraphlets(G(q, q), L(q), n);

        % now, orient those graphlets
        for j = 1 : length(g)
            edges = sum(sum(g{j}.G));

            minscore = 1e100; % something really large

            % graphlet to be added
            next = length(graphlets) + 1;

            for k = 1 : length(BaseGraphs{n})
                bg = BaseGraphs{n}{k}.G;
                %bps = BaseGraphs{n}{k}.P;
                if edges == sum(sum(bg))
                    %for l = 1 : size(bps, 3)
                    for l = 1 : no_perm
                        %t = bps(:, :, l)' * bg * bps(:, :, l);
                        R = eye(n);
                        R = R(:, [P(l, :)]);
                        %t = R' * bg * R;
                        t = R * bg * R';
                        if sum(sum(abs(g{j}.G - t))) == 0
                            
                            % DZ : Pre-allocate prm
                            % -----
                            prm = [];
                            % -----
                            %prm = zeros(1 , size(bg, 1));
                            % END DZ : Pre-allocate prm

                            for m = 1 : size(bg, 1)
                                %prm(m) = find(bps(m, :, l) == 1);
                                %prm(m) = find(R(m, :) == 1);
                                prm(m) = find(R(:, m) == 1);
                            end

                            %prm = BaseGraphs{n}{k}.p(:, :, l);
                            label = g{j}.L(prm);

                            %lbl = g{j}.l(prm);

                            % something wrong here with generating p
                            %label = g{j}.L(BaseGraphs{n}{k}.p(:, :, l));
                            score = getscore(label, alphabet);
                            %score = getscore(lbl, alphabet);

                            % if graph is identical, and label minimized
                            if score < minscore
                                minscore = score;
                                %graphlets{next}.G = t;
                                graphlets{next}.G = bg; % not needed
                                graphlets{next}.L = label;
                                graphlets{next}.S = score;
                                graphlets{next}.B = k;
                            end
                        end
                    end
                end
            end
            if minscore == 1e100
                error('Problem with counting');
            end
        end
    end

    % v = sparse(1, length(BaseGraphs{n}) * length(alphabet) ^ n);
    % sparse_indices = zeros(length(graphlets));
    sparse_indices = [];
    for i = 1 : length(graphlets)
        position = graphlets{i}.S + (graphlets{i}.B - 1) * length(alphabet) ^ n;
        sparse_indices = [sparse_indices position];
    end
    countvector = [countvector sparse(1, sparse_indices, 1, 1, length(BaseGraphs{n}) * length(alphabet) ^ n)];
end

return
