function num = numConnComp(g)
    % Compute the Laplacian matrix
    L = laplacian(g);

    % Compute the eigenvectors and eigenvalues
    eigenvalues = eig(L);

    % Get the diagonal elements of the eigenmatrix
    % eigenvalues = diag(D);

    % Sort the eigenvalues in decreasing order
    % sortedEigenvalues = sort(eigenvalues, 'descend');

    % Count the number of elements less than 10^(-5)
    num = sum(eigenvalues < 10^(-5));
end

% function num = numConnComp(g)
%     num = length(unique(conncomp(g)));
% end