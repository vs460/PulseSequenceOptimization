function P = perm_mat_to_make_block_diag(A)
    % Make the undirected adjacency graph for A.
    nonzero = (A ~= 0);
    adjacency = nonzero + nonzero';
    G = graph(adjacency);

    % Perform breadth-first search.
    v_1 = 1;
    V_prime = bfsearch(G, v_1);
    n = length(A);
    if length(V_prime) == n
        error('Input is not permutation-similar to a block-diagonal matrix.');
    end

    % Make the permutation matrix.
    i = (1:n)';
    V_prime_complement = setdiff(i, V_prime);
    j = [V_prime; V_prime_complement];
    P = sparse(i, j, ones(n, 1), n, n);
end

