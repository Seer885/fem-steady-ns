% MAIN FUNCTION
%
% num_nodes (Natural number, i.e. integer)
%       Number of nodes after 0 for x and y, e.g. 2 gives x,y={0 , 0.5 , 1}
%           Note this is assuming P1 Mesh (i.e. for element partitions)
% m (Natural number, i.e. integer)
%       Constant that affects how many previous terms are saved in A.A.
% mu (Real number greater than 0)
%       Viscosity constant needed to construct f

function lid_cavity(num_nodes, m, mu)
    tic;

    % ====================================================================
    %% 0 - Numerical functions for later reference
    % ====================================================================
    % Solution initial guess
    u0_x = @(x,y) 0;
    u0_y = @(x,y) 0;

    % P1 basis functions
    eta1 = @(x,y) x;
    eta2 = @(x,y) 1-x-y;
    eta3 = @(x,y) y;


    % P2 basis functions and gradients
    zeta1 = @(x,y) x .* (2*x - 1);
    zeta1i = @(x,y) 4 * x .* (1-x-y);
    zeta2 = @(x,y) (1-2*x-2*y) .* (1-x-y);
    zeta2i = @(x,y) 4 * (1-x-y) .* y;
    zeta3 = @(x,y) y .* (2*y - 1);
    zeta3i = @(x,y) 4 * x .* y;
    g_zeta1 = @(x,y) [4*x-1 ; 0];
    g_zeta1i = @(x,y) 4 * [1-2*x-y ; -x];
    g_zeta2 = @(x,y) [4*x+4*y-3 ; 4*x+4*y-3];
    g_zeta2i = @(x,y) 4 * [-y ; 1-x-2*y];
    g_zeta3 = @(x,y) [0 ; 4*y-1];
    g_zeta3i = @(x,y) 4 * [y ; x];


    % Assemble quadrature weights and nodes
    Qweights = (1/2) * [-9/16, 25/48, 25/48, 25/48];
    Qx_nodes = [1/3, 0.6, 0.2, 0.2];
    Qy_nodes = [1/3, 0.2, 0.6, 0.2];


    % Pre-compute outputs of Eta, Zeta, and Gradient Zeta at quadrature nodes
    QE = zeros(6,3);
    QZ = zeros(6,6);
    QGZ = zeros(12,6);
    for q = 1:4
        xq = Qx_nodes(q);
        yq = Qy_nodes(q);
        QE(q,:) = [eta1(xq,yq), eta2(xq,yq), eta3(xq,yq)];
        % Row is the current quadrature index, column is the appropriate
        % basis function (in order from 1, 1i, ..., m
        QZ(q,:) = [zeta1(xq,yq), zeta1i(xq,yq), zeta2(xq,yq), ...
            zeta2i(xq,yq), zeta3(xq,yq), zeta3i(xq,yq)];
        % Similar setup for gradients, except you jump 2 rows to switch to
        % a new quadrature index
        QGZ(2*q-1:2*q , :) = [g_zeta1(xq,yq), g_zeta1i(xq,yq), g_zeta2(xq,yq), ...
                        g_zeta2i(xq,yq), g_zeta3(xq,yq), g_zeta3i(xq,yq)];
    end



    % ====================================================================
    %% 1 - P1/P2 Setups
    % ====================================================================
    % Total number of rectangular elements for *both* P1 and P2
    num_elem = 2 * num_nodes^2;  
    num_nodes_all1 = num_nodes+1;  % Total nodes per variable in P1
    num_nodes_all2 = 2*num_nodes+1;  % Total nodes per variable in P2

    xy_vals2 = linspace( 0 , 1 , num_nodes_all2);  % List of all distinct values in P2
    [X2,Y2] = meshgrid(xy_vals2 , xy_vals2);  % Set up mesh grid for P2

    % Assign a unique global index to each unique (i,j) in the above meshes
    % Will count up from top to bottom and then from left to right
    total_nodes1 = (num_nodes_all1)^2;  % Total nodes period in P1
    total_nodes2 = (num_nodes_all2)^2;  % Total nodes period in P2
    glob_idx1 = reshape( linspace(1, total_nodes1, total_nodes1), ...
        [num_nodes_all1 , num_nodes_all1] );  % P1
    glob_idx2 = reshape( linspace(1, total_nodes2, total_nodes2), ...
        [num_nodes_all2 , num_nodes_all2] );  % P2


    % Create local connectivity arrays
    connections1 = zeros(3 , num_elem);
    connections2 = zeros(6 , num_elem);
    
    % P1: (1,2,3) --> (g1,g2,g3)
    elem = 0;
    for i = 1:num_nodes
        for j = 1:num_nodes
            % Starting at (xi,yj), 00 > 10 > 11
            elem = elem + 1;
            connections1(:,elem) = [glob_idx1(i , j); glob_idx1(i , j+1); glob_idx1(i+1 , j+1)];

            % Starting at (xi,yj), 00 > 01 > 11
            elem = elem + 1;
            connections1(:,elem) = [glob_idx1(i , j); glob_idx1(i+1 , j); glob_idx1(i+1 , j+1)];
        end
    end

    % P2: (1,1i,2,2i,3,3i) --> (g1,g2,g3,g4,g5,g6)
    elem = 0;
    for i = 1 : 2 : (num_nodes_all2 - 2)  % Set up P2
        for j = 1 : 2 : (num_nodes_all2 - 2)
            % 00 > 10 > 20 > 21 > 22 > 11
            elem = elem + 1;
            connections2(:,elem) = [glob_idx2(i , j); glob_idx2(i , j+1); 
                glob_idx2(i , j+2); glob_idx2(i+1 , j+2); 
                glob_idx2(i+2 , j+2); glob_idx2(i+1 , j+1)];

            % 00 > 01 > 02 > 12 > 22 > 11
            elem = elem + 1;
            connections2(:,elem) = [glob_idx2(i , j); glob_idx2(i+1 , j); 
                glob_idx2(i+2 , j); glob_idx2(i+2 , j+1); 
                glob_idx2(i+2 , j+2); glob_idx2(i+1 , j+1)];
        end
    end
    



    % ====================================================================
    %% 2 - Set up constant stiffness matrices and vectors
    % ====================================================================
    % Set up our h
    h = 1 / num_nodes;

    % Recurring arrays for later assembly
    D1 = sparse(total_nodes1,total_nodes1);
    D2 = sparse(total_nodes2,total_nodes2);
    


    % We have 5 global stiffness matrices and 2 global load vectors:
    %          A, B, Cx_odd/Cx_even, and Cy_odd/Cy_even
    % The local stiffness matrices for A, Cx, and Cy are consistent thanks to our uniform distribution, so we can define them as so:
    lA = (1/6) * [3,-4,1,0,0,0;-4,16,-4,0,0,-8;1,-4,6,-4,1,0;0,0,-4,16,-4,-8;0,0,1,-4,3,0;0,-8,0,-8,0,16];
    lCx_odd = (h/6) * [1,0,0;-1,1,0;0,-1,0;-1,-1,-2;0,0,0;1,1,2];
    lCx_even = (h/6) * [0,0,0;2,1,1;0,1,0;0,-1,1;0,0,-1;-2,-1,-1];
    lCy_odd = lCx_even;
    lCy_even = lCx_odd;

    A = sparse(total_nodes2,total_nodes2);
    Cx = sparse(total_nodes2,total_nodes1);   Cy = Cx;

    for elem = 1:num_elem
        % Obtain P1 and P2 global indices
        g_ind1 = connections1(: , elem);
        g_ind2 = connections2(: , elem);

        % Update values of global stiffness matrices and load vectors
        A(g_ind2,g_ind2) = A(g_ind2,g_ind2) + lA;
        if rem(elem,2) == 1
            Cx(g_ind2,g_ind1) = Cx(g_ind2,g_ind1) + lCx_odd;
            Cy(g_ind2,g_ind1) = Cy(g_ind2,g_ind1) + lCy_odd;
        else
            Cx(g_ind2,g_ind1) = Cx(g_ind2,g_ind1) + lCx_even;
            Cy(g_ind2,g_ind1) = Cy(g_ind2,g_ind1) + lCy_even;
        end
    end


    % Create 1D array of all P2 global indices on the boundary
    % HOWEVER, using the these indices alone would only count for u_x's
    % boundary conditions, so to also factor in u_y's conditions, we
    % would need to the same indices translated the length of A+B,
    % which happens to be num_nodes_all2.
    bound_idx = [glob_idx2( 1 , : ) , ...
        glob_idx2( 1 , : ) + total_nodes2 , ...
        glob_idx2( num_nodes_all2 , : ) , ...
        glob_idx2( num_nodes_all2 , : ) + total_nodes2 , ...
        glob_idx2( : , 1 )' , ...
        glob_idx2( : , 1 )' + total_nodes2 , ...
        glob_idx2( : , num_nodes_all2 )' , ...
        glob_idx2( : , num_nodes_all2 )' + total_nodes2];
    top_idx = glob_idx2( num_nodes_all2 , : );

    % Set up vector of appropriate weights for the discrete mean sum later
    weights = h^2 * ones(1,total_nodes1);
    % Edge weights need to be divided by 2
    edge_idx = [glob_idx1( 1 , : ) , glob_idx1( num_nodes_all1 , : ) , glob_idx1( : , 1 )' , glob_idx1( : , num_nodes_all1 )'];
    weights(edge_idx) = weights(edge_idx) / 2;
    % Corner weights need to be divided by 4, though since they're also
    % edge nodes, they were already divided by 2 once, so we just divide by
    % 2 again once we've isolated them.
    corner_idx = [glob_idx1(1,1) , glob_idx1(1,num_nodes_all1), ...
        glob_idx1(num_nodes_all1,1),glob_idx1(num_nodes_all1,num_nodes_all1)];
    weights(corner_idx) = weights(corner_idx) / 2;




    % ====================================================================
    %% 3 - Get u0, u1, and G(u1) set up for Anderson Acceleration
    % ====================================================================
    % Create ux_0 and uy_0 arrays for our iteration
    uk_x = zeros(total_nodes2,1);   uk_y = uk_x;
    for i = 1:num_nodes_all2
        for j = 1:num_nodes_all2
            uk_x((i-1)*num_nodes_all2+j) = u0_x(xy_vals2(i),xy_vals2(j));
            uk_y((i-1)*num_nodes_all2+j) = u0_y(xy_vals2(i),xy_vals2(j));
        end
    end
    % Combine them into a single array
    uk = [uk_x ; uk_y];

    % Using G(u) as shorthand for our process that takes u(k-1) and outputs
    % u(k), before we can even start our Anderson iteration process, we
    % need to store the following arrays (given u1 = G(u0)):
    %    stored_uk = [u0 , u1]  
    %    stored_res = [u1 - u0 , G(u1) - u1]
    stored_uk = [];  stored_res = [];

    for k = 1:2
        stored_uk = [stored_uk , uk];

        B = sparse(total_nodes2,total_nodes2);
        for elem = 1:num_elem
            % Set up local stiffness matrices for B
            lB = sparse(zeros(6,6)); 
    
            % Explicitly set up functions for ux_0 and uy_0 based on
            % the current element K, which given global nodes n1,n1i, is:
            %   u(n1) * zeta_1 + u(n1i) * zeta_1i + u(n2) * zeta_2 +...
            UB_x = QZ * uk_x(connections2(:,elem));
            UB_y = QZ * uk_y(connections2(:,elem));
    
    
            % Set up B matrix using Gaussian Quadratures
            if rem(elem,2) == 1
                for q = 1:4  % Quadrature index for weight/node
                    % Note we need [-z_x,z_y]
                    lB(:,:) = lB(:,:) + Qweights(q) * QZ(q,:)' * [UB_x(q) , UB_y(q)] * (QGZ([2*q-1 , 2*q] , :) .* [-1;1]);
                end       
            else
                for q = 1:4  % Quadrature index for weight/node
                    % Note we need [z_y,-z_x]
                    lB(:,:) = lB(:,:) + Qweights(q) * QZ(q,:)' * [UB_x(q) , UB_y(q)] * (QGZ([2*q , 2*q-1] , :) .* [1;-1]);
                end
            end 
            lB = h * lB;
    
    
            % Obtain P1 and P2 global indices
            g_ind2 = connections2(: , elem);
    
            % Update values of the global stiffness matrices for B
            B(g_ind2,g_ind2) = B(g_ind2,g_ind2) + lB;
        end
    
    
        % Combined stiffness matrix K will be of the form:
        K = [mu*A+B ,  D2, Cx;
            D2 , mu*A+B , Cy;
            Cx' , Cy' , D1];
        % Meanwhile the combined load vector g will be:
        g = zeros(2*total_nodes2+total_nodes1, 1);
    
    
        % Set all rows corresponding to bound indices to 0, and
        % then set corresponding rows of g for top indices for ux to 1
        %K(bound_idx , :) = 0;   K(: , bound_idx) = 0;  g(top_idx) = 1;
        K(bound_idx , :) = 0;  g(top_idx) = 1;
        % Update corresponding diagonals of K to 1
        for idx = bound_idx(:)'
            K(idx , idx) = 1;
        end
    
        % Adjustment to K to enforce the mean of ph needing to equal 0
        last = 2*total_nodes2 + total_nodes1;
        K(last,1:2*total_nodes2) = zeros(1,2*total_nodes2);
        % Using the above setup, we enforce sum(ph_i * w_i) = 0.
        K(last,2*total_nodes2+1:last) = weights;  g(last) = 0;
    
    
        % Solve the system Kz=g.
        z = K \ g;
        uk_x = z(1:total_nodes2);  uk_y = z(total_nodes2+1:2*total_nodes2);
        uk = z(1:2*total_nodes2);

        stored_res = [stored_res , uk - stored_uk(:,k)];  
    end




    % ====================================================================
    %% 4 - Run Anderson Acceleration
    % ====================================================================
    % Set up tolerance, maximum iterations, and constant m
    ToL = 1e-6;   max_Iter = 1000;

    % To keep our computations succinct, Anderson Acceleration has us
    % create special arrays Resk and Uk.
    % First we'll set up the initial versions of these arrays.
    Resk = stored_res(:,2) - stored_res(:,1);  
    Uk = stored_uk(:,2) - stored_uk(:,1);

    % Run the following iteration until either the maximum number of
    % iterations has been reached, or both of our most recent residuals are
    % less than our specified tolerance
    while k < (max_Iter) && norm( stored_res(:,k))  > ToL
        mk = min(m,k);  % Set up mk for the current loop

        % Solve our two least-squares problems to get gamma_k
        [QQ,RR] = qr(Resk, "econ");  % QR-factorize the current Resk 
        % Solve for our gamma using Cholesky
        gamma_k = RR \ (QQ' * stored_res(:,k));  

        % Find and store the actual uk using our AA formula
        stored_uk(:,k+1) = stored_uk(:,k) + stored_res(:,k) - (Uk + Resk) * gamma_k;
        
        % Run our FEM solver to get G(uk_x) and G(uk_y) to update our residual matrix
        uk_x = stored_uk(1:total_nodes2,k+1);
        uk_y = stored_uk(total_nodes2+1:2*total_nodes2,k+1);
        B = sparse(total_nodes2,total_nodes2);
        for elem = 1:num_elem
            % Set up local stiffness matrices for B
            lB = sparse(zeros(6,6)); 
    
            % Explicitly set up functions for ux_0 and uy_0 based on
            % the current element K, which given global nodes n1,n1i, is:
            %   u(n1) * zeta_1 + u(n1i) * zeta_1i + u(n2) * zeta_2 +...
            UB_x = QZ * uk_x(connections2(:,elem));
            UB_y = QZ * uk_y(connections2(:,elem));
    
    
            % Set up B matrix using Gaussian Quadratures
            if rem(elem,2) == 1
                for q = 1:4  % Quadrature index for weight/node
                    % Note we need [-z_x,z_y]
                    lB(:,:) = lB(:,:) + Qweights(q) * QZ(q,:)' * [UB_x(q) , UB_y(q)] * (QGZ([2*q-1 , 2*q] , :) .* [-1;1]);
                end       
            else
                for q = 1:4  % Quadrature index for weight/node
                    % Note we need [z_y,-z_x]
                    lB(:,:) = lB(:,:) + Qweights(q) * QZ(q,:)' * [UB_x(q) , UB_y(q)] * (QGZ([2*q , 2*q-1] , :) .* [1;-1]);
                end
            end 
            lB = h * lB;
    
    
            % Obtain P1 and P2 global indices
            g_ind2 = connections2(: , elem);
    
            % Update values of the global stiffness matrices for B
            B(g_ind2,g_ind2) = B(g_ind2,g_ind2) + lB;
        end
    
    
        % Combined stiffness matrix K will be of the form:
        K = [mu*A+B ,  D2, Cx;
            D2 , mu*A+B , Cy;
            Cx' , Cy' , D1];
        % Meanwhile the combined load vector g will be:
        g = zeros(2*total_nodes2+total_nodes1, 1);
    
    
        % Set all rows and columns corresponding to bound indices to 0, and
        % then set corresponding rows of g for top indices for ux to 1
        %K(bound_idx , :) = 0;   K(: , bound_idx) = 0;  g(top_idx) = 1;
        K(bound_idx , :) = 0;  g(top_idx) = 1;
        % Update corresponding diagonals of K to 1
        for idx = bound_idx(:)'
            K(idx , idx) = 1;
        end
    
        % Adjustment to K to enforce the mean of ph needing to equal 0
        last = 2*total_nodes2 + total_nodes1;
        K(last,1:2*total_nodes2) = zeros(1,2*total_nodes2);
        % Using the above setup, we enforce sum(ph_i * w_i) = 0.
        K(last,2*total_nodes2+1:last) = weights;  g(last) = 0;
    
    
        % Solve the system Kz=g, and save G(u_{k-1})
        z = K \ g;
        Gu = z(1:2*total_nodes2);


        % Update our residual arrays with this result
        stored_res = [stored_res , Gu - stored_uk(:,k+1)];  

        % Update Resk and Uk arrays
        Resk = [Resk , stored_res(:,k+1) - stored_res(:,k)]; 
        Uk = [Uk , stored_uk(:,k+1) - stored_uk(:,k)];

        % However, the minimum index of both Resk and Uk is k-mk, so we 
        % must trim off values as needed
        n = size(Resk, 2);
        if n > mk
            Resk = Resk(: , n-mk+1 : end);
            Uk = Uk(: , n-mk+1 : end);
        end

        % Increment k for the next iteration
        k = k + 1;
    end
    fprintf('Converged at iteration %d\n', k);
    
    uh_x = z(1:total_nodes2);  uh_y = z(total_nodes2+1:2*total_nodes2);
    %ph = z(2*total_nodes2+1:length(z));

    UX = reshape( uh_x, [num_nodes_all2 , num_nodes_all2] );
    UY = reshape( uh_y, [num_nodes_all2 , num_nodes_all2] );

    toc  % Output runtime
    fprintf('\n')
    
    % Set starting locations for paths in upcoming streamline plot
    [start_X , start_Y] = meshgrid([linspace(0,0.1,3) , linspace(0.1,0.9,8) , linspace(0.9,1,10)] , ...
                                    [linspace(0,0.1,3) , linspace(0.1,1,10)]);

    % Generate streamline plot to visualize our vector field
    streamline(X2,Y2,UX,UY,start_X,start_Y, color="white")

    %semilogy(vecnorm(stored_res, 1, 1), 'DisplayName', sprintf('Picard Anderson (%d)', m));
end