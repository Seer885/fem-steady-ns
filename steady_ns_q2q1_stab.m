% MAIN FUNCTION
%
% num_nodes (Natural number, i.e. integer)
%       Number of nodes after 0 for x and y, e.g. 2 gives x,y={0 , 0.5 , 1}
%           Note this is assuming Q1 Mesh (i.e. for element partitions)
% mu (Real number greater than 0)
%       Viscosity constant needed to construct f

function steady_ns_q2q1_stab(num_nodes,mu)
    tic;

    % ====================================================================
    %% 0 - Numerical functions for later reference
    % ====================================================================
    % Solutions, external force f, and velocity field beta
    u_x = @(x,y) pi * sin(pi * x).^2 .* sin(2 * pi * y);
    u_y = @(x,y) -pi * sin(2 * pi * x) .* sin(pi * y).^2;
    u_x_dx = @(x,y) pi^2 * sin(2 * pi * x) .* sin(2 * pi * y);
    u_x_dy = @(x,y) 2 * pi^2 * sin(pi * x) .^ 2 .*  cos(2 * pi * y);
    u_y_dx = @(x,y) -2 * pi^2 * cos(2 * pi * x) .* sin(pi * y) .^ 2;
    u_y_dy = @(x,y) -pi^2 * sin(2 * pi * x) .* sin(2 * pi * y);

    p = @(x,y) sin(pi * x) .* sin(pi * y) - (4 / pi^2);
    p_dx = @(x,y) pi * cos(pi * x) .* sin(pi * y);
    p_dy = @(x,y) pi * sin(pi * x) .* cos(pi * y);
 
    f_x = @(x,y) mu * 2 * pi^3 * (1-2 * cos(2 * pi * x)) .* sin(2 * pi * y) + u_x(x,y) * u_x_dx(x,y) + u_y(x,y) * u_x_dy(x,y) + p_dx(x,y);
    f_y = @(x,y) mu * 2 * pi^3 * sin(2 * pi * x) .* (2 * cos(2 * pi * y)-1) + u_x(x,y) * u_y_dx(x,y) + u_y(x,y) * u_y_dy(x,y) + p_dy(x,y);

    u0_x = @(x,y) sin(pi * x).^2 .* sin(2 * pi * y);
    u0_y = @(x,y) -sin(2 * pi * x) .* sin(pi * y).^2;


    % Q1 basis functions and gradients
    eta1 = @(x,y) (1-x) .* (1-y);
    eta2 = @(x,y) x .* (1-y);
    eta3 = @(x,y) x .* y;
    eta4 = @(x,y) (1-x) .* y;
    g_eta1 = @(x,y) [y-1 ; x-1];
    g_eta2 = @(x,y) [1-y ; -x];
    g_eta3 = @(x,y) [y ; x];
    g_eta4 = @(x,y) [-y ; 1-x];


    % Q2 basis functions and gradients, constructed from zero lines
    z1 = @(x) x;   z2 = @(x) 1-2*x;   z3 = @(x) 1-x;
    z4 = @(y) y;   z5 = @(y) 1-2*y;   z6 = @(y) 1-y;
    % Basic functions
    zeta1 = @(x,y)       z2(x) .* z3(x) .* z5(y) .* z6(y);
    zeta1i = @(x,y) 4 *  z1(x) .* z3(x) .* z5(y) .* z6(y);
    zeta2 = @(x,y)  -    z1(x) .* z2(x) .* z5(y) .* z6(y);
    zeta2i = @(x,y) -4*  z1(x) .* z2(x) .* z4(y) .* z6(y);
    zeta3 = @(x,y)       z1(x) .* z2(x) .* z4(y) .* z5(y);
    zeta3i = @(x,y) -4 * z1(x) .* z3(x) .* z4(y) .* z5(y);
    zeta4 = @(x,y)  -    z2(x) .* z3(x) .* z4(y) .* z5(y);
    zeta4i = @(x,y) 4 *  z2(x) .* z3(x) .* z4(y) .* z6(y);
    zetam = @(x,y) 16 *  z1(x) .* z3(x) .* z4(y) .* z6(y);
    % Gradients of basis functions
    g_zeta1 = @(x,y)       [z5(y) .* z6(y) .* (-2*z3(x) - z2(x)) ; z2(x) .* z3(x) .* (-2*z6(y) - z5(y))];
    g_zeta1i = @(x,y) 4 *  [z5(y) .* z6(y) .* (z3(x) - z1(x)) ; z1(x) .* z3(x) .* (-2*z6(y) - z5(y))];
    g_zeta2 = @(x,y)  -    [z5(y) .* z6(y) .* (z2(x) -2* z1(x)) ; z1(x) .* z2(x) .* (-2*z6(y) - z5(y))];
    g_zeta2i = @(x,y) -4 * [z4(y) .* z6(y) .* (z2(x) -2* z1(x)) ; z1(x) .* z2(x) .* (z6(y) - z4(y))];
    g_zeta3 = @(x,y)       [z4(y) .* z5(y) .* (z2(x) -2* z1(x)) ; z1(x) .* z2(x) .* (z5(y) -2* z4(y))];
    g_zeta3i = @(x,y) -4 * [z4(y) .* z5(y) .* (z3(x) - z1(x)) ; z1(x) .* z3(x) .* (z5(y) -2* z4(y))];
    g_zeta4 = @(x,y)  -    [z4(y) .* z5(y) .* (-2*z3(x) - z2(x)) ; z2(x) .* z3(x) .* (z5(y) -2* z4(y))];
    g_zeta4i = @(x,y) 4 *  [z4(y) .* z6(y) .* (-2*z3(x) - z2(x)) ; z2(x) .* z3(x) .* (z6(y) - z4(y))];
    g_zetam = @(x,y)  16 * [z4(y) .* z6(y) .* (z3(x) - z1(x)) ; z1(x) .* z3(x) .* (z6(y) - z4(y))];
    % Laplacians of basis functions
    l_zeta1 = @(x,y) 4 * ( z5(y) .* z6(y) + z2(x) .* z3(x) );
    l_zeta1i = @(x,y) 8 * ( 2 * z1(x) .* z3(x) - z5(y) .* z6(y) );
    l_zeta2 = @(x,y) 4 * ( z5(y) .* z6(y) - z1(x) .* z2(x) );
    l_zeta2i = @(x,y) 8 * ( 2 * z4(y) .* z6(y) + z1(x) .* z2(x) );
    l_zeta3 = @(x,y) -4 * ( z4(y) .* z5(y) + z1(x) .* z2(x) );
    l_zeta3i = @(x,y) 8 * ( z4(y) .* z5(y) + 2 * z1(x) .* z3(x) );
    l_zeta4 = @(x,y) 4 * ( z2(x) .* z3(x) - z4(y) .* z5(y) );
    l_zeta4i = @(x,y) 8 * ( 2 * z4(y) .* z6(y) - z2(x) .* z3(x) );
    l_zetam = @(x,y) -32 * ( z4(y) .* z6(y) + z1(x) .* z3(x) );



    % Recurring floats saved for optimization
    a = 25/324;   b = 10/81;   c = (1-sqrt(3/5)) / 2;   d = (1+sqrt(3/5)) / 2;
    % Assemble quadrature weights and nodes
    Qweights = [a,b,a,b, 16/81 ,b,a,b,a];
    Qx_nodes = [c,c,c,1/2,1/2,1/2,d,d,d];
    Qy_nodes = [c,1/2,d,c,1/2,d,c,1/2,d];

    % Pre-compute outputs of Eta, Zeta, and Gradient Zeta at quadrature nodes
    QE = zeros(9,4);
    QZ = zeros(9,9);
    QLZ = zeros(9,9);
    QGE = zeros(18,4);
    QGZ = zeros(18,9);
    for q = 1:9
        xq = Qx_nodes(q);
        yq = Qy_nodes(q);
        QE(q,:) = [eta1(xq,yq), eta2(xq,yq), eta3(xq,yq), eta4(xq,yq)];
        % Row is the current quadrature index, column is the appropriate
        % basis function (in order from 1, 1i, ..., m
        QZ(q,:) = [zeta1(xq,yq), zeta1i(xq,yq), zeta2(xq,yq), ...
            zeta2i(xq,yq), zeta3(xq,yq), zeta3i(xq,yq), ...
            zeta4(xq,yq), zeta4i(xq,yq), zetam(xq,yq)];
        QLZ(q,:) = [l_zeta1(xq,yq), l_zeta1i(xq,yq), l_zeta2(xq,yq), ...
                    l_zeta2i(xq,yq), l_zeta3(xq,yq), l_zeta3i(xq,yq), ...
                    l_zeta4(xq,yq), l_zeta4i(xq,yq), l_zetam(xq,yq)];
        % Similar setup for gradients, except you jump 2 rows to switch to
        % a new quadrature index
        QGE(2*q-1:2*q , :) = [g_eta1(xq,yq), g_eta2(xq,yq), g_eta3(xq,yq), g_eta4(xq,yq)];
        QGZ(2*q-1:2*q , :) = [g_zeta1(xq,yq), g_zeta1i(xq,yq), g_zeta2(xq,yq), ...
                        g_zeta2i(xq,yq), g_zeta3(xq,yq), g_zeta3i(xq,yq), ...
                        g_zeta4(xq,yq), g_zeta4i(xq,yq), g_zetam(xq,yq)];
    end


    % ====================================================================
    %% 1 - Q1/Q2 Setups
    % ====================================================================
    % Total number of rectangular elements for *both* Q1 and Q2
    num_elem = num_nodes^2;  
    num_nodes_all1 = num_nodes+1;  % Total nodes per variable in Q1
    num_nodes_all2 = 2*num_nodes+1;  % Total nodes per variable in Q2

    xy_vals2 = linspace( 0 , 1 , num_nodes_all2);  % List of all distinct values in Q2
    [X2,Y2] = meshgrid(xy_vals2 , xy_vals2);  % Set up mesh grid for Q2

    % Assign a unique global index to each unique (i,j) in the above meshes
    % Will count up from top to bottom and then from left to right
    total_nodes1 = (num_nodes_all1)^2;  % Total nodes period in Q1
    total_nodes2 = (num_nodes_all2)^2;  % Total nodes period in Q2
    glob_idx1 = reshape( linspace(1, total_nodes1, total_nodes1), ...
        [num_nodes_all1 , num_nodes_all1] );  % Q1
    glob_idx2 = reshape( linspace(1, total_nodes2, total_nodes2), ...
        [num_nodes_all2 , num_nodes_all2] );  % Q2


    % Create local connectivity arrays
    connections1 = zeros(4 , num_elem);
    connections2 = zeros(9 , num_elem);
    
    % Q1: (1,2,3,4) --> (g1,g2,g3,g4)
    elem = 0;
    for i = 1:num_nodes
        for j = 1:num_nodes
            elem = elem + 1;

            % Starting at (xi,yj), 00 > 10 > 11 > 01
            connections1(:,elem) = [glob_idx1(i , j); glob_idx1(i , j+1);
                glob_idx1(i+1 , j+1); glob_idx1(i+1 , j)];
        end
    end

    % Q2: (1,1i,2,2i,3,3i,4,4i,m) --> (g1,g2,g3,g4,g5,g6,g7,g8,g9)
    elem = 0;
    for i = 1 : 2 : (num_nodes_all2 - 2)  % Set up Q2
        for j = 1 : 2 : (num_nodes_all2 - 2)
            elem = elem + 1;

            % 00 > 10 > 20 > 21 > 22 > 12 > 02 > 01 > 11
            connections2(:,elem) = [glob_idx2(i , j); glob_idx2(i , j+1); 
                            glob_idx2(i , j+2); glob_idx2(i+1 , j+2); 
                            glob_idx2(i+2 , j+2); glob_idx2(i+2 , j+1); 
                            glob_idx2(i+2 , j); glob_idx2(i+1 , j); 
                            glob_idx2(i+1 , j+1)];
        end
    end
    



    % ====================================================================
    %% 2 - Set up constant stiffness matrices and vectors
    % ====================================================================
    % Set up our consistent BK and its results
    h = 1 / num_nodes;   BK = [h , 0; 0 , h];   det_BK = h^2;
    %======================================
    % Constant for constructing tau later
    tau_helper = h / (2 * mu);
    %======================================

    % Recurring arrays for later assembly
    D1 = sparse(total_nodes1,total_nodes1);
    D2 = sparse(total_nodes2,total_nodes2);

    % We have 7 global stiffness matrices and 2 global load vectors:
    %          A, B, Cx, and Cy; Fx and Fy
    % The local stiffness matrices for A, Cx, and Cy are consistent thanks to our uniform distribution, so we can define them as so:
    lA = (-1/45) * [-28,9,1.5,-5,1,-5,1.5,9,16;9,-88,9,16,-5,0,-5,16,48;1.5,9,-28,9,1.5,-5,1,-5,16;-5,16,9,-88,9,16,-5,0,48;1,-5,1.5,9,-28,9,1.5,-5,16;-5,0,-5,16,9,-88,9,16,48;1.5,-5,1,-5,1.5,9,-28,9,16;9,16,-5,0,-5,16,9,-88,48;16,48,16,48,16,48,16,48,-256];
    lCx = h*(1/36) * [5,1,0,0;-4,4,0,0;-1,-5,0,0; -2,-10,-10,-2;0,0,-5,-1;0,0,4,-4;0,0,1,5;10,2,2,10;-8,8,8,-8];
    lCy = h*(1/36) * [5,0,0,1;10,10,2,2;0,5,1,0;0,-4,4,0;0,-1,-5,0;-2,-2,-10,-10;-1,0,0,-5;-4,0,0,4;-8,-8,8,8];

    % Meanwhile, we'll have to construct F the long way, but the value of
    % u-bar has no bearing on its elements
    A = sparse(total_nodes2,total_nodes2);
    Cx = sparse(total_nodes2,total_nodes1);   Cy = Cx;
    Fx = zeros(total_nodes2,1);   Fy = Fx;

    for elem = 1:num_elem
        % Set up load vectors for F
        lFx = sparse(zeros(9,1));   lFy = lFx;

        % Grab the corresponding 1st local node
        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx2 == connections2(1,elem),1);
        bK = [X2(row,col) ; Y2(row,col)];

        % Set up F vector using Gaussian Quadratures
        for q = 1:9  % Quadrature index for weight/node
            Xg = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];

            lFx = lFx + Qweights(q) * f_x(Xg(1) , Xg(2)) * QZ(q,:)';
            lFy = lFy + Qweights(q) * f_y(Xg(1) , Xg(2)) * QZ(q,:)';
        end
        lFx = det_BK * lFx;   lFy = det_BK * lFy;


        % Obtain Q1 and Q2 global indices
        g_ind1 = connections1(: , elem);
        g_ind2 = connections2(: , elem);

        % Update values of global stiffness matrices and load vectors
        A(g_ind2,g_ind2) = A(g_ind2,g_ind2) + lA;
        Cx(g_ind2,g_ind1) = Cx(g_ind2,g_ind1) + lCx;
        Cy(g_ind2,g_ind1) = Cy(g_ind2,g_ind1) + lCy;
        Fx(g_ind2) = Fx(g_ind2) + lFx;
        Fy(g_ind2) = Fy(g_ind2) + lFy;
    end


    % Create 1D array of all Q2 global indices on the boundary
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
    %% 3 - Run the Loop Iterations
    % ====================================================================
    % Create u_bar_x and u_bar_y arrays for our iteration
    u_bar_x = zeros(total_nodes2,1);   u_bar_y = u_bar_x;
    for i = 1:num_nodes_all2
        for j = 1:num_nodes_all2
            u_bar_x((i-1)*num_nodes_all2+j) = u0_x(xy_vals2(i),xy_vals2(j));
            u_bar_y((i-1)*num_nodes_all2+j) = u0_y(xy_vals2(i),xy_vals2(j));
        end
    end
    % Set up u_(k-1) and u_(k-2)
    u_k1_x = u_bar_x;   u_k1_y = u_bar_y;
    u_k2_x = u_bar_x;   u_k2_y = u_bar_y;


    % Tolerance for ending loop early, max iterations, and storage constant
    ToL = 1e-6;   maxIter = 100;   rho = 1;

    for iter = 1:maxIter
        B = sparse(total_nodes2,total_nodes2);
        %============================================
        % New arrays for stabilization
        %============================================
        As = B;   Bs = B;
        Cxs = sparse(total_nodes2,total_nodes1);   Cys = Cxs;
        Fxs = zeros(total_nodes2,1);   Fys = Fxs;
        %============================================
        
        for elem = 1:num_elem
            % Set up local stiffness matrices for B
            lB = sparse(9,9); 

            % Explicitly set up functions for u_bar_x and u_bar_y based on
            % the current element K, which given global nodes n1,n1i, is:
            %   u(n1) * zeta_1 + u(n1i) * zeta_1i + u(n2) * zeta_2 +...
            UB_x = QZ * u_bar_x(connections2(:,elem));
            UB_y = QZ * u_bar_y(connections2(:,elem));
            UB = [UB_x , UB_y];

            %=============================================
            % New local arrays for stabilization, plus tau
            %=============================================
            lAs = lB;   lBs = lB;
            lCxs = sparse(9,4);   lCys = lCxs;
            lFxs = zeros(9,1);   lFys = lFxs;

            % Value of u-bar at the center of the current element
            u_mid = [u_bar_x(connections2(9,elem)) ; u_bar_y(connections2(9,elem))];
            perch = tau_helper * u_mid;
            tau = mu * (coth(perch(1)) * perch(1) + coth(perch(2)) * perch(2) - 2) / (norm(u_mid)^2);

            % Compute all dot product results of u-bar * gradient zeta_i
            dot_prods = sparse(9,9);
            for k = 1:9
                dot_prods(k,:) = UB(k,:) * QGZ(2*k-1 : 2*k , :);
            end

            [row,col] = find(glob_idx2 == connections2(1,elem),1);
            bK = [X2(row,col) ; Y2(row,col)];
            %=============================================

            % Set up our matrices/arrays using Gaussian Quadratures
            for q = 1:9  % Quadrature index for weight/node
                lB = lB + Qweights(q) * QZ(q,:)' * UB(q,:) * QGZ(2*q-1:2*q , :);
                %====================================================
                % Set up our local stiffness matrices
                %====================================================
                lAs = lAs + Qweights(q) * dot_prods(q,:)' * QLZ(q,:);
                lBs = lBs + Qweights(q) * dot_prods(q,:)' * dot_prods(q,:);
                lCxs = lCxs + Qweights(q) * dot_prods(q,:)' * QGE(2*q-1,:);
                lCys = lCys + Qweights(q) * dot_prods(q,:)' * QGE(2*q,:);

                Xg = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];
                lFxs = lFxs + Qweights(q) * dot_prods(q,:)' * f_x(Xg(1) , Xg(2));
                lFys = lFys + Qweights(q) * dot_prods(q,:)' * f_y(Xg(1) , Xg(2));
                %====================================================
            end
            lB = h * lB;
            %====================================================
            lAs = (1/h) * tau * lAs;
            lBs = tau * lBs;
            lCxs = tau * lCxs;
            lCys = tau * lCys;
            lFxs = h * tau * lFxs;
            lFys = h * tau * lFys;
            %====================================================


            % Obtain Q1 and Q2 global indices
            g_ind1 = connections1(: , elem);
            g_ind2 = connections2(: , elem);

            % Update values of the global stiffness matrices for B
            B(g_ind2,g_ind2) = B(g_ind2,g_ind2) + lB;
            %====================================================
            As(g_ind2,g_ind2) = As(g_ind2,g_ind2) + lAs;   
            Bs(g_ind2,g_ind2) = Bs(g_ind2,g_ind2) + lBs;
            Cxs(g_ind2,g_ind1) = Cxs(g_ind2,g_ind1) + lCxs;   
            Cys(g_ind2,g_ind1) = Cys(g_ind2,g_ind1) + lCys;
            Fxs(g_ind2) = Fxs(g_ind2) + lFxs;   
            Fys(g_ind2) = Fys(g_ind2) + lFys;
            %====================================================
        end


        %=====================================================
        % Combined stiffness matrix K will be of the form:
        K = [mu*(A-As)+B+Bs ,  D2, Cx+Cxs;
            D2 , mu*(A-As)+B+Bs , Cy+Cys;
            Cx' , Cy' , D1];
        % Meanwhile the combined load vector g will be:
        g = [Fx + Fxs; Fy + Fys; zeros(total_nodes1, 1)];
        %=====================================================
    
    
        % Set all rows and columns corresponding to bound indices to 0
        K(bound_idx , :) = 0;   K(: , bound_idx) = 0;  g(bound_idx) = 0;
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
        %z0 = [u_bar_x ; u_bar_y ; zeros(total_nodes1,1)];
        % z0 = zeros(2*total_nodes2 + total_nodes1,1);
        % M = [K(1:2*total_nodes2 , 1:2*total_nodes2) , zeros(2*total_nodes2,total_nodes1);
        %     zeros(total_nodes1,2*total_nodes2) , Q];
        %z = gmres(K, g, 20, 1e-6, 100, M);
        %z = new_gmres(K, z0, g, M, 20, 100, 1e-6 );
        uh_x = z(1:total_nodes2);  uh_y = z(total_nodes2+1:2*total_nodes2);
        ph = z(2*total_nodes2+1:length(z));


        % End the loop early if norms of u_(k-1)-u_k is < ToL
        if norm(uh_x - u_k1_x) < ToL && norm(uh_y - u_k1_y) < ToL
            fprintf('Converged at iteration %d\n', iter);  break;
        end


        % Otherwise, update u(k-2), u(k-1), and u-bar 
        % u(k-2) = u(k-1)
        u_k2_x = u_k1_x;   u_k2_y = u_k1_y;
        % u(k-1) = u_k
        u_k1_x = uh_x;   u_k1_y = uh_y;
        % u-bar = p u(k-1) + (1-p) u(k-2)
        u_bar_x = rho * u_k1_x + (1-rho) * u_k2_x;  
        u_bar_y = rho * u_k1_y + (1-rho) * u_k2_y;
        %disp(norm([u_k1_x ; u_k1_y] - [u_k2_x ; u_k2_y]))
    end




    % ====================================================================
    %% 4 - Compute L2 error estimates by element
    % ====================================================================
    % This'll be the result of summing the L2 errors over the elements,
    % so we first set the total accumulated errors to 0
    ux_l2 = 0;   uy_l2 = 0;   p_l2 = 0;
    % Small optimization since going forward, we need det_BK * Qweights(q)
    Qweights = det_BK * Qweights;


    % Now we loop over each element to obtain individual L2 errors
    for elem = 1:num_elem
        % On the reference element:
        %   uh = uh(n1) * zeta1 + uh(n1i) * zeta1i + ...
        UH_x = QZ * uh_x(connections2(:,elem));
        UH_y = QZ * uh_y(connections2(:,elem));

        % Grab the corresponding 1st local node (same for both Q2 and Q1)
        % This is needed for X = bK + BK X_hat
        [row,col] = find(glob_idx2 == connections2(1,elem),1);
        bK = [X2(row,col) ; Y2(row,col)];

        % Ditto, but for the Q1 nodes (relevant to p)
        PH = QE * ph(connections1(:,elem));

        for q = 1:9  % Quadrature index for weight/node
            Xg = bK + BK * [Qx_nodes(q) ; Qy_nodes(q)];  % Global nodes

            % Update L2 norms of ux, uy, and p
            ux_l2 = ux_l2 + Qweights(q) * (UH_x(q) - u_x(Xg(1),Xg(2)))^2;
            uy_l2 = uy_l2 + Qweights(q) * (UH_y(q) - u_y(Xg(1),Xg(2)))^2;
            p_l2 = p_l2 + Qweights(q) * (PH(q) - p(Xg(1),Xg(2)))^2;
        end
    end


    % Update our errors since they're the square roots of the above sums
    l2s = sqrt([ux_l2 , uy_l2 , p_l2]);

    % Print the results
    fprintf('ux L2 Norm = %12.6f\n', l2s(1));
    fprintf('uy L2 Norm = %12.6f\n', l2s(2));
    fprintf('p L2 Norm = %12.6f\n', l2s(3));
    toc  % Output runtime
    fprintf('\n') 
end