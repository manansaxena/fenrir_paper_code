functions {  
  void print_vector(vector v) {
        for (i in 1:size(v)) {
            print("v[", i, "] = ", v[i]);
        }
    }
  void print_array(array[] int v) {
        for (i in 1:size(v)) {
            print("v[", i, "] = ", v[i]);
        }
    }
  void print_matrix(matrix m) {
      for (i in 1:rows(m)) {
          for (j in 1:cols(m)) {
              print("m[", i, ", ", j, "] = ", m[i, j]);
          }
      }
    }
  matrix make_symmetric(matrix X){
        return (X+X')/2;
    }

  real gmdlm_lpdf(matrix eta, array[]vector F, array[] matrix G, array[] matrix W,
                array[] real gamma, array[] matrix M0, array[] matrix C0, matrix Xi0, real upsilon0, 
                array[] int observed_TT, array[] int N_total_list, int num_timeseries){
    
        int TT = size(observed_TT); // total number of timepoints (includes missing values + all timeseries)
        int p = rows(eta);
        int q = rows(G[1]);
        
        // system matricies for each iteration
        vector[q] F_t;
        matrix[q,q] G_t;
        matrix[q,q] W_t;
        matrix[p,p] Xi;
        real upsilon;
        real gamma_t;
        
        // internals to filter
        matrix[q,p] M;
        matrix[q,q] C;
        matrix[q,p] A; 
        matrix[q,q] R;
        vector[q] S;
        vector[p] f;
        real qq;
        vector[p] e;
        
        real log_prob_accum = 0;
        int pos = 1;
        int t = 1;

        Xi = Xi0;
        upsilon = upsilon0;

        for(timeseries in 1:num_timeseries){
        M = M0[timeseries];
        C = C0[timeseries];
        for(i in 1:N_total_list[timeseries]){
            F_t = F[t];
            G_t = G[t];
            W_t = W[t];
            gamma_t = gamma[t];

            A = G_t*M;
            R = quad_form(C, G_t') + W_t;
            
            if (observed_TT[t] == 1){
                f = A'*F_t;
                qq = quad_form(R, F_t) + gamma_t;
                log_prob_accum += multi_student_t_lpdf(eta[,pos]  | upsilon, f, (qq*Xi)/upsilon);

                S = R*F_t/qq;
                e = eta[,pos] - f;
                M = A + S*e';
                C = R - tcrossprod(to_matrix(S))*qq;
                upsilon = upsilon + 1;
                Xi = Xi + (e*e')/qq;
                pos = pos + 1;
            }
            else {
                M = A;
                C = R;
                upsilon = upsilon;
                Xi = Xi;
            }
            t = t + 1;
            } 
        }
        return(log_prob_accum);
    }

  matrix matrix_normal_rng(matrix C, matrix M, matrix Sigma){
        int p = cols(M);
        int q = rows(M);
        matrix[q, p] X;
    
        for (i in 1:q){
            for (j in 1:p){
                X[i,j] = normal_rng(0,1);
            }
        }
        matrix[q,p] theta = M + cholesky_decompose(C)*X*cholesky_decompose(Sigma)';
        return theta;
    }
  array[] matrix gmdlm_smoothing_rng(matrix eta, array[] vector F, array[] matrix G, array[] matrix W,
                                                array[] real gamma,array[] matrix M0, array[] matrix C0, matrix Xi0, real upsilon0,
                                                array[] int observed_TT, array[] int N_total_list, int num_timeseries){
        int TT = size(observed_TT);
        int p = rows(eta);
        int q = rows(G[1]);

        vector[q] F_t;
        matrix[q,q] G_t;
        matrix[q,q] W_t;
        array[TT+1] matrix[p,p] Xi;
        array[TT+1] real upsilon;
        real gamma_t;

        array[TT+num_timeseries] matrix[q,p] M;
        array[TT+num_timeseries] matrix[q,q] C;
        array[TT] matrix[q,p] A; 
        array[TT] matrix[q,q] R;
        vector[q] S;
        vector[p] f;
        real qq;
        vector[p] e;
        matrix[q,q] Z;

        array[TT] matrix[q,p] theta;
        int pos = 1;
        int t = 1;

        Xi[1] = Xi0;
        upsilon[1] = upsilon0;
        
        for(timeseries in 1:num_timeseries){
        M[t + timeseries - 1] = M0[timeseries];
        C[t + timeseries - 1] = C0[timeseries];
        for (i in 1:N_total_list[timeseries]){
            F_t = F[t];
            G_t = G[t];
            W_t = W[t];
            gamma_t = gamma[t];

            A[t] = G_t*M[t + timeseries - 1];
            R[t] = quad_form(C[t + timeseries - 1], G_t') + W_t;

            if(observed_TT[t] == 1){
                f = A[t]'*F_t;
                qq = quad_form(R[t], F_t) + gamma_t;
                
                S = R[t]*F_t/qq;
                e = eta[,pos] - f;
                M[t+timeseries] = A[t] + S*e';
                C[t+timeseries] = R[t]- tcrossprod(to_matrix(S))*qq;
                upsilon[t+1] = upsilon[t] + 1;
                Xi[t+1] = Xi[t] + (e*e')/qq;
                pos = pos + 1;
            }
            else {
                M[t+timeseries] = A[t];
                C[t+timeseries] = R[t];
                upsilon[t+1] = upsilon[t];
                Xi[t+1] = Xi[t];
            }
            t = t + 1;
        }
        }

        Xi[TT+1] = make_symmetric(Xi[TT+1]);
        matrix[p,p] Sigma = inv_wishart_rng(upsilon[TT+1], Xi[TT+1]);

        int k = 0;
        int reset_flag = 1;

        array[num_timeseries] int array_from_1_to_N;
        for (i in 1:num_timeseries) {
            array_from_1_to_N[i] = i;
        }
        array_from_1_to_N = reverse(array_from_1_to_N);

        for(timeseries in array_from_1_to_N){
            reset_flag = 1;
            for (i in 1:N_total_list[timeseries]){
                if(reset_flag == 1){
                    C[TT - k + timeseries] = make_symmetric(C[TT - k + timeseries]);
                    theta[TT - k] = matrix_normal_rng(C[TT - k + timeseries], M[TT - k + timeseries], Sigma);
                    reset_flag = 0;
                }else {
                    G_t = G[TT - k + 1];
                    Z = C[TT + timeseries - k] * G_t' * inverse(R[TT - k + 1]);
                    C[TT + timeseries - k] = C[TT + timeseries - k] - quad_form(R[TT - k + 1], Z');
                    C[TT + timeseries - k] = make_symmetric(C[TT + timeseries - k]);
                    M[TT + timeseries - k] = M[TT + timeseries - k] + Z * (theta[TT - k + 1] - A[TT - k + 1]);
                    theta[TT - k] = matrix_normal_rng(C[TT + timeseries - k], M[TT + timeseries - k], Sigma);
                }
                k = k + 1;
            }
        }
        return(theta);
    }

  vector softmax_id(vector alpha) {
        vector[num_elements(alpha) + 1] alphac1;
        for (k in 1:num_elements(alpha))
            alphac1[k] = alpha[k];
            alphac1[num_elements(alphac1)] = 0;
        return softmax(alphac1);
    }
}