/*
	This function solves a system of linear equations, given the augmented matrix
	in 'a', the solutions to the respected equations in 'y_vals', the square dimention
	'n' of the augmented matrix, and the shared input 'type' of 'a' and 'y_vals'. It
	requires that 'type' be either 'i' or 'f' for integers or floats respectively,
	returning NULL otherwise.
 */
float *solveLinearSystem(void *a, void *y_vals, int n, char type){
    int size = sizeof(float);
    int matrix_size = n*n*size;
    float *coef = (float*)malloc(matrix_size);
    float *coef_temp = (float*)malloc(matrix_size);
    float *e = (float*)malloc(matrix_size);
    float *y = (float*)malloc(n*size);
    float *y_temp = (float*)malloc(n*size);
    float *x = (float*)malloc(n*size);
    float *identity = (float*)malloc(n*n*sizeof(float));
    
    for(int i = 0; i<n; i++){
        for(int j = 0; j<n; j++){
            if(i==j) *(identity + i*n + j) = 1;
            else *(identity + i*n + j) = 0;
        }
    }
    
    //Copy matrix 'a' into the new augmented matrix (stretch integer input into float matrix if needed)
    for(int i = 0; i < n; i++){
        for(int j = 0; j <n; j++){
            if(type == 'f') *(float*)((char*)coef + sizeof(float)*(i*n+j)) = *(float*)((char*)a + size*(i*n+j));
            else if (type == 'i') *(float*)((char*)coef + sizeof(float)*(i*n+j)) = *(int*)((char*)a + size*(i*n+j));
            else return NULL;
        }
    }
    
    //Copy matrix 'y_vals' into the new y matrix (stretch integer input into float matrix if needed)
    for(int i = 0; i < n; i++){
        if(type == 'f') *(float*)((char*)y + sizeof(float)*i) = *(float*)((char*)y_vals + size*i);
        else if (type == 'i') *(float*)((char*)y + sizeof(float)*i) = *(int*)((char*)y_vals + size*i);
        else return NULL;
    }
    
    for(int row = 0; row < n; row++){
        //Copy identity into e matrix
        memcpy(e, identity, matrix_size);
        *(float*)((char*)e + row*n*size + row*size) = (1/(*(float*)((char*)coef + row*n*size + row*size)));
        
        //Multiply by augmented matrix
        float sum;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                sum = 0;
                for(int k = 0; k < n; k++){
                    sum += *(float*)((char*)e + i*n*size + k*size) * *(float*)((char*)coef + k*n*size + j*size);
                }
                
                *(float*)((char*)coef_temp + i*n*size + j*size) = sum;
            }
        }
        
        //Copy temp into coef
        memcpy(coef, coef_temp, matrix_size);
        
        //multiply by y
        for(int i = 0; i < n; i++){//each row of e mat
            float sum = 0;
            for(int j = 0; j < n; j++){ //index of e column/yrow
                float first = *(float*)((char*)e + i*n*size + j*size);
                float second = *(float*)((char*)y+j*size);
                sum += first * second;
            }
            *(float*)((char*)y_temp + i*size) = sum;
        }
        memcpy(y, y_temp, n*size);
        
        //Zero out positions below diagonal
        if(row + 1< n){
            for(int j = row+1; j < n; j++){
                //Copy identity into e matrix
                memcpy(e, identity, matrix_size);
                //Change e matrix to negative sign
                *(float*)((char*)e + j*n*size + row*size) = -*(float*)((char*)coef + j*n*size + row*size);
                //Multiply by augmented matrix
                float sum;
                for(int i = 0; i < n; i++){
                    for(int j = 0; j < n; j++){
                        sum = 0;
                        for(int k = 0; k < n; k++){
                            sum += *(float*)((char*)e + i*n*size + k*size) * *(float*)((char*)coef + k*n*size + j*size);
                        }
                        
                        *(float*)((char*)coef_temp + i*n*size + j*size) = sum;
                    }
                }
                
                //Copy temp into augmented matrix
                memcpy(coef, coef_temp, matrix_size);
                
                //Multiply by y
                for(int i = 0; i < n; i++){//each row of e mat
                    float sum = 0;
                    for(int j = 0; j < n; j++){ //index of e column/yrow
                        float first = *(float*)((char*)e + i*n*size + j*size);
                        float second = *(float*)((char*)y+j*size);
                        sum += first * second;
                    }
                    *(float*)((char*)y_temp + i*size) = sum;
                }
                //Copy y_temp into y
                memcpy(y, y_temp, n*size);
            }
        }
    }
    
    printf("Final echelon form:\n");
    for(int i = 0; i < n*n; i++){
        printf("%1.3f ", *(float*)((char*)coef + i*size));
        if(i%n == n-1) printf("\n");
    }
    printf("\nFinal Y:\n");
    for(int i = 0; i < n; i++){
        printf("%1.3f\n", *(float*)((char*)y+i*size));
    }
    //Gaussian  elimination
    for(int i = n-1; i >= 0; i--){
        float temp = 0;
        for(int j = i+1; j<n; j++){
            temp += *(float*)((char*)coef + i*n*size +j*size) * *(float*)((char*)x+j*size);
        }
        *(float*)((char*)x+i*size) = *(float*)((char*)y+i*size) - temp;
    }
    
    //Linear System Solved----------------------------------------------------------------------------------
    
    printf("\nReturning X:\n");
    for(int i = 0; i < n; i++){
        printf("%1.3f\n", *(float*)((char*)x+i*size));
    }
    
    return x;
}
