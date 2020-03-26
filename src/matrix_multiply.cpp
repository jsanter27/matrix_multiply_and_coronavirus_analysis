// Justin Santer
// AMS 326 HW 2
// Problem 2.1
// March 20, 2020

#include <random>
#include <time.h>


const int SIZE_1 = pow(2, 10);
const int SIZE_2 = pow(2, 12);
const double MIN_RANGE = -1.0;
const double MAX_RANGE = 1.0;

unsigned long long opcount = 0;


double rand_double(double a, double b){
    double random = ((double) rand()) / (double) RAND_MAX;
    double difference = b - a;
    return random * difference + a;
}


double **generate_matrix(int size, bool fill){
    // Declare and define new matrix
    double **A = new double*[size];
    for (int i = 0; i < size; i++){
        A[i] = new double[size];
    }

    if (fill){
        // Loop through elements to randomly insert elements
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                A[i][j] = rand_double(MIN_RANGE, MAX_RANGE);
            }
        }
    }

    return A;   
}


void delete_matrix(double **A, int size){
    for (int i = 0; i < size; i++){
        delete A[i];
    }
    delete A;
}


void print_matrix(double **A, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            printf("%8.4f ", A[i][j]);
        }
        printf("\n");
    }
}


bool equal_matrices(double **A, double **B, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            if (A[i][j] != B[i][j]){
                return false;
            }
        }
    }
    return true;
}


void add_matrices(double **A, double **B, double **C, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return;
}


void subtract_matrices(double **A, double **B, double**C, int size){
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }

    return;
}


void naive_helper(double **A, double **B, double **C, int size){
    double element;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            element = 0;
            for (int k = 0; k < size; k++){
                element += A[i][k] * B[k][j];
                opcount += 2;
            }
            C[i][j] = element;
        }
    }

    return;
}


double **naive_method(double **A, double **B, int size){
    double **C = generate_matrix(size, false);
    
    naive_helper(A, B, C, size);

    return C;
}


void strassen_helper(double *A[], double *B[], double *C[], int size, int levels){
    if (size == 2 || levels == 0){
        naive_helper(A, B, C, size);
        
        return;
    }

    int next_size = size / 2;

    //GET THE SUB MATRICES FOR A

    double **A_11 = generate_matrix(next_size, false);
    double **A_12 = generate_matrix(next_size, false);
    double **A_21 = generate_matrix(next_size, false);
    double **A_22 = generate_matrix(next_size, false);

    for (int i = 0; i < next_size; i++){
        for (int j = 0; j < next_size; j++){
            A_11[i][j] = A[i][j];
            A_12[i][j] = A[i][next_size+j];
            A_21[i][j] = A[next_size+i][j];
            A_22[i][j] = A[next_size+i][next_size+j];
        }
    }

    //GET THE SUB MATRICES FOR B
    double **B_11 = generate_matrix(next_size, false);
    double **B_12 = generate_matrix(next_size, false);
    double **B_21 = generate_matrix(next_size, false);
    double **B_22 = generate_matrix(next_size, false);

    for (int i = 0; i < next_size; i++){
        for (int j = 0; j < next_size; j++){
            B_11[i][j] = B[i][j];
            B_12[i][j] = B[i][next_size+j];
            B_21[i][j] = B[next_size+i][j];
            B_22[i][j] = B[next_size+i][next_size+j];
        }
    }

    //GENERATE M MATRICES
    double **M_1 = generate_matrix(next_size, false);
    double **M_2 = generate_matrix(next_size, false);
    double **M_3 = generate_matrix(next_size, false);
    double **M_4 = generate_matrix(next_size, false);
    double **M_5 = generate_matrix(next_size, false);
    double **M_6 = generate_matrix(next_size, false);
    double **M_7 = generate_matrix(next_size, false);

    double **temp1 = generate_matrix(next_size, false);
    double **temp2 = generate_matrix(next_size, false);


    //SEVEN SUB MULTIPLICATIONS
    add_matrices(A_11, A_22, temp1, next_size);
    add_matrices(B_11, B_22, temp2, next_size);
    strassen_helper(temp1, temp2, M_1, next_size, levels-1);

    add_matrices(A_21, A_22, temp1, next_size);
    strassen_helper(temp1, B_11, M_2, next_size, levels-1);

    subtract_matrices(B_12, B_22, temp1, next_size);
    strassen_helper(A_11, temp1, M_3, next_size, levels-1);

    subtract_matrices(B_21, B_11, temp1, next_size);
    strassen_helper(A_22, temp1, M_4, next_size, levels-1);

    add_matrices(A_11, A_12, temp1, next_size);
    strassen_helper(temp1, B_22, M_5, next_size, levels-1);

    subtract_matrices(A_21, A_11, temp1, next_size);
    add_matrices(B_11, B_12, temp2, next_size);
    strassen_helper(temp1, temp2, M_6, next_size, levels-1);

    subtract_matrices(A_12, A_22, temp1, next_size);
    add_matrices(B_21, B_22, temp2, next_size);
    strassen_helper(temp1, temp2, M_7, next_size, levels-1);

    delete_matrix(A_11, next_size);
    delete_matrix(A_12, next_size);
    delete_matrix(A_21, next_size);
    delete_matrix(A_22, next_size);

    // CALCULATE SUB C MATRICES AND COMBINE THEM TO GREATER C
    double **C_11 = B_11;
    add_matrices(M_1, M_4, C_11, next_size);
    subtract_matrices(C_11, M_5, C_11, next_size);
    add_matrices(C_11, M_7, C_11, next_size);

    double **C_12 = B_12;
    add_matrices(M_3, M_5, C_12, next_size);

    double **C_21 = B_21;
    add_matrices(M_2, M_4, C_21, next_size);

    double **C_22 = B_22;
    subtract_matrices(M_1, M_2, C_22, next_size);
    add_matrices(C_22, M_3, C_22, next_size);
    add_matrices(C_22, M_6, C_22, next_size);

    for (int i = 0; i < next_size; i++){
        for (int j = 0; j < next_size; j++){
            C[i][j] = C_11[i][j];
            C[i][j+next_size] = C_12[i][j];
            C[i+next_size][j] = C_21[i][j];
            C[i+next_size][j+next_size] = C_22[i][j];
        }
    }

    delete_matrix(M_1, next_size);
    delete_matrix(M_2, next_size);
    delete_matrix(M_3, next_size);
    delete_matrix(M_4, next_size);
    delete_matrix(M_5, next_size);
    delete_matrix(M_6, next_size);
    delete_matrix(M_7, next_size);
    delete_matrix(temp1, next_size);
    delete_matrix(temp2, next_size);
    delete_matrix(B_11, next_size);
    delete_matrix(B_12, next_size);
    delete_matrix(B_21, next_size);
    delete_matrix(B_22, next_size);
}


double **strassen_method(double **A, double **B, int size, int levels){
    double **C = generate_matrix(size, false);

    strassen_helper(A, B, C, size, levels);

    return C;
}


int main(){
    // Set random seed
    srand(time(0));

    for (int i = 0; i < 2; i++){
        int size;
        if (i == 0){
            size = SIZE_1;
        }
        else{
            size = SIZE_2;
        }
        // Generate two random matrices
        printf("Generating first random %d x %d matrix...\n", size, size);
        double **A = generate_matrix(size, true);
        // print_matrix(A, SIZE_1);
        printf("Generating second random %d x %d matrix...\n", size, size);
        double **B = generate_matrix(size, true);
        // print_matrix(B, SIZE_1);

        double start_time, end_time;

        // Multiply the two matrices using NAIVE METHOD and print results
        printf("\nMultiplying the two matrices using Naive Multiplication...\n");
        start_time = time(0);
        double **naive = naive_method(A, B, size);
        end_time = time(0);
        double naive_time = (end_time - start_time) / 60.0;

        // Print the amount of operations that NAIVE METHOD used
        // print_matrix(naive, SIZE_1);
        printf("Number of operations when using Naive Multiplication: %llu\n", opcount);
        printf("Time taken using Naive Multiplication: %.2f minutes\n\n", naive_time);
        opcount = 0;

        // Multiply using STRASSEN METHOD
        printf("Multiplying the two matrices using Strassen's algorithm for three levels...\n");
        start_time = time(0);
        double **strassen = strassen_method(A, B, size, 3); 
        end_time = time(0);
        double strassen_time = (end_time - start_time) / 60.0;

        // Print results
        // print_matrix(strassen, SIZE_1);
        printf("Number of operations when using Strassen's algorithm: %llu\n", opcount);
        printf("Time taken using Strassen's algorithm: %.2f minutes\n\n", strassen_time);

        printf(equal_matrices(naive, strassen, size) ? "The matrices are NOT EQUAL\n\n" : "The matrices are EQUAL\n\n");

        // Delete generated matrices
        delete_matrix(A, size);
        delete_matrix(B, size);
        delete_matrix(naive, size);
        delete_matrix(strassen, size);
    }

    return 0;
}
