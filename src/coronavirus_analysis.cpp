// Justin Santer
// AMS 326 HW 2
// Problem 2.2
// March 27, 2020

#include <random>
#include <time.h>

const int INIT_PATIENTS = 10;
const int BLOCK_OF_DAYS = 45;
const double GROWTH_MEAN = 0.18;
const double GROWTH_STDDEV = 0.08;
const double DECAY_MEAN = -0.24;
const double DECAY_STDDEV = 0.04;
const int NUM_SELECTED_POINTS = 5;

void generate_values(double *values, int size){

    std::default_random_engine generator1;
    generator1.seed(time(0));
    std::normal_distribution<double> distr1(GROWTH_MEAN, GROWTH_STDDEV);

    double percent_change;

    for (int i = 2; i < size/2 + 1; i++){
        percent_change = distr1(generator1);
        values[i] = (percent_change + 1.0) * values[i-1];
    }

    std::default_random_engine generator2;
    generator2.seed(time(0));
    std::normal_distribution<double> distr2(DECAY_MEAN, DECAY_STDDEV);

    for (int i = size/2 + 1; i < size; i++){
        percent_change = distr2(generator2);
        values[i] = (percent_change + 1.0) * values[i-1];
    }

    return;
}

void print_values(double *values, int size){
    for (int i = 1; i < size; i++){
        printf("DAY %2d: %.0f\n", i, values[i]);
    }

    return;
}

void print_system_matrix(double **matrix){
    for(int i = 0; i < NUM_SELECTED_POINTS; i++){
        for(int j = 0; j < NUM_SELECTED_POINTS+1; j++){
            if (j != NUM_SELECTED_POINTS){
                printf("%14.0f ", matrix[i][j]);
            }
            else{
                printf("%14f ", matrix[i][j]);
            }
        }
        printf("\n");
    }
}

double **construct_system_matrix(double *selected_xvalues, double *selected_yvalues){
    double **A = new double*[NUM_SELECTED_POINTS];
    for (int i = 0; i < NUM_SELECTED_POINTS; i++){
        A[i] = new double[NUM_SELECTED_POINTS+1];
    }

    for (int i = 0; i < NUM_SELECTED_POINTS; i++){
        int n = NUM_SELECTED_POINTS;
        for (int j = 0; j < NUM_SELECTED_POINTS+1; j++){
            if (j != NUM_SELECTED_POINTS){
                A[i][j] = pow(selected_xvalues[i], n-1-j);
            }
            else{
                A[i][j] = selected_yvalues[i];
            }
        }
    }

    return A;
}

void delete_system_matrix(double **A){
    for (int i = 0; i < NUM_SELECTED_POINTS; i++){
        delete A[i];
    }
    delete A;
}

double *gaussian_elimination(double **matrix){
    double *a = new double[NUM_SELECTED_POINTS];

    double b;

    for(int j = 0; j < NUM_SELECTED_POINTS; j++){
        for(int i = 0; i < NUM_SELECTED_POINTS; i++){
            if(i != j){
                b = matrix[i][j] / matrix[j][j];
                for(int k = 0; k < NUM_SELECTED_POINTS+1; k++){
                    matrix[i][k] = matrix[i][k] - b * matrix[j][k];
                }
            }
        }
    }

    for(int i = 0; i < NUM_SELECTED_POINTS; i++){
        a[i] = matrix[i][NUM_SELECTED_POINTS] / matrix[i][i];
        matrix[i][NUM_SELECTED_POINTS] = a[i];
        matrix[i][i] = 1.0;
    }

    return a;
}

double mean_of_xvalues(int starting_day, int size){
    double sum = 0;

    int count = 0;
    for (int i = starting_day; i < size; i++){
        count++;
        sum += (double) i;
    }

    return (sum / (double) count);
}

double mean_of_lnyvalues(double *values, int starting_day, int size){
    double sum = 0;

    int count = 0;
    for (int i = starting_day; i < size; i++){
        count++;
        sum += log(values[i]);
    }

    return (sum / (double) count);
}

void fit_line(double *values, int starting_day, int size, double *a, double *b){
    double sum_of_x = 0.0;
    double sum_of_x2 = 0.0;
    double sum_of_xy = 0.0;

    double x_mean = mean_of_xvalues(starting_day, size);
    double y_mean = mean_of_lnyvalues(values, starting_day, size);

    int n = 0;
    for (int i = starting_day; i < size; i++){
        sum_of_x += (double)i;
        sum_of_x2 += pow((double)i, 2);
        sum_of_xy += ((double) i * log(values[i]));
        n++;
    }

    double lna;
    *b = (sum_of_xy - n * x_mean * y_mean) / (sum_of_x2 - n * pow(x_mean, 2));
    lna = (y_mean * sum_of_x2 - x_mean * sum_of_xy) / (sum_of_x2 - n * pow(x_mean, 2));
    *a = exp(lna);
}

int main(){

    double function_values[2*BLOCK_OF_DAYS+1];

    function_values[0] = 0.0;
    function_values[1] = INIT_PATIENTS;

    generate_values(function_values, 2*BLOCK_OF_DAYS+1);

    print_values(function_values, 2*BLOCK_OF_DAYS+1);

    double step = (double) BLOCK_OF_DAYS / NUM_SELECTED_POINTS;
    double selected_xvalues[NUM_SELECTED_POINTS];
    for (int i = 0; i < NUM_SELECTED_POINTS; i++){
        selected_xvalues[i] = step + i * step;
    }
    double selected_yvalues[NUM_SELECTED_POINTS];
    for (int i = 0; i < NUM_SELECTED_POINTS; i++){
        selected_yvalues[i] = function_values[(int)selected_xvalues[i]];
    }

    double **matrix = construct_system_matrix(selected_xvalues, selected_yvalues);

    double *coefficents = gaussian_elimination(matrix);


    printf("\nFirst Half Polynomial Interpolation:\n");
    printf("f(x) = %f", coefficents[NUM_SELECTED_POINTS-1]);
    int j = 1;
    for(int i = NUM_SELECTED_POINTS-2; i >= 0; i--){
        printf(" + (%f)x^%d", coefficents[i], j);
        j++;
    }
    printf("\n\n");
    delete coefficents;

    double a, b;
    fit_line(function_values, BLOCK_OF_DAYS, 2*BLOCK_OF_DAYS+1, &a, &b);

    printf("Second Half Exponential Line Fitting:\n");
    printf("f(x) = %fexp(%fx)\n", a, b);
    delete_system_matrix(matrix);
    return 0;
}