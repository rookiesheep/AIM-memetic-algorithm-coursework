//6518568 Ziyang Cheng

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

/* parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 5;
static int POP_SIZE = 100; //global parameters
int MAX_NUM_OF_GEN = 1000; //max number of generations
int MAX_TIME = 60;  //max amount of time permited (in sec)
float CROSSOVER_RATE = 0.8;
float MUTATION_RATE = 0.1;

struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1
float rand_01() {
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max) {
    int div = max-min+1;
    int val =rand() % div + min;
    // printf("rand_range= %d \n", val);
    return val;
}

struct item_struct {
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
};

struct problem_struct {
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
};

void free_problem(struct problem_struct* prob) {
    if(prob!=NULL) {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL) {
            for(int j=0; j<prob->n; j++) {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

void init_problem(int n, int dim, struct problem_struct** my_prob, char filename) {
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}

//memetic algorithm
int MA(struct problem_struct* prob) {
    struct solution_struct curt_pop[POP_SIZE];
    struct solution_struct new_pop[POP_SIZE];
    init_population(prob, curt_pop);
    init_population(prob, new_pop);
    int gen=0;
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    while(gen<MAX_NUM_OF_GEN && time_spent < MAX_TIME) {
        cross_over(curt_pop, new_pop);
        mutation(new_pop);
        feasibility_repair(new_pop);
        local_search_first_descent(new_pop);
        replacement(curt_pop, new_pop);
        gen++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }
    
    update_best_solution(curt_pop);
    
    
    free_population(curt_pop, POP_SIZE);
    free_population(new_pop, POP_SIZE);
    
    return 0;
}

int main(int argc, const char * argv[]){




}