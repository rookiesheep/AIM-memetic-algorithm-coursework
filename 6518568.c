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
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    // printf("rand_range= %d \n", val);
    return val;
}

struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
};

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
};

void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}

struct problem_struct** load_problems(const char *filename)
{
    FILE *fp;
    int i, j, k;
    int num_of_problems, num_of_variables, num_of_constraints,optimal_solution_value;
    int counter = 0;

    fp = fopen(filename, "r");
    fscanf(fp,"%d", &num_of_problems);
    
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k = 0; k < num_of_problems; k++)
    {
        fscanf(fp,"%d", &num_of_variables);
        fscanf(fp,"%d", &num_of_constraints);
        fscanf(fp,"%d", &optimal_solution_value);
        // printf("%d\n", num_of_variables);
        // printf("%d\n", num_of_constraints);
        // printf("%d\n", optimal_solution_value);
        int n = num_of_variables, dim = num_of_constraints;
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j = 0; j < n; j++) {   
            int p_value;
            fscanf(fp,"%d", &p_value);
            // printf("%d\n", p_value);
            my_problems[k]->items[j].dim = dim;
            my_problems[k]->items[j].p = p_value;
            counter++;
        }
        // printf("p: %d\n", counter);
        j = 0;
        for(j =0; j < dim; j++) {
            for(i = 0; i < n; i++) {
                int size;
                fscanf(fp,"%d", &size);
                // printf("%d\n", size);
                my_problems[k]->items[i].size[j] = size;
                counter++;
            }
            // printf("p: %d\n", counter);
        }
        i = 0;
        for(i = 0; i < dim; i++) {
            int capacities;
            fscanf(fp,"%d", &capacities);
            // printf("%d\n", capacities);
            my_problems[k]->capacities[i] = capacities;
            counter++;
        }
        // printf("p: %d\n", counter);
    }
    // printf("total：%d", counter);
    printf("initialzed successfully!");
    fclose(fp);
    return my_problems;
}

void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; 
    sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++) {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++) {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1; //exceeding capacity
                return;
            }
        }
    }
    if(sln->feasibility>0) {
        for(int j=0; j<sln->prob->n; j++) {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, const char* out_file) {
    // printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
    printf("%f",sln->objective);
    int number_of_items = (int)sln->prob->n;
    for(int i = 0; i < number_of_items; i++) {
        printf("%d ", sln->x[i]);
    }
}

//intialise the population with random solutions
void init_population(struct problem_struct* prob, struct solution_struct* pop)
{
    for(int p=0; p<POP_SIZE; p++)
    {
        pop[p].prob = prob;
        pop[p].x = malloc(sizeof(int)*prob->n);
        pop[p].cap_left = malloc(sizeof(int)*prob->dim);
        for(int j=0; j<prob->n; j++) {
            pop[p].x[j] = 0;
        }
        for(int i=0; i<prob->dim; i++) {
            pop[p].cap_left[i]=prob->capacities[i];
        }
        /* create a random initial x that is feasible */
        int j=rand_int(0, prob->n-1);
        while(true) {
            while(pop[p].x[j]==1) {
                j=rand_int(0, prob->n-1); //select an unpacked item randomly
            }
            pop[p].x[j]=1;
            bool can_pack=true;
            for(int i=0; i< prob->dim; i++) {
                pop[p].cap_left[i] -= prob->items[j].size[i];
                if(pop[p].cap_left[i] <0) {
                    can_pack=false;
                }
            }
            if(!can_pack) {   
                //unpack item i
                //printf("packing item %d failed. random initialisation stoped.\n", j);
                pop[p].x[j]=0;
                for(int i=0; i< prob->dim; i++)
                    pop[p].cap_left[i] += prob->items[j].size[i];
                break;
            }
        }
        evaluate_solution (&pop[p]);
    }
    // for(int p=0; p<POP_SIZE; p++)
    // {
    //     output_solution(&pop[p], NULL);
    // }
}

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}

void free_population(struct solution_struct* pop, int size)
{
    for(int p=0; p<size; p++)
    {
        free(pop[p].x);
        free(pop[p].cap_left);
        pop[p].objective = 0;
        pop[p].feasibility = 0;
        free_problem(pop[p].prob);
    }
}

//generate a new population by crossover
//crossover a fragment of chromosome
void cross_over(struct solution_struct* curt_pop, struct solution_struct* new_pop)
{
    //exchange from the node x
    int pair = POP_SIZE/2;
    for(int p = 0; p < pair; p++) {
        float exchange_rate = rand_01();
        int chromosome_length = new_pop[p].prob->n-1;
        if(exchange_rate > CROSSOVER_RATE){
            int exchange_node = rand_int(0, chromosome_length);
            for(int i = exchange_node; i < chromosome_length; i++) {
                int temp = new_pop[p].x[i];
                new_pop[p].x[i] = new_pop[p+pair].x[i];
                new_pop[p+pair].x[i] = temp;
            }
            evaluate_solution(&new_pop[p]);
            evaluate_solution(&new_pop[p+pair]);
        }
    }
    printf("cross_over successfully!");
}

//apply mutation to a population
//every chromosome can mutate randomly
void mutation(struct solution_struct* pop)
{
    int mutationRate = 0;
    for(int p = 0; p < POP_SIZE; p++) {
        int j = pop[p].prob->n;
        for(int i = 0; i < j; i++) {
            int mutationRate = rand_01();
            if(mutationRate > MUTATION_RATE) {
                if(pop[p].x[i] ==1) {
                    pop[p].x[i] = 0;
                }
                else {
                    pop[p].x[i] = 1;
                }
            }
        }
        evaluate_solution(&pop[p]);
    }
    printf("mutation successfully!");
}

//modify the solutions that violate the capacity constraints
//by delecting p
void feasibility_repair(struct solution_struct* pop) {
    //repair the feasibility and objective of the solution
    for(int p = 0; p < POP_SIZE; p++) {
        evaluate_solution (&pop[p]);
        if(pop[p].feasibility < 0) {
            for(int i = 0; i < pop[p].prob->n; i++) {
                if(pop[p].x[i] == 1) {
                    pop[p].x[i] = 0;
                    for(int j = 0; j < pop[p].prob->dim; j++) {
                        pop[p].cap_left[j] += pop[p].prob->items[i].size[j];
                    }   
                }
                evaluate_solution (&pop[p]);
                if(pop[p].feasibility > 0) {
                    break;
                }
            }
        }
    }
    printf("feasibility_repair successfully!");
}

//local search
//search the first x0 and change to x1 if the capacity is enough
void local_search_first_descent(struct solution_struct* pop)
{
    for(int p = 0; p < POP_SIZE; p++) {
        for(int i = 0; i < pop[p].prob->n; i++) {
            if(pop[p].x[i] == 0) {
                int* cap_left_test = pop[p].cap_left; 
                int x_capacity = 1;
                for(int j = 0; j < pop[p].prob->dim; j++) {
                    cap_left_test[j] -= pop[p].prob->items[i].size[j];
                }
                for(int k = 0; k < pop[p].prob->dim; k++) {
                    if(cap_left_test < 0){
                        int x_capacity = 0;
                    }
                }
                if(x_capacity == 1){
                    pop[p].x[i] = 1;
                    break;
                }
            }
        }
        evaluate_solution(&pop[p]);
    }
        printf("local_search_first_descent successfully!");
}

//replacement
void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop)
{
    //todo
    //sorting the top 100 population from curt_pop and new_pop
    //replace the top 100 population to new_pop

    struct solution_struct rep_pop[POP_SIZE*2];
    struct solution_struct temp_pop[0];
    for(int p = 0; p<POP_SIZE*2; p++){
     rep_pop[p].prob = curt_pop[0].prob;
    rep_pop[p].x = malloc(sizeof(int)*curt_pop[0].prob->n);
    rep_pop[p].cap_left = malloc(sizeof(int)*curt_pop[0].prob->dim);
    for(int j=0; j<curt_pop[0].prob->n; j++)    rep_pop[p].x[j] = 0;
    for(int i=0; i<curt_pop[0].prob->dim; i++)  rep_pop[p].cap_left[i]=curt_pop[0].prob->capacities[i];
    }
    
    temp_pop[0].prob = curt_pop[0].prob;
    temp_pop[0].x = malloc(sizeof(int)*curt_pop[0].prob->n);
    temp_pop[0].cap_left = malloc(sizeof(int)*curt_pop[0].prob->dim);
    for(int jt=0; jt<curt_pop[0].prob->n; jt++)    temp_pop[0].x[jt] = 0;
    for(int it=0; it<curt_pop[0].prob->dim; it++)  temp_pop[0].cap_left[it]=curt_pop[0].prob->capacities[it];

    for(int i = 0; i < POP_SIZE; i++) {
        copy_solution(&rep_pop[i], &curt_pop[i]);
        //if not printout error
    }
    for(int j = 0; j < POP_SIZE; j++) {
        copy_solution(&rep_pop[POP_SIZE+j], &new_pop[j]);
        //if not printout error
    }
    //insertion sort
    for(int k = 1; k < POP_SIZE*2; k++) {
        copy_solution(&temp_pop[0], &rep_pop[k]);
        for(int i = k; i > 0 && rep_pop[k-1].objective < temp_pop[0].objective; i--) {
            copy_solution(&rep_pop[k], &rep_pop[k-1]);
            copy_solution(&rep_pop[k-1], &temp_pop[0]);
        }
    }
    //replacing the top100 to the new_pop
    for(int l = 0; l < POP_SIZE; l++) {
        copy_solution(&new_pop[l], &rep_pop[l]);
    }
    free_population(rep_pop, POP_SIZE*2);
    // free_population(&temp_pop[0], 1);
}

//update global best solution with best solution from pop if better
void update_best_solution(struct solution_struct* pop)
{
    best_sln = pop[0];
    output_solution(&best_sln, NULL);
}

//memetic algorithm
int MA(struct problem_struct* prob)
{
    struct solution_struct curt_pop[POP_SIZE];
    struct solution_struct new_pop[POP_SIZE];
    init_population(prob, curt_pop);
    init_population(prob, new_pop);
    int gen=0;
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    while(gen<MAX_NUM_OF_GEN && time_spent < MAX_TIME)
    {
        cross_over(curt_pop, new_pop);
        mutation(new_pop);
        feasibility_repair(new_pop);
        local_search_first_descent(new_pop);
        replacement(curt_pop, new_pop);
        gen++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }
    
    // update_best_solution(curt_pop);
    
    
    // free_population(curt_pop, POP_SIZE);
    // free_population(new_pop, POP_SIZE);
    
    return 0;
}

int main(int argc, const char * argv[]){
    FILE *fp;
    int num_of_problems;
    fp = fopen(argv[2], "r");
    fscanf(fp,"%d", &num_of_problems);

    int maxRunningTime = atoi(argv[6]);
    printf("maxRunningTime: %d\n", maxRunningTime);
    fp = fopen(argv[2], "r");
    fscanf(fp,"%d", &num_of_problems);
    fclose(fp);


    struct problem_struct** my_problems = load_problems(argv[2]);
    for(int k=0; k<num_of_problems; k++) {
        printf("This is the solution of problem %d",k+1);
        for(int run=0; run<NUM_OF_RUNS; run++) {
            srand(RAND_SEED[run]);
            MA(my_problems[k]); //call MA
            output_solution(&best_sln, argv[4]);
        }
        free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global


    return 0;
}