//6518568 Ziyang Cheng

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

/* parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 5;
static int POP_SIZE = 100; //global parameters
int MAX_NUM_OF_GEN = 1000; //max number of generations
int MAX_TIME = 60;  //max amount of time permited (in sec)
float CROSSOVER_RATE = 0.8;
float MUTATION_RATE = 0.3;

struct solution_struct best_sln;  //global best solution
struct solution_struct temp_best; // temporary best

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
    if(prob != NULL)
    {
        if(prob->capacities !=NULL) {
            free(prob -> capacities);
        }
        if(prob->items!=NULL) {
            for(int j = 0; j < prob->n; j++) {
                if(prob->items[j].size != NULL) {
                    free(prob->items[j].size);
                }
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
    int counterer = 0;

    fp = fopen(filename, "r");
    fscanf(fp,"%d", &num_of_problems);
    
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k = 0; k < num_of_problems; k++)
    {
        fscanf(fp,"%d", &num_of_variables);
        fscanf(fp,"%d", &num_of_constraints);
        fscanf(fp,"%d", &optimal_solution_value);
        int n = num_of_variables, dim = num_of_constraints;
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory
        for(j = 0; j < n; j++) {   
            int p_value;
            fscanf(fp,"%d", &p_value);
            my_problems[k]->items[j].dim = dim;
            my_problems[k]->items[j].p = p_value;
            counterer++;
        }
        j = 0;
        for(j =0; j < dim; j++) {
            for(i = 0; i < n; i++) {
                int size;
                fscanf(fp,"%d", &size);
                my_problems[k]->items[i].size[j] = size;
                counterer++;
            }
        }
        i = 0;
        for(i = 0; i < dim; i++) {
            int capacities;
            fscanf(fp,"%d", &capacities);
            my_problems[k]->capacities[i] = capacities;
            counterer++;
        }
    }
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
    printf("%f\n",sln->objective);
    int number_of_items = (int)sln->prob->n;
    for(int i = 0; i < number_of_items; i++) {
        printf("%d ", sln->x[i]);
    }
    printf("\n");
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
}

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
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

void free_population(struct solution_struct* pop, int size) {
    for(int p=0; p<size; p++) {   
        if(pop[p].x != NULL && pop[p].cap_left != NULL) {
        free(pop[p].x);
        free(pop[p].cap_left);
        }
    }
}

//generate a new population by crossover
//crossover a fragment of chromosome
//heuristic method
//Wheel disk gambling selection operator
void cross_over(struct solution_struct* curt_pop, struct solution_struct* new_pop) {
    for(int p = 0; p < POP_SIZE/2; p++) {
        float crossover_rate = rand_01();
        if(crossover_rate < CROSSOVER_RATE){
            int i; 
            int counter1 = 0;
            int counter2 = 0;
            long int total_p = 0;
            int chromosome_length = curt_pop[p].prob->n-1;
            for(i = 0; i < POP_SIZE; i++) {
                // printf("obj: %f\n", curt_pop[i].objective);
                total_p += (int)curt_pop[i].objective;
                // printf("%ld\n", total_p);
            }
            long int exchange_node1 = (int)total_p * rand_01();
            // printf("e1: %ld ",exchange_node1);
            while(exchange_node1 > 0) {
                exchange_node1 -= curt_pop[counter1].objective;
                counter1++;
            }
            exchange_node1 = counter1-1;
            long int exchange_node2 = (int)total_p * rand_01();
            // printf("e2: %ld\n", exchange_node2);
            while(exchange_node2 > 0) {
                exchange_node2 -= curt_pop[counter2].objective;
                counter2++;
            }
            exchange_node2 = counter2-1;
            // printf("exchange_node: %ld %ld\n", exchange_node1, exchange_node2);
            int exchange_node = rand_int(0, chromosome_length);
            for(i = exchange_node; i < new_pop[p].prob->n; i++) {
                new_pop[exchange_node1].x[i] = curt_pop[exchange_node2].x[i];
            }
            evaluate_solution(&new_pop[p]);
        }
    }

    //normal crossover method
    // for(int p = 0; p < POP_SIZE; p++) {
    //     float exchange_rate = rand_01();
    //     int chromosome_length = new_pop[p].prob->n-1;
    //     if(exchange_rate < CROSSOVER_RATE){
    //         int exchange_pair = rand_int(0, POP_SIZE-1);
    //         int exchange_node = rand_int(0, chromosome_length);
    //         // int exchange_length = rand_int(exchange_node, chromosome_length);
    //         for(int i = exchange_node; i < new_pop[p].prob->n; i++) {
    //             new_pop[p].x[i] = curt_pop[exchange_pair].x[i];
    //         }
    //         evaluate_solution(&new_pop[p]);
    //     }
    // }    
}

//apply mutation to a population
//every chromosome can mutate randomly
//using dyadic mutation method
//using XRL
void mutation(struct solution_struct* new_pop) {
    float mutationRate = 0;
    //dyadic mutation method
    for(int p = 0; p < POP_SIZE; p++) {
        for(int i = 0; i < new_pop[p].prob->n; i++) {
            mutationRate = rand_01();
            if(mutationRate <= MUTATION_RATE) {
            int mutation_pair = rand_int(0, POP_SIZE-1);
                if(new_pop[p].x[i] == new_pop[mutation_pair].x[i]) {
                    if(new_pop[p].x[i] == 0) {
                        new_pop[p].x[i] = 1;
                    }
                    else {
                        new_pop[p].x[i] = 0;
                    }
                    evaluate_solution(&new_pop[p]);
                }
            }
        }
    }
    //normal mutation method
    // for(int p = 0; p < POP_SIZE; p++) {
    //     mutationRate = rand_01();
    //     if(mutationRate <= MUTATION_RATE) {
    //         int mutation_node = rand_int(0, new_pop[p].prob->n-1);
    //         if(new_pop[p].x[mutation_node] == 1) {
    //             new_pop[p].x[mutation_node] = 0;
    //         }
    //         else {
    //             new_pop[p].x[mutation_node] = 1;
    //         }
    //         evaluate_solution(&new_pop[p]);
    //     }
    // }
}

//modify the solutions that violate the capacity constraints
//by delecting p
//p/v is also in consider, but the results is almost the same
void feasibility_repair(struct solution_struct* pop) {
    //repair the feasibility and objective of the solution
    for(int p = 0; p < POP_SIZE; p++) {
        if(pop[p].feasibility < 0) {
            int counter, i;
            for(i = 0; i < pop[p].prob->n; i++) {
                if(pop[p].x[i] == 1){
                    counter++;
                }
            } 
            int p_value[counter];
            int p_index[counter];
            memset(p_value,0,sizeof(int)*counter);
            memset(p_index,0,sizeof(int)*counter);
            // for(int j = 0; j < counter; j++) {
            //     printf("value: %d, index: %d\n",pvalue[j],pindex[j]);
            // }
            counter = 0;
            for(i = 0; i < pop[p].prob->n; i++) {
                if(pop[p].x[i] == 1){
                    int value = pop[p].prob->items[i].p;
                    int index = i;
                    p_value[counter] = value;
                    p_index[counter] = index;
                    counter++;
                }
            }
            //insertion sort
            for(i = 1; i < counter; i++) {
                int temp_index, temp_value,j;
                temp_index = p_index[i];
                temp_value = p_value[i];
                for(j = i-1; j >= 0 && p_value[j] > temp_value; j-- ) {
                    p_value[j+1] = p_value[j];
                    p_index[j+1] = p_index[j];
                }
                p_value[j+1] = temp_value;
                p_index[j+1] = temp_index;
            }
        
            int remove_node = 0;
            while(pop[p].feasibility < 0) {
                pop[p].x[p_index[remove_node]] = 0;
                for(int i = 0; i < pop[p].prob->dim; i++) {
                    pop[p].cap_left[i] += pop[p].prob->items[p_index[remove_node]].size[i];
                }
                evaluate_solution (&pop[p]);
                remove_node++;
            }
        }   
    }
}


//local search
//search the first x0 and change to x1 if the capacity is enough
//what if there it is already the best solution?
void local_search_first_descent(struct solution_struct* pop) {
    int x_capacity = 1;
    //greedy search
    // for(int p = 0; p < POP_SIZE; p++) {
    //     int counter = 0;
    //     // printf("1:%f %d %d %d %d %d %d", pop[p].objective,pop[p].cap_left[0],pop[p].cap_left[1],pop[p].cap_left[2],pop[p].cap_left[3],pop[p].cap_left[4],counter);
        
    //     for(int i = 0; i < pop[p].prob->n; i++) {
    //         if(pop[p].x[i] == 0) {
    //             int* cap_left_test = pop[p].cap_left;
    //             for(int j = 0; j < pop[p].prob->dim; j++) {
    //                 cap_left_test[j] -= pop[p].prob->items[i].size[j];
    //                 if(cap_left_test[j] < 0) {
    //                     x_capacity = 0;
    //                 }
    //             }
    //             if(x_capacity == 0) {
    //                 x_capacity = 1;
    //                 continue;
    //             }
    //             else {
    //                 pop[p].x[i] = 1;
    //                 evaluate_solution(&pop[p]);
    //                 continue;
    //             }
    //         }
    //     }
    // }
    /*-------------------------------------*/
    // hill climbing
    for(int p = 0; p < POP_SIZE; p++) {
        int search_point = (int)pop[p].prob->n*rand_01();
        if(search_point - 40 > 0 && pop[p].prob->n - search_point > 40) {
            int idx = 0;
            int searchnode[41];
            memset(searchnode,0,sizeof(int)*41);
            int i;
            for(i = 0; i <= 20; i++) {
                int x_capacity =1;
                if(pop[p].x[search_point-i] == 0) {
                    int* cap_left_test = pop[p].cap_left; 
                    for(int j = 0; j < pop[p].prob->dim; j++) {
                        cap_left_test[j] -= pop[p].prob->items[i].size[j];
                    }
                    for(int k = 0; k < pop[p].prob->dim; k++) {
                        if(cap_left_test[k] < 0) {
                            x_capacity = 0;
                        }
                    }
                    if(x_capacity == 1) {
                        searchnode[i] = pop[p].prob->items[search_point-i].p;
                        continue;
                    }
                    x_capacity = 1;
                }
            }
            for(i = 1; i <= 20; i++) {
                int x_capacity =1;
                if(pop[p].x[search_point+i] == 0) {
                    int* cap_left_test = pop[p].cap_left; 
                    for(int j = 0; j < pop[p].prob->dim; j++) {
                        cap_left_test[j] -= pop[p].prob->items[i].size[j];
                    }
                    for(int k = 0; k < pop[p].prob->dim; k++) {
                        if(cap_left_test[k] < 0) {
                            x_capacity = 0;
                        }
                    }
                    if(x_capacity == 1) {
                        searchnode[i+20] = pop[p].prob->items[search_point+i].p;
                        continue;
                    }
                    x_capacity = 1;
                }
            }
            for(i = 0; i < 41; i++) {
                int max_p = 0;
                if(searchnode[i] > max_p) {
                    max_p = searchnode[i];
                    idx = i;
                }
            }
            pop[p].x[idx] = 1;
        }
        else if(search_point - 10 < 0) {
            int idx = 0;
            int searchnode[41];
            memset(searchnode,0,sizeof(int)*41);
            int i;
            for(i = 0; i <= 40; i++) {
                int x_capacity =1;
                if(pop[p].x[search_point+i] == 0) {
                    int* cap_left_test = pop[p].cap_left; 
                    for(int j = 0; j < pop[p].prob->dim; j++) {
                        cap_left_test[j] -= pop[p].prob->items[i].size[j];
                    }
                    for(int k = 0; k < pop[p].prob->dim; k++) {
                        if(cap_left_test[k] < 0) {
                            x_capacity = 0;
                        }
                    }
                    if(x_capacity == 1) {
                        searchnode[i] = pop[p].prob->items[search_point+i].p;
                        continue;
                    }
                    x_capacity = 1;
                }
            }
            for(i = 0; i < 41; i++) {
                int max_p = 0;
                if(searchnode[i] > max_p) {
                    max_p = searchnode[i];
                    idx = i;
                }
            }
            pop[p].x[idx] = 1;
        }
        else if(pop[p].prob->n - search_point < 40) {
            int idx = 0;
            int searchnode[41];
            memset(searchnode,0,sizeof(int)*41);
            int i;
            for(i = 0; i <= 40; i++) {
                int x_capacity =1;
                if(pop[p].x[search_point-i] == 0) {
                    int* cap_left_test = pop[p].cap_left; 
                    for(int j = 0; j < pop[p].prob->dim; j++) {
                        cap_left_test[j] -= pop[p].prob->items[i].size[j];
                    }
                    for(int k = 0; k < pop[p].prob->dim; k++) {
                        if(cap_left_test[k] < 0) {
                            x_capacity = 0;
                        }
                    }
                    if(x_capacity == 1) {
                        searchnode[i] = pop[p].prob->items[search_point-i].p;
                        continue;
                    }
                    x_capacity = 1;
                }
                for(i = 0; i < 41; i++) {
                int max_p = 0;
                if(searchnode[i] > max_p) {
                    max_p = searchnode[i];
                    idx = i;
                }
            }
            pop[p].x[idx] = 1;
            }
        }
    }
}

// replacement
void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop)
{
    //todo
    //sorting the top 100 population from curt_pop and new_pop
    //replace the top 100 population to new_pop

    struct solution_struct rep_pop[POP_SIZE*2];
    struct solution_struct temp_pop;
    int i,j,k;
 
    for(i = 0; i < POP_SIZE; i++) {
        copy_solution(&rep_pop[i], &curt_pop[i]);
    }
    for(j = 0; j < POP_SIZE; j++) {
        copy_solution(&rep_pop[POP_SIZE+j], &new_pop[j]);
    }
    //insertion sort
    for(k = 1; k < POP_SIZE*2; k++) {
        copy_solution(&temp_pop, &rep_pop[k]);
        int i = k-1;
        while(i > -1 && rep_pop[i].objective < temp_pop.objective) {
            free(rep_pop[i+1].cap_left);
            free(rep_pop[i+1].x);
            copy_solution(&rep_pop[i+1], &rep_pop[i]);
            i--;
            // free(rep_pop[i-1].cap_left);
            // free(rep_pop[i-1].x);
            // copy_solution(&rep_pop[i-1], &temp_pop);
        }
        free(rep_pop[i+1].cap_left);
        free(rep_pop[i+1].x);
        copy_solution(&rep_pop[i+1], &temp_pop);
    }

    //replacing the top100 to the new_pop
    for(int l = 0; l < POP_SIZE; l++) {
        free(new_pop[l].cap_left);
        free(new_pop[l].x);
        copy_solution(&new_pop[l], &rep_pop[l]);
        free(curt_pop[l].cap_left);
        free(curt_pop[l].x);
        copy_solution(&curt_pop[l], &rep_pop[l]);
    }
    free_population(rep_pop, POP_SIZE*2);
    if(temp_pop.x!=NULL && temp_pop.cap_left!=NULL) { 
        free(temp_pop.cap_left); 
        free(temp_pop.x);
    }
    // if(temp_best.objective < curt_pop[0].objective) {
    //     copy_solution(&temp_best, &curt_pop[0]);
    // }
    // else {
    //     copy_solution(&new_pop[0], &temp_best);
    //     copy_solution(&curt_pop[0], &temp_best);
    // }
    // printf("%f\n",temp_best.objective);
    // printf("%f\n", curt_pop[0].objective);
}

// update global best solution with best solution from pop if better
void update_best_solution(struct solution_struct* pop)
{
    // copy_solution(&best_sln, &temp_best);
    // for(int p = 0; p< POP_SIZE; p++) {
    //     if(pop[p].objective > best_sln.objective) {
    //         copy_solution(&best_sln, &pop[p]);
    //     }
    // }
    copy_solution(&best_sln, pop);
    // wait to be solved
    // output_solution(&best_sln, NULL);
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
    FILE *fp;
    int num_of_problems;
    fp = fopen(argv[2], "r");
    fscanf(fp,"%d", &num_of_problems);

    fp = fopen(argv[2], "r");
    fscanf(fp,"%d", &num_of_problems);
    fclose(fp);
    char data_file[50]={}, solution_file[50]={}, maxtime[50]={};
    if(argc < 3) {
        printf("Illegal numbers of arguments!");
        return 1;
    }
    else if(argc > 7) {
        printf("Too many arguments!");
        return 2;
    }
    else {
        for(int i = 1; i < argc; i= i+2) {
            if(strcmp(argv[i],"-s")==0) {
                strcpy(data_file, argv[i+1]);
            }
            else if(strcmp(argv[i],"-o")==0) {
                strcpy(solution_file, argv[i+1]);
            }
            else if(strcmp(argv[i],"-t")==0) {
                strcpy(maxtime, argv[i+1]);
                MAX_TIME =(int)atoi(maxtime);
            }
        }
    }
    //out put the file to the giving path
    // stdout = freopen(solution_file, "w", stdout);
    // printf("%s, %s, %d\n", data_file, solution_file, MAX_TIME);
    struct problem_struct** my_problems = load_problems(data_file);
    
    for(int k=0; k<num_of_problems; k++) {
        printf("This is the solution of problem %d\n",k);
        for(int run=0; run<NUM_OF_RUNS; run++) {
            // srand((unsigned)time(NULL));
            srand(RAND_SEED[run]);
            MA(my_problems[k]); //call MA
        }
        output_solution(&best_sln, NULL);
        free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL) { 
        free(best_sln.cap_left); 
        free(best_sln.x);
    } //free global

    return 0;
}