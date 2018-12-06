/* input is of form
 *
 * n m
 * axi bxj ... cxk = <min/max>
 * 
 * axi bxj ... cxk = C
 * axi bxj ... cxk = C
 * .
 * .
 * .
 *
 * where
 * i,j,k in Z+
 * a, b, c in R
 * C in R
 * n is a number of variables
 * m is a number of constraints
 *
 * i.e. all variables are implicitly non-negative
 *      coefficients of one MUST be defined
 *      variable indexes are integers (starting from one!)
 */

#include <stdio.h>  //printf, scanf, sscanf
#include <stdlib.h> //exit, malloc
#include <string.h> //strncmp
#include <ctype.h>  //tolower

#define BUFFER_SIZE 265
#define FACTOR_SIZE 8

struct FACTOR{
    int coefficient;
    int variable;
};

static double **readTable(int * const min, int * const m1, int * const n1);
static struct FACTOR readFactor(const char * const factor);
static void printTable(const double ** const t, const int n, const int m);
static void readLine(double ** const table, const int line);
static int isMinimization(void);
static void invertRow(double * const row, const int n);
static int *findBase(const double ** const table, const int m, const int n);
static int printVector(const int * const vector, int len);
static void unitize(double ** const table, const int m, const int n,
		    const int * const base);

void main(void){
    int min, m, n;
    double **table = readTable(&min, &m, &n);
    printTable((const double ** const)table, m, n);
    printVector(findBase((const double ** const)table, m, n), n);
    unitize(table, m,n, findBase((const double ** const) table, m,n));
    printTable((const double ** const) table, m,n);
}

static int printVector(const int * const vector, int len){
    for(int i = 0; i < len; i++)
	printf("%i\t", vector[i]);
    printf("\n");
}

/*
 * returns binary vector indicating variables in base
 *
 * presumes there is a valid base and m linearly independent colums
 *
 * if those conditions does not apply null pointer is returned instead
 */
static int *findBase(const double ** const table, const int m, const int n){
    int *base = (int *) malloc(sizeof(int) * n);
    int basic = 0;
    int nonzero;
    
    for(int j = 0; j < n; j++){
	nonzero = 0;
	for(int i = 1; i < m+1; i++){
	    
	    if(table[i][j] != 0) nonzero++;
	    if(nonzero > 1)      break;
	}
	if(nonzero <= 1) base[j] = 1;
	else             base[j] = 0;
    }
    for(int i = 0; i < n; i++)
	if(base[i] == 1)
	    basic++;
    
    if(basic == m) return base;
    else           return NULL;
}

void multiplyRow(double * const row, double coefficient, int len){
    for(int i = 0; i < len+1; i++)
	row[i] *= coefficient;
}

/*
 * if there is a base column with value other than one, corresponding 
 * row is unitized
 */
static void unitize(double ** const table, const int m, const int n,
		    const int * const base){
    for(int j =  0; j < n; j++){
	if(!base[j]) continue;
	for(int i = 1; i < m+1; i++)
	    if(table[i][j] != 0 && table[i][j] != 1)
		multiplyRow(table[i], 1/table[i][j], n);
    }
}

static double **readTable(int * const min, int * const m, int * const n){
    double **table;
    int l;
 
    //dimensions
    scanf("%d %d\n", m, n);

    //initialize two dimensional array
    table = (double **) malloc(sizeof(double *) * (*m+1));         //+1 for objective
    for(int i=0; i < *m +1; i++)
	table[i] = (double *) malloc(sizeof(double) * (*n+1));     //+1 for solutions

    //objective
    readLine(table, 0);
    invertRow(table[0], *n);
    *min = isMinimization();

    //input must specify min or max
    if(*min < 0) exit(-1);
    
    //constraints
    for(l = 1; l < *m+1; l++){
	readLine(table, l);
	if(l < *m+1 )    scanf("%lf\n", &table[l][*n]);
    }

    return table;
}

static void invertRow(double * const row, const int n){
    for(int i = 0; i < n+1; i++)
	row[i] = -row[i];
}

static int isMinimization(void){
    char buffer[BUFFER_SIZE];
    char *cp;

    scanf("%s", buffer);
    cp = buffer;

    for(; *cp; cp++) *cp = tolower(*cp);

    if(strncmp(buffer, "max", 3) == 0) return 0;
    if(strncmp(buffer, "min", 3) == 0) return 1; 

    return -1;
}

static void readLine(double ** const table, const int line){
    char factor[FACTOR_SIZE];
    struct FACTOR f;
    
    	while(scanf("%s", factor) == 1 && factor[0] != '='){
    	    f = readFactor(factor);
    	    table[line][f.variable-1] = (double)f.coefficient;
    	}
}

static struct FACTOR readFactor(const char * const factor){
    struct FACTOR f;
    sscanf(factor, "%ix%i", &f.coefficient, &f.variable);
    return f;
}

static void printTable(const double ** const t, const int m, const int n){
    for(int i = 0; i < m+1; i++){
	for(int j = 0; j < n+1; j++)
	    printf("%lf\t", t[i][j]);

	printf("\n");
    }
}
