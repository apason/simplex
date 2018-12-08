/* input is of form
 *
 * m n
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
#include <float.h>  //DBL_MAX

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
static void unitizeBase(double ** const table, const int m, const int n,
		    const int * const base);
static int findEntering(const double * const values, const int * const base, int n, int min);
static int findLeaving(const double ** const table, int entering, const int m, const int n);
static void pivotTable(table, m, n, a, b);
static void unitizeRowByColumn(double * const row, int n, int col);
static void multiplyRow(double * const row, double coefficient, int len);
static void addRow(double * const source, double * const subject, int n, double coefficient);
static double *getObjective(const double * const objective, int n);
static void printResult(const double * const objective, const double ** const table, int m, int n, int min);
static int isOptimal(const double * const x0, int n);
static void printObjective(const double * const objective, int n);
static void printValues(const double ** const table, int m, int n);

void main(void){
    int entering, leaving;                          // pivot column and pivot row
    int m, n;                                       // amount of constraints and variables
    int min;                                        // is minimization

    double  **table    = readTable(&min, &m, &n);   // simplex table 
    double  *objective = getObjective(table[0], n); // original objective function

    int *base;                                      // binary vector indicating whether base[i] is in base or not

    
    /* For the purpose of the simplex table, the objective */
    /* function is expressed in equation form              */
    invertRow(table[0], n);
    
    /* If the base is not unit length it must be unitized */
    unitizeBase(table, m,n, findBase((const double ** const) table, m,n));

    printf("Initial simplex table:\n");
    printTable((const double ** const)table, m, n);	

    /*
     * Main loop of simplex algorithm
     * 
     * Finds pivoting element, unitizes it and handles all row operations
     * neccessary to construct a new base
     */
    while(!isOptimal(table[0], n)){
	entering = findEntering(table[0], findBase((const double ** const)table, m, n), n, min);
	leaving = findLeaving((const double ** const) table,
			      findEntering(table[0], findBase((const double ** const)table, m, n), n, min)
			      , m, n);

	printf("Pivoting table[%d][%d]:\n\n", entering, leaving);
	pivotTable(table, m, n, entering, leaving);
	printTable((const double ** const)table, m, n);	

    }

    printf("Done!\n\n");
    printResult(objective, table, m, n, min);
}

/*
 * Prints left side expression of given equation (row).
 * Ignores coefficients of zero
 *
 * Precision of double may need adjustment
 */
static void printObjective(const double * const objective, int n){
    int first = 1;
    
    for(int j = 0; j < n; j++)
	if(objective[j] != 0)
	    if(first){
		printf("%.0lfx%i ",  objective[j], j+1);
		first = 0;
	    }
	    else
		printf("%+.0lfx%i ", objective[j], j+1);
}

/*
 * Prints the final results of the algoritm including
 * optimized result, and values of every variable
 */
static void printResult(const double * const objective, const double ** const table, int m, int n, int min){
    int * base = findBase(table, m, n);

    printObjective(objective, n);
    printf("= %.2lf, when\n\n", table[0][n]);
    printValues(table, m, n);
}

/* 
 * Prints value of every variable in format:
 * x1 = <value>
 */
static void printValues(const double ** const table, int m, int n){
    int *base = findBase((const double ** const)table, m, n);

    /* First print variables in the base */
    for(int j = 0; j < n; j++){
	if (base[j] == 1)
	    for(int i = 1; i < m+1; i++)
		if(table[i][j] == 1)
		    printf("x%d = %.2lf\n", j+1, table[i][n]);
    }
    
    /* Zero valued variables */
    for(int j = 0; j < n; j++)
	if(base[j] == 0)
	    printf("x%d = 0\n", j+1);
}

/*
 * Checks if the optimal value is reached
 */
static int isOptimal(const double * const x0, int n){
    
    for(int j = 0; j < n; j++)
	if(x0[j] < 0)
	    return 0;
    
    return 1;
}

/*
 * Copies given row of the table to a different memory region.
 * 
 * New memory region is returned.
 */
static double *getObjective(const double * const objective, int n){
    double * original = (double *)  malloc(sizeof(double) * n);

    for(int j = 0; j < n; j++)
	original[j] = objective[j];

    return original;
}


/*
 * Multiplies source row by coefficient and adds the result to subject.
 *
 * Row operation.
 */
static void addRow(double * const source, double * const subject, int n, double coefficient){
  
    for(int j = 0; j < n+1; j++)
	subject[j] += coefficient * source[j];
}

/*
 * Unitizes the (base) column of given row
 */
static void unitizeRowByColumn(double * const row, int n, int col){
    double coefficient = row[col];
    
    multiplyRow(row, 1/coefficient, n+1);
}

/*
 * Pivots the table by element table[a][b].
 */
static void pivotTable(const double ** const table, const int m, const int n,
		       const int a, const int b){
    
    unitizeRowByColumn(table[b], n, a);
    
    for(int i = 0; i < m+1; i++){
    	if(i == b)
	    continue;
	
	addRow(table[b], table[i], n, -table[i][a]);
    }
}

/*
 * Finds the leaving variable by taking minimum ratio of table[imin][n]/table[imin][entering].
 *
 * Returns ROW index of the leaving variable. (NOT the normal variable (name) index)
 */
static int findLeaving(const double ** const table, int entering, const int m, const int n){
    int imin;
    
    /* Initialize to first candidate */
    for(int i = 1; i < m; i++)
	if(table[i][entering] > 0)
	   imin = i;

    /* Find smallest ratio */
    for(int i = 1; i < m+1; i++){
	if(table[i][entering] <= 0)
	    continue;
	if(table[i][n]/table[i][entering] < table[imin][n]/table[imin][entering])
	    imin = i;
    }

    /* imin now corresponds to imin:th nonzero variable in base and imin:th row in simplex table */
    /* Return value is relative to the simplex row (and not the variable column)                 */
    return imin;
}

/*
 * Finds the entering variable by taking most negative/positive (in maximization/minimization)
 * variable which is not in the base. If there is no feasible leaving variable this method
 * may not work.
 *
 * Returns column index of the leaving variable
 */
static int findEntering(const double * const values, const int * const base, int n, int min){
    int index;

    /* Initialize to first nonbasic variable */
    for(int i = 0; i < n; i++)
	if(!base[i])
	    index = i;

    /* Find minumum/maximum */
    for(int i = 0; i < n; i++){
	if(base[i]) continue;
	if(!min){
	    /* Maximization */
	    if(values[i] < values[index])
		index = i;
	}
	else
	    /* Minimization */
	    if(values[i] > values[index])
		index = i;
    }
	
    return index;
}

/*
 * Prints the given (row) vector of the simplex table
 */
static int printVector(const int * const vector, int len){
    for(int i = 0; i < len; i++)
	printf("%i\t", vector[i]);
    printf("\n");
}

/*
 * Returns binary vector indicating variables in base.
 * Presumes there is a valid base and m linearly independent colums.
 * If those conditions does not apply null pointer is returned instead.
 */
static int *findBase(const double ** const table, const int m, const int n){
    int *base = (int *) malloc(sizeof(int) * n);
    int basic = 0;
    int nonzero;

    /* Scan every (column) vector corresponding to a variable. */
    /* If there is more than one non-zero element, the         */
    /* corresponding colum can not be in the base.             */
    for(int j = 0; j < n; j++){
	nonzero = 0;
	for(int i = 1; i < m+1; i++){
	    
	    if(table[i][j] != 0) nonzero++;
	    if(nonzero > 1)      break;
	}
	
	/* Set the column corresponding to a variable to zero  */
	/* or one depending whether the base condition is met. */
	if(nonzero <= 1) base[j] = 1;
	else             base[j] = 0;
    }

    /* Count the base to ensure the base is large enough.      */
    for(int i = 0; i < n; i++)
	if(base[i] == 1)
	    basic++;

    if(basic == m) return base;
    else           return NULL;
}

/*
 * Multiplies row by coefficient.
 */
static void multiplyRow(double * const row, double coefficient, int len){

    for(int i = 0; i < len+1; i++)
	row[i] *= coefficient;
}

/*
 * if there is a base column with value other than one, corresponding 
 * row is unitized
 */
static void unitizeBase(double ** const table, const int m, const int n, const int * const base){

    for(int j =  0; j < n; j++){
	if(!base[j]) continue;
	for(int i = 1; i < m+1; i++)
	    if(table[i][j] != 0 && table[i][j] != 1)
		multiplyRow(table[i], 1/table[i][j], n);
    }
}

/*
 * Reads the input and constructs a simplex table.
 * Format of the input is described in the beginning of this file.
 *
 * Returns the generated simplex table.
 */
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

/*
 * Multiplies row by -1.
 */
static void invertRow(double * const row, const int n){

    for(int i = 0; i < n+1; i++)
	row[i] = -row[i];
}

/*
 * Parses stdio for the type of the objective function. (min or max)
 *
 * Returns 1 for minimization, 0 for maximization and -1 otherwise.
 */
static int isMinimization(void){
    char buffer[BUFFER_SIZE];
    char *cp;

    scanf("%s", buffer);
    cp = buffer;

    /* Lowercases all letters in buffer */
    for(; *cp; cp++) *cp = tolower(*cp);

    if(strncmp(buffer, "max", 3) == 0) return 0;
    if(strncmp(buffer, "min", 3) == 0) return 1; 

    return -1;
}

/*
 * Reads one expression from the input. This can be either the
 * objective function or left side of the constraint equation.
 */
static void readLine(double ** const table, const int line){
    char factor[FACTOR_SIZE];
    struct FACTOR f;
    
    while(scanf("%s", factor) == 1 && factor[0] != '='){
	f = readFactor(factor);
	table[line][f.variable-1] = (double)f.coefficient;
    }
}

/*
 * Reads one factor (in form of axi, where a is coefficient 
 * and i is the variable index) from stdio.
 *
 * Returns struct FACTOR corresponding this input.
 */
static struct FACTOR readFactor(const char * const factor){
    struct FACTOR f;
    
    sscanf(factor, "%ix%i", &f.coefficient, &f.variable);
    
    return f;
}

/* 
 * Prints the simplex table.
 * Precision of the output may need adjustment
 */
static void printTable(const double ** const t, const int m, const int n){

    for(int i = 0; i < m+1; i++){
	for(int j = 0; j < n+1; j++)
	    printf("%.2lf\t", t[i][j]);

	printf("\n");
    }
    printf("\n");
}
