/* simplex.c - implements the simplex algorithm */

/* Author: Arttu Kilpinen                       */

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
#include <stdlib.h> //exit, malloc, free
#include <string.h> //strncmp
#include <ctype.h>  //tolower

#define BUFFER_SIZE 265
#define FACTOR_SIZE 8

#define MINIMIZE 1

#define EPSILON 0.000000001

struct FACTOR{
    int coefficient;
    int variable;
};

/* Input functions */
static double * const * const readTable (int * const min, int * const m1, int * const n1);
static void readLine(double * const * const table, const int line);
static struct FACTOR readFactor(const char * const factor);
static int isMinimization(void);

/* Output functions */
static void printTable(const double * const * const t, const int n, const int m, const double const * original_objective);
static void printValues(const double * const * const table, const int m, const int n);
static void printObjective(const double * const objective, const int n);
static int printVector(const int * const vector, const int len);
static void printResult(const double * const objective, const double * const * const table,
			const int m, const int n, const int min);

/* Main simplex functions */
static void pivotTable(const double * const * const table, const int m, const int n, const int a, const int b);
static void addArtificialVariables(double * const * const table, const int m, const int n, const int extra);
static double * const * firstPhase(double ** const table, const int m, const int n, const int extra);
static int findEntering(const double * const values, const int * const base, const int n, const int min);
static int findLeaving(const double * const * const table, const int entering, const int m, const int n);
static int lex(const double * const x1, const double * const x2, const int n, const int entering);
static int simplex(double * const * table, const int m, const int n, const int min, double * const original);
static int * findBase(const double * const * const table, const int m, const int n);
static int isOptimal(const double * const x0, const int n, const int min);
static void useOnlyNonBaseInObjective(const double * const * const table, const int * const base,
				      const int m, const int n);

/* Helper functions */
static void addRow(const double * const source, double * const subject, const int n, const double coefficient);
static void unitizeBase(double * const * const table, const int m, const int n, const int * const base);
static void multiplyRow(double * const row, const double coefficient, const int len);
static void unitizeRowByColumn(double * const row, const int n, const int col);
static double * getCopyOfRow(const double * const objective, const int n);
static int baseSize(const int * const base, const int n);
static void invertRow(double * const row, const int n);
static void freeTable(double ** table, const int m);
static void copyToOriginal(double * const * const table, const double * const * const new_table,
			   const int m, const int n, const int extra);


int main(char *argc[], int argv){
    int m, n;                                                     // amount of constraints and variables
    int min;                                                      // is minimization

    double * const * table     = readTable(&min, &m, &n);         // simplex table 
    double * const   objective = getCopyOfRow(table[0], n);       // original objective function

    const int *base;                                              // binary vector indicating base variables

    /* For the purpose of the simplex table, the objective */
    /* function is expressed in equation form.             */
    invertRow(table[0], n); 
    
    /* If the base is not unit length it must be unitized. */
    base = findBase((const double * const * const)table, m, n);
    unitizeBase(table, m,n, base);
    
    /* If there is no suitable base in original simplex    */
    /* table, two phase method is used.                    */
    if(baseSize(base, n) < m){
	table = firstPhase((double ** const)table, m, n, m-baseSize(base, n));

	/* If first phase does not find a valib base the   */
	/* original instance is infeasible.                */
	if(table == NULL){
	    free((void *)base);
	    free(objective);    
	    return 1;
	}
    }

    if(simplex(table, m, n, min, NULL))
	printResult(objective, (const double * const * const)table, m, n, min);

    free((void *)base);
    freeTable((double **)table, m);
    free(objective);    
}

/* *********************************************** */
/* *******    Main simplex functions    ********** */
/* *********************************************** */


/*
 * Handles the initialization of the first phase of two phase simplex:
 * Creates new (bigger) simplex table and copies the old information 
 * into it. Modifies new colums and objective so it can be solved with
 * "normal" simplex algorithm.
 */
static double * const * firstPhase(double ** const table, const int m, const int n, const int extra){
    double ** const new_table          = (double ** const) malloc (sizeof(double *) * (m+1));
    double *  const original_objective = (double * const) malloc(sizeof(double) * (n + extra + 1));
    
    /* Create the first row and allocate it to zero. Row length is original + extra +1. */
    new_table[0] = (double *) malloc(sizeof(double) * (n + extra +1));
    memset(new_table[0], 0, sizeof(double) * (n + extra +1));

    /* Create copy of the first row of the original simplex table.                      */
    for(int i = 0; i < n; i++) original_objective[i] = table[0][i];
    original_objective[n+extra]                      = table[0][n];

    /* Create other rows corresponding to base variables. Copy old data to a new table. */
    for(int i = 1; i < m+1; i++){
	new_table[i] = (double *) malloc(sizeof(double) * (n + extra +1));
	memset(new_table[i], 0, sizeof(double) * (n + extra +1));

	/* Copy variables and solution column. The space for extra is filled later.     */
	for(int j = 0; j < n; j++)
	    new_table[i][j]   = table[i][j];
	new_table[i][n+extra] = table[i][n];
    }

    /* Prepare simplex table and use regular one phase simplex to find starting bfs.    */
    addArtificialVariables(new_table, m, n, extra);
    printf("Initial first phase simplex table:\n\n");
    printTable((const double ** const)new_table, m, n+extra, original_objective);
    simplex(new_table, m, n+extra, MINIMIZE, original_objective);

    /* Check if original table has a feasible solution.                                 */
    if(new_table[0][n+extra] >  EPSILON ||              
       new_table[0][n+extra] < -EPSILON  ){
	
	freeTable((double **)new_table, m);
	freeTable((double **)    table, m);
	free(original_objective);
	printf("Infeasible!\n");
	return NULL;
    }

    free(new_table[0]);
    new_table[0] = original_objective;

    /* Copy data back to original table and free the new one with artificial variables. */
    copyToOriginal(table, (const double * const * const)new_table, m, n, extra);
    freeTable(new_table, m);
    
    printf("\n\nInitial simplex table:\n\n");
    printTable((const double ** const)table, m, n, NULL);

    return table;
}

/*
 * Regular "one phase" simplex.
 */
static int simplex(double * const * table, const int m, const int n, const int min, double * const oo){
    int entering, leaving;                          // pivot column and pivot row
    int *base;                                      // binary vector indicating whether base[i] is in base or not
    
    /* If the base is not unit length it must be unitized  */
    base = findBase((const double ** const) table, m, n);
    unitizeBase(table, m,n, base);

    /* For two phase simplex */
    useOnlyNonBaseInObjective((const double ** const)table, base, m, n);
    free(base);

    /* Print initial table.  */
    printf("\n\nCost row modified to use only nonbasis variables:\n\n");
    printTable((const double ** const)table, m, n, oo);

    /*
     * Main loop of simplex algorithm
     * 
     * Finds pivoting element, unitizes it and handles all row operations
     * neccessary to construct a new base
     */
    while(!isOptimal(table[0], n, min)){
	base = findBase((const double ** const)table, m, n);
	entering = findEntering(table[0], base, n, min);
	leaving = findLeaving((const double ** const) table, findEntering(table[0], base, n, min),  m, n);

	if(leaving == 0){
	    free(base);
	    printf("Unbounded!\n");
	    return 0;
	}

	printf("\n\nPivoting table[%d][%d]:\n\n", entering, leaving);
	pivotTable((const double ** const)table, m, n, entering, leaving);
	if(oo) addRow(table[0], oo, n, -oo[entering]);

	printTable((const double ** const)table, m, n, oo);
	free(base);
    }

    return 1;
}

/*
 * Returns binary vector indicating variables in base.
 * Presumes there is a valid base and m linearly independent colums.
 * Base size should be checked after a call to findBase!
 */
static int *findBase(const double * const * const table, const int m, const int n){
    int *base = (int *) malloc(sizeof(int) * n);
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
	if(nonzero <= 1){
	    /* Find the element and ensure it is positive      */
	    for(int i = 1; i < m+1; i++)
		if(table[i][j] > 0)
		    base[j] = 1;
	}
	else        base[j] = 0;
    }

    return base;
}

/*
 * Finds the leaving variable by taking minimum ratio of table[imin][n]/table[imin][entering].
 * If there is two equal minimum ratios, lexicographic method is used.
 *
 * Returns ROW index of the leaving variable. (NOT the normal variable (name) index)
 * If the instance is unbounded, imin is set to zero.
 */
static int findLeaving(const double * const * const table, const  int entering, const int m, const int n){
    int imin = 0;  //initialize to invalid value
    
    /* Initialize to first candidate */
    for(int i = 1; i < m; i++)
	if(table[i][entering] > 0)
	   imin = i;

    /* Find smallest ratio */
    if(imin > 0)
	for(int i = 1; i < m+1; i++){
	    if(table[i][entering] <= 0)
		continue;
	    if(table[i][n]/table[i][entering] < table[imin][n]/table[imin][entering] - EPSILON)      // smaller
		imin = i;
	    else if(table[i][n]/table[i][entering] < table[imin][n]/table[imin][entering] + EPSILON) // equal
		imin = lex(table[i], table[imin], n, entering) == 0 ? i : imin;
	}

    /* imin now corresponds to imin:th nonzero variable in base and imin:th row in simplex table */
    /* Return value is relative to the simplex row (and not the variable column)                 */
    return imin;
}

/*
 * Finds the entering variable by taking most negative/positive (in maximization/minimization)
 * variable which is not in the base.
 *
 * Returns column index of the leaving variable
 */
static int findEntering(const double * const values, const int * const base, const int n, const int min){
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
 * Artificial variables comes to base. This function calculates the unit (column) vectors
 * for those variables. (the place where to put one in that column)
 * Reserved vector describes rows where new base units can be put to.
 */
static void addArtificialVariables(double * const * const table, const int m, const int n, const int extra){
    int * base            = findBase((const double ** const)table, m, n+extra);
    int * const reserved  = (int * const) malloc(sizeof(int) * m);          

    memset(reserved, 0, sizeof(int) * m);    

    /* Find rows with unit element corresponding to a base variabel         */
    for(int j = 0; j < n+extra; j++)
    	if(base[j])
    	    for(int i = 1; i < m+1; i++)
    		if(table[i][j])
    		    reserved[i] = 1;

    /* Distribute units corresponding to a artificial elements to free fows */
    for(int j = n; j < n+extra; j++)
    	for(int i = 1; i < m+1; i++)
    	    if(!reserved[i]){
    		table[i][j] = 1;
    		table[0][j] = -1;
    		reserved[i] = 1;
    		break;
    	    }

    free(base);
    free(reserved);
}

/*
 * Divides rows with entering column and calculates lexicographical ordering for resultung rows.
 * If the first row is lexicographically smaller then 0 is returned, 1 otherwise.
 */
int lex(const double * const x1, const double * const x2, const int n, const int entering){
    /* Copy rows from the original simplex table.        */
    double * const xa = getCopyOfRow(x1, n);
    double * const xb = getCopyOfRow(x2, n);

    /* "Pivot" rows                                      */
    unitizeRowByColumn(xa, n, entering);
    unitizeRowByColumn(xb, n, entering);

    /* Find lexicographical orde and return if possible. */
    for(int j = 0; j < n; j++)
	if(xa[j] < xb[j] - EPSILON)                      // truly smaller
	    return 0;
	else if (xa[j] > xb[j] + EPSILON)                // truly greater
	    return 1;

    free(xa);
    free(xb);
    
    return 1;                                            // return value here does not matter
}

/*
 * Pivots the table by element table[a][b].
 */
static void pivotTable(const double * const * const table, const int m, const int n, const int a, const int b){

    unitizeRowByColumn((double * const)table[b], n, a);

    for(int i = 0; i < m+1; i++){
    	if(i == b)
	    continue;

	addRow((double * const)table[b], (double * const)table[i], n, -table[i][a]);
    }
}

/*
 * This function  guarantees that there is no 
 * base variables in the objective function.
 */
static void useOnlyNonBaseInObjective(const double * const * const table, const int * const base, const int m,
				      const int n){
    for(int j = 0; j < n; j++)
	if(base[j])
	    if(table[0][j])
		for(int i = 1; i < m+1; i++)
		    if(table[i][j])
			addRow((double * const)table[i], (double * const)table[0], n, -table[0][j]);
}

/*
 * Checks if the optimal value is reached.
 */
static int isOptimal(const double * const x0, const int n, const int min){

    for(int j = 0; j < n; j++)
	if(min ? x0[j] > EPSILON : x0[j] < -EPSILON)
	    return 0;
    
    return 1;
}

/* *********************************************** */
/* ***********     Input functions     *********** */
/* *********************************************** */

/*
 * Reads the input and constructs a simplex table.
 * Format of the input is described in the beginning of this file.
 *
 * Returns the generated simplex table.
 */
static double * const * const readTable(int * const min, int * const m, int * const n){
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
static void readLine(double * const * const table, const int line){
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

/* *********************************************** */
/* ***********    Output functions     *********** */
/* *********************************************** */


/* 
 * Prints value of every variable in format:
 * x1 = <value>
 *
 * Precision of double may need adjustment.
 */
static void printValues(const double * const * const table, const int m, const int n){
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

    free(base);
}

/*
 * Prints left side expression of given equation (row).
 * Ignores coefficients of zero
 *
 * Precision of double may need adjustment
 */
static void printObjective(const double * const objective, const int n){
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
 * optimized result and values of every variable.
 *
 * Precision of double may need adjustment.
 */
static void printResult(const double * const objective, const double * const * const table,
			const int m, const int n, const int min){
    
    int * base = findBase(table, m, n);

    printObjective(objective, n);
    printf("= %.2lf, when\n\n", table[0][n]);
    printValues(table, m, n);

    free(base);
}

/* 
 * Prints the simplex table.
 * Precision of the output may need adjustment
 */
static void printTable(const double * const * const t, const int m, const int n, const double * const original_objective){

    int * const base = findBase(t, m, n);
    int bi = 0;

    if(original_objective != NULL){
	printf(" z | ");
	for(int j = 0; j < n+1; j++){
	    if(j == n) printf(" | ");
	    printf("%.2lf\t", original_objective[j]);
	}

	printf("\n");
	for(int j = 0; j < (n+1)*9+1; j++)
	    printf("-");
	printf("\n");
    }

    printf(" w | ");
    for(int j = 0; j < n+1; j++){
	if(j == n) printf(" | ");
	printf("%.2lf\t", t[0][j]);
    }

    printf("\n");
    for(int j = 0; j < (n+1)*9+1; j++)
	printf("-");
    printf("\n");
    
    for(int i = 1; i < m+1; i++){
	if(i > 0){
	    for(; !base[bi]; bi++);
	    printf("x%d", ++bi);
	}


	for(int j = 0; j < n+1; j++){
	    if(j == n || j == 0) printf(" | ");
	    printf("%.2lf\t", t[i][j]);
	}

	printf("\n");
    }
    printf("\n");

    free(base);
}

/*
 * Prints the given (row) vector of the simplex table
 */
static int printVector(const int * const vector, const int len){
    for(int i = 0; i < len; i++)
	printf("%i\t", vector[i]);
    printf("\n");
}



/* *********************************************** */
/* ***********    Helper functions     *********** */
/* *********************************************** */


/*
 * Copies phase1 simplex table (larger) to original (smaller) simplex table.
 * The first row is copied back in the main function.
 */
void copyToOriginal(double * const * const table, const double * const * const new_table,
		    const int m, const int n, const int extra){

    for(int i = 0; i < m+1; i++){
	for(int j = 0; j < n; j++)
	    table[i][j] = new_table[i][j];
	table[i][n] = new_table[i][n+extra];
    }
}
	    
/*
 * Returns the amount of ones in the given binary vector.
 */
static int baseSize(const int * const base, const int n){
    int sum = 0;

    for(int j = 0; j < n; j++)
	if(base[j] == 1)
	    sum++;

    return sum;
}

/*
 * Frees the memory allocated to the simplex table
 */
static void freeTable(double ** table, const int m){

    for(int i = 0; i < m+1; i++)
	free(table[i]);

    free(table);
}

/*
 * Copies given row of the table to a different memory region.
 * 
 * New memory region is returned.
 */
static double * getCopyOfRow(const double * const objective, const int n){
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
static void addRow(const double * const source, double * const subject, const int n, const double coefficient){
  
    for(int j = 0; j < n+1; j++)
	subject[j] += coefficient * source[j];
}

/*
 * Unitizes the (base) column of given row.
 * Make sure that row[col] is nonzero!
 */
static void unitizeRowByColumn(double * const row, const int n, const int col){
    double coefficient = row[col];
    
    multiplyRow(row, 1/coefficient, n+1);
}

/*
 * Multiplies row by coefficient.
 */
static void multiplyRow(double * const row, const double coefficient, const int len){

    for(int i = 0; i < len; i++)
	row[i] *= coefficient;
}

/*
 * if there is a base column with value other than one, corresponding 
 * row is unitized
 */
static void unitizeBase(double * const * const table, const int m, const int n, const int * const base){

    for(int j =  0; j < n; j++){
	if(!base[j]) continue;
	for(int i = 1; i < m+1; i++)
	    if(table[i][j] != 0 && table[i][j] != 1)
		multiplyRow(table[i], 1/table[i][j], n);
    }
}

/*
 * Multiplies row by -1.
 */
static void invertRow(double * const row, const int n){

    for(int i = 0; i < n+1; i++)
	row[i] = -row[i];
}

