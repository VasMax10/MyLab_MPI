#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int ProcNum = 0;      // Number of available processes 
int ProcRank = 0;     // Rank of current process


// A structure to represent a Set for union-find
struct Set
{
	int parent;
	int rank;
};

// Function prototypes for union-find (These functions are defined after boruvkaMST() )
int find(struct Set sets[], int i);
void Union(struct Set sets[], int x, int y);



// Function for random definition of matrix and vector elements
void DataInitialization(int* pMatrix, int V) {//, int E, int SizeR, int SizeC) {
	int i, j;  // Loop variables

	int count = 0;
	srand(time(NULL));
	for (i = 0; i < V; i++)
		for (j = i + 1; j < V; j++)
		{
			int a, b, c;
			a = i;
			b = j;
			//c = b + a;
			c = rand() % V + 1;
			pMatrix[3 * count] = a;
			pMatrix[3 * count + 1] = b;
			pMatrix[3 * count + 2] = c;
			count++;
		}
}

// Function for memory allocation and data initialization
void ProcessInitialization(int*& pMatrix, int*& pResult, int*& pProcRows, int*& pProcResult, int& SizeR, int& SizeC , int& V, int& E, int& RowNum) {
	
	int RestRows; // Number of rows, that haven’t been distributed yet
	int i;        // Loop variable

	setvbuf(stdout, 0, _IONBF, 0);

	if (ProcRank == 0) {
		do {
			printf("\nEnter size of the initial objects: ");
			scanf_s("%d", &V);
			E = V * (V - 1) / 2;
			SizeR = E;
			SizeC = 3;
			if (SizeR < ProcNum) {
				printf("Size of the objects must be greater than number of processes! \n ");
			}
		} while (SizeR < ProcNum);
	}
	MPI_Bcast(&V, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&E, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&SizeR, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&SizeC, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Determine the number of matrix rows stored on each process
	RestRows = SizeR;
	for (i = 0; i < ProcRank; i++)
		RestRows = RestRows - RestRows / (ProcNum - i);
	RowNum = RestRows / (ProcNum - ProcRank);

	// Memory allocation
	pResult = new int[SizeR];
	pProcRows = new int[RowNum * SizeC];
	pProcResult = new int[RowNum];

	// Obtain the values of initial objects elements
	if (ProcRank == 0) {
		// Initial matrix exists only on the pivot process
		pMatrix = new int[SizeR * SizeC];
		// Values of elements are defined only on the pivot process
		DataInitialization(pMatrix, V); // , E, SizeR, SizeC);
	}
}

// Function for formatted matrix output
void PrintMatrix(int* pMatrix, int RowCount, int ColCount) {
	int i, j; // Loop variables
	for (i = 0; i < RowCount; i++) {
		for (j = 0; j < ColCount; j++)
			printf("%d ", pMatrix[i * ColCount + j]);
		printf("\n");
	}
}


// Function for distribution of the initial objects between the processes
void DataDistribution(int* pMatrix, int* pProcRows, int SizeR, int SizeC, int RowNum) {
	int* pSendNum; // The number of elements sent to the process
	int* pSendInd; // The index of the first data element sent to the process
	int RestRows = SizeR; // Number of rows, that haven’t been distributed yet

	// Alloc memory for temporary objects
	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];

	// Define the disposition of the matrix rows for current process
	RowNum = (SizeR / ProcNum);
	pSendNum[0] = RowNum * SizeC;
	pSendInd[0] = 0;
	
	for (int i = 1; i < ProcNum; i++) 
	{
		RestRows -= RowNum;
		RowNum = RestRows / (ProcNum - i);
		pSendNum[i] = RowNum * SizeC;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
	}
	// Scatter the rows
	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT, pProcRows, pSendNum[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
	
	// Free the memory
	delete[] pSendNum;
	delete[] pSendInd;
}

void TestDistribution(int* pMatrix, int* pProcRows, int SizeR, int SizeC, int RowNum) {
	if (ProcRank == 0) {
		printf("Initial Matrix: \n");
		PrintMatrix(pMatrix, SizeR, SizeC);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < ProcNum; i++) {
		if (ProcRank == i) {
			printf("\nProcRank = %d \n", ProcRank);
			printf(" Matrix Stripe:\n");
			PrintMatrix(pProcRows, RowNum, SizeC);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

// Function for computational process termination
void ProcessTermination(int* pMatrix, int* pResult, int* pProcRows, int* pProcResult) {
	if (ProcRank == 0)
		delete[] pMatrix;
	delete[] pResult;
	delete[] pProcRows;
	delete[] pProcResult;
}


void ParallelBoruvkaMST(int* pMatrix, int V, int E, int* pProcRows, int RowNum)
{
	struct Set* sets = new Set[V];

	// масив для збереження індексів ребер найменшої ваги.
	int* cheapest = new int[V];
	int* pProcCheap;
	int CheapNum;

	int rest = V;
	for (int i = 0; i < ProcRank; i++)
		rest = rest - rest / (ProcNum - i);
	CheapNum = rest / (ProcNum - ProcRank);

	// Memory allocation
	pProcCheap = new int[CheapNum];
	//int* pProcCheapest = new int;
	//int CheapestNum;
	// Create V subsets with single elements
	for (int v = 0; v < V; v++)
	{
		sets[v].parent = v;
		sets[v].rank = 0;
		cheapest[v] = -1;
	}

	// Початково маємо V різних дерев.
	// Фінально отримаємо одне дерево, що і буде MST.
	int numTrees = V;
	double MSTweight = 0;

	MPI_Bcast(&numTrees, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&MSTweight, 1, MPI_INT, 0, MPI_COMM_WORLD);

	DataDistribution(pMatrix, pProcRows, E, 3, RowNum);
	//DataDistribution(cheapest, pProcCheap, V, 1, CheapNum);

	// Keep combining components (or sets) until all compnentes are not combined into single MST.
	while (numTrees > 1)
	{
		if (ProcRank == 0) 
			for (int v = 0; v < V; v++)
			{
				cheapest[v] = -1;
			}
		//MPI_Bcast(&cheapest[0], V, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < RowNum; i++)
		{
			int set1 = find(sets, pProcRows[3 * i]);
			int set2 = find(sets, pProcRows[3 * i + 1]);

			// If two corners of current edge belong to same set, ignore current edge. 
			if (set1 == set2)
				continue;
			// Else check if current edge is closer to previous cheapest edges of set1 and set2
			else
			{
				if (cheapest[set1] == -1 || pProcRows[3 * cheapest[set1] + 2] > pProcRows[3 * i + 2])
					cheapest[set1] = i;
				if (cheapest[set2] == -1 || pProcRows[3 * cheapest[set2] + 2] > pProcRows[3 * i + 2])
					cheapest[set2] = i;
			}
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		
		// Consider the above picked cheapest edges and add them to MST
		if (ProcRank == 0)
		{
			for (int v = 0; v < V; v++)
			{
				// Check if cheapest for current set exists
				if (cheapest[v] != -1)
				{
					int set1 = find(sets, pProcRows[3 * cheapest[v]]);
					int set2 = find(sets, pProcRows[3 * cheapest[v] + 1]);

					if (set1 == set2)
						continue;
					MSTweight += pProcRows[3 * cheapest[v] + 2];

					Union(sets, set1, set2);
					
					numTrees--;
				}
			}
			for (int i = 1; i < ProcNum; i++)
				MPI_Send(&numTrees, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Status status;
		if (ProcRank != 0)
		{
			MPI_Recv(&numTrees, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		}
	}
	if (ProcRank == 0)
	{
		printf("Weight of MST is %f\n", MSTweight);
	}
}


// A utility function to find set of an element i (uses path compression technique)
int find(struct Set sets[], int i)
{
	// find root and make root as parent of i
	// (path compression)
	if (sets[i].parent != i)
		sets[i].parent = find(sets, sets[i].parent);
	return sets[i].parent;
}

void Union(struct Set sets[], int x, int y)
{
	int xroot = find(sets, x);
	int yroot = find(sets, y);

	// Attach smaller rank tree under root of high
	// rank tree (Union by Rank)
	if (sets[xroot].rank < sets[yroot].rank)
		sets[xroot].parent = yroot;
	else if (sets[xroot].rank > sets[yroot].rank)
		sets[yroot].parent = xroot;
	// If ranks are same, then make one as root and
	// increment its rank by one
	else
	{
		sets[yroot].parent = xroot;
		sets[xroot].rank++;
	}
}


int main(int argc, char* argv[])
{
	int* pMatrix;  // The first argument - initial matrix
	int* pResult;  // Result vector for matrix-vector multiplication 
	int V, E;
	int SizeR, SizeC;		    // Sizes of initial matrix and vector
	int* pProcRows;   // Stripe of the matrix on the current process
	int* pProcResult; // Block of the result vector on the current process
	int RowNum;          // Number of rows in the matrix stripe
	double Start, Finish, Duration;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		printf("Parallel Minimum Spanning Tree Program\n\n");
	}

	// Memory allocation and data initialization
	ProcessInitialization(pMatrix, pResult, pProcRows, pProcResult, SizeR, SizeC, V, E, RowNum);

	Start = MPI_Wtime();

	ParallelBoruvkaMST(pMatrix, V, E, pProcRows, RowNum);

	Finish = MPI_Wtime();

	Duration = Finish - Start;

	if (ProcRank == 0)
	{
		printf("\n\nTime of execution = %f\n", Duration);
	}
	// Process termination
	ProcessTermination(pMatrix, pResult, pProcRows, pProcResult);

	MPI_Finalize();

	return 0;
}