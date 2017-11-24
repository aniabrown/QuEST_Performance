#include <stdio.h>
#include <stdlib.h>

#define N_QUBITS 28
#define TRIAL 1

int main(int argc, char **argv){
	FILE *fp;
	char filename[100];
	int rand1, rand2, rand3, randbool;
	int i;
	sprintf(filename, "q%d/c%d.txt", N_QUBITS, TRIAL);
	fp = fopen(filename, "w");

	for (i=0; i<150; i++){
		rand1 = (rand() % N_QUBITS) + 1;
		rand2 = (rand() % N_QUBITS) + 1;
		while(rand2==rand1){
			rand2 = (rand() % N_QUBITS) + 1;
		}
		rand3 = rand() % N_QUBITS + 1;
		while(rand3==rand2 || rand3==rand1){
			rand3 = rand() % N_QUBITS + 1;
		}
		randbool = rand() % 2;
		if (randbool == 1) rand1 *= -1;
		randbool = rand() % 2;
		if (randbool == 1) rand2 *= -1;
		randbool = rand() % 2;
		if (randbool == 1) rand3 *= -1;
		fprintf(fp, "%d %d %d\n", rand1, rand2, rand3);
	}

	fclose(fp);
	return 1;
}
