#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SEQNUM 1e4
#define KMERSIZE 5

char **matrix;
int *kmers;

int readfirstN(char *filename, int maxseq){
	FILE   *seqfile;
	char   *line       = NULL;
	size_t  linelen    = 0;
	ssize_t readed;
	int     seqcounter = 0;
	int     linecount  = 0;
	int     returncode = 0;
	int     mod;

	seqfile = fopen(filename, "r");
	if(seqfile == NULL){
		fprintf(stderr, "Fail to read sequence file\n");
		return(-1);
	}

	while( (readed = getline(&line, &linelen, seqfile)) != -1 ){
		/* Determiny the file type from the first line */
		if(linecount == 0){
			/* This is a fastq file */
			if(line[0] == '@'){
				mod = 4;
			}
			/* This is a fasta file */
			else if(line[0] == '>'){
				mod = 2;
			}
			/* I do not know what kind of file it is */
			else{
				returncode = -2;
				goto end;
			}
		}

		/* This is a sequence line */
		if(linecount % mod == 1){

			matrix[seqcounter] = malloc(sizeof(char) * readed);
			strncpy(matrix[seqcounter], line, readed - 1);
			matrix[seqcounter][readed-1] = '\0';

			if(matrix[seqcounter] == NULL){
				returncode = -3;
				goto end;
			}
			seqcounter++;
			if(seqcounter == maxseq){
				goto end;
			}
		}

		linecount++;
	}

end:
	fclose(seqfile);
	free(line);
	return(returncode);
}

int nuc2int(char nucl){
	switch(nucl){
		case 'A':
		case 'a':
			return(0);
		case 'C':
		case 'c':
			return(1);
		case 'T':
		case 't':
			return(2);
		case 'G':
		case 'g':
			return(3);
	};
	return(-1);
}

void int2nuc(int index, char *nuc){
	int num = index;
	int m;

	while(num != 0){
		m = num % 4;
		num = num / 4;
		switch(m){
			case 0:
				printf("A");
				break;
			case 1:
				printf("C");
				break;
			case 2:
				printf("T");
				break;
			case 3:
				printf("G");
		};
	}
	printf(" ");
}

void calckmers(){
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;
	int r;
	int index;

	for(i = 0; i < SEQNUM; i++){
		for(j = 0; j < strlen(matrix[i]) - KMERSIZE; j++){
			index = 0;
			for(k = j, l = 1; k < j + KMERSIZE; k++, l *= 4){
				r = nuc2int(matrix[i][k]);
				if(r < 0){
					break;
				}
				index += r * l;
			}
			if(r >= 0){
				kmers[index] += 1;
			}
		}
	}
}

int main(int argc, char **argv){
	int r;

	if(argc < 2){
		fprintf(stderr, "Please specify the file name\n");
		return(1);
	}

	matrix = malloc(sizeof(char*) * SEQNUM);
	if(matrix == NULL){
		return(3);
	}

	r = readfirstN(argv[1], SEQNUM);
	if(r < 0){
		fprintf(stderr, "Something wrong happened during file processing\n");
		return(2);
	}

	kmers = calloc(pow(4, KMERSIZE), sizeof(int));
	if(kmers == NULL){
		fprintf(stderr, "Not enough memory for this kmer size\n");
		return(3);
	}
	calckmers();

	free(kmers);

	return(EXIT_SUCCESS);
}
