#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SEQNUM 1e4
#define KMERSIZE 5

typedef struct _Kmer {
	int count;
	char seq[KMERSIZE+1];
} Kmer;

char **matrix;
Kmer *kmers;

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
	int i;

	for(i = 0; i < KMERSIZE; i++){
		nuc[i] = 'A';
	}
	nuc[KMERSIZE] = '\0';
	i = 0;

	while(num != 0){
		m = num % 4;
		num = num / 4;
		switch(m){
			case 0:
				nuc[i] = 'A';
				break;
			case 1:
				nuc[i] = 'C';
				break;
			case 2:
				nuc[i] = 'T';
				break;
			case 3:
				nuc[i] = 'G';
		};
		i++;
	}
}

void fillkmers(){
	int i;

	for(i = 0; i < pow(4, KMERSIZE); i++){
		int2nuc(i, kmers[i].seq);
		kmers[i].count = 0;
	}
}

void countkmers(){
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;
	int r;
	int index;

	for(i = 0; i < SEQNUM; i++){
		for(j = 0; j < strlen(matrix[i]) - KMERSIZE + 1; j++){
			index = 0;
			for(k = j, l = 1; k < j + KMERSIZE; k++, l *= 4){
				r = nuc2int(matrix[i][k]);
				if(r < 0){
					break;
				}
				index += r * l;
			}
			if(r >= 0){
				kmers[index].count += 1;
			}
		}
	}
}

void swapkmers(Kmer *kmers, int first, int second){
	int swapcount;
	char swapseq[KMERSIZE+1];

	swapcount = kmers[first].count;
	strncpy(swapseq, kmers[first].seq, KMERSIZE);

	kmers[first].count = kmers[second].count;
	strncpy(kmers[first].seq, kmers[second].seq, KMERSIZE);

	kmers[second].count = swapcount;
	strncpy(kmers[second].seq, swapseq, KMERSIZE);
}

void quicksort(Kmer *kmers, int left, int right){
	int i = left;
	int j = right;
	int pivot = kmers[(left + right) / 2].count;

	while(i <= j){
		while(kmers[i].count < pivot) i++;
		while(kmers[j].count > pivot) j--;
		if(i <= j){
			swapkmers(kmers, i, j);
			i++;
			j--;
		}
	}

	if(left < j) quicksort(kmers, left, j);
	if(right > i) quicksort(kmers, i, right);
}

void buildadaptor(){
	int i, j, k;
	int *pos;
	int found;
	int fcount;
	char nuc;
	int dist[4];

	pos = malloc(sizeof(int) * SEQNUM);

	/* Find the position of the most frequent kmer in the sequence set */
	for(i = 0; i < SEQNUM; i++){
		fcount = 0;
		for(j = 0; j < strlen(matrix[i]) - KMERSIZE; j++){
			found = 1;
			for(k = 0; k < KMERSIZE; k++){
				if(matrix[i][j + k] != kmers[0].seq[k]){
					found = 0;
					break;
				}
			}
			if(found){
				pos[i] = j;
				fcount++;
			}
		}
		if(fcount != 1){
			pos[i] = -1;
		}
	}

	while(1){

		dist[0] = 0;
		dist[1] = 0;
		dist[2] = 0;
		dist[3] = 0;

		for(i = 0; i < SEQNUM; i++){
			if(pos[i] > 0){
				nuc = matrix[i][pos[i] - 1];
				switch(nuc){
					case 'A':
						dist[0]++;
						break;
					case 'T':
						dist[1]++;
						break;
					case 'G':
						dist[2]++;
						break;
					case 'C':
						dist[3]++;
				};
			}
		}

		if(decision(dist)){
			break;
		}
	}

	free(pos);
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

	kmers = malloc(pow(4, KMERSIZE) * sizeof(Kmer));
	if(kmers == NULL){
		fprintf(stderr, "Not enough memory for this kmer size\n");
		return(3);
	}
	fillkmers();
	countkmers();
	quicksort(kmers, 0, pow(4, KMERSIZE)-1);
	buildadaptor();
	for(r = 0; r < pow(4, KMERSIZE); r++){
		printf("%s %d\n", kmers[r].seq, kmers[r].count);
	}
	free(kmers);

	return(EXIT_SUCCESS);
}
