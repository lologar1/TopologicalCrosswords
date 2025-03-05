#include "usfio.h"
#include "usfstring.h"
#include <time.h>
#include <omp.h>

#define WILDCARD '*'

_Atomic uint64_t solutions;
uint64_t formatlen;
int **format, **query[256], declen;
char **words;
usf_dynarr *wtree[256];
FILE *outstream;

void start(char *);
void fit(int, char *);

int main(int args, char *argv[]) {
	uint64_t i, j, n, k;
	uint64_t rawlen, wordlen;
	char **rawformat, *declaration, *instruction, *word;
	char possibleLengths[256] = {0}; /* To avoid constructing tree with unused words */
	usf_dynarr *matches, **treenode;
	struct timespec begin, end;

	if (args < 3) {
		fprintf(stderr, "Syntax: topwords [FORMAT] [WORDLIST] ?[OUTFILE]\n");
		exit(1);
	}

	fprintf(stderr, "Starting TopologicalCrosswords V3.0 with format %s and wordlist %s.\n", argv[1], argv[2]);
	fprintf(stderr, "IMPORTANT: wordlists must end with linefeed *only* and be in lowercase !\n");
	fprintf(stderr, "IMPORTANT: wordlists and formats must only contain ASCII a-z chars, except for the declaration which is only limited to ASCII characters.\n");

	/* Initialization */

	/* Read format file */
	rawformat = usf_ftot(argv[1], "r", &rawlen);

	/* Len of format without declaration */
	formatlen = rawlen - 1;

	/* Read wordlist */
	words = usf_ftot(argv[2], "r", &wordlen);

	/* Output stream is stdout if not specified */
	outstream = args > 3 ? fopen(argv[3], "w") : stdout;

	/* Error handling */

	if (rawformat == NULL) {
		perror("Couldn't open format file !\n");
		exit(2);
	}

	if (words == NULL) {
		perror("Couldn't open word list !\n");
		exit(2);
	}

	if (rawlen < 2) {
		perror("Unintelligible format !\n");
		exit(2);
	}

	printf("Constructing format... ");

	/* Remove all \n endings from word list and format */
	for (i = 0; i < wordlen; i++) {
		declaration = words[i];
		declaration[strlen(declaration) - 1] = '\0';
	}

	for (i = 0; i < rawlen; i++) {
		declaration = rawformat[i];
		declaration[strlen(declaration) - 1] = '\0';
	}

	/* First line declares all chars to be used */
	declaration = rawformat[0];
	declen = strlen(declaration);

	if (declen >= 256) {
		perror("Format declaration too long !\n");
		exit(2);
	}

	/* Convert step -> affected slots, ending with -1*/
	format = malloc(sizeof(int *) * formatlen);

	for (i = 0; i < formatlen; i++) {
		/* +1 to omit declaration */
		instruction = rawformat[i + 1];

		/* Length of this format instruction */
		n = strlen(instruction);

		/* Make words of this length mandatory */
		possibleLengths[n] = 1;

		/* +1 to account for terminating -1 */
		format[i] = malloc(sizeof(int) * (n + 1));

		/* Index of char in declaration is its slot number */
		for (j = 0; j < strlen(instruction); j++)
			format[i][j] = strchr(declaration, instruction[j]) - declaration;

		/* Terminating -1 */
		format[i][j] = -1;
	}

	/* DEBUG
	for (i = 0; i < formatlen; i++) {
		printf("For format %lu: ", i);

		for (j = 0; j < strlen(rawformat[i + 1]); j++) {
			printf("%d ", format[i][j]);
		}
		printf("\n");
	} */

	/* Convert affected slot -> array of words to be fitted
	 * Array of pointers is terminated by 0 (NULL pointer)
	 * Words end with -1 as per format */

	for (i = n = 0; i < (unsigned) declen; i++) {
		/* Overallocating; no impact on performance */
		query[i] = malloc(sizeof(int *) * (formatlen + 1));

		/* For every instruction, if it features the
		 * character at index i in the declaration,
		 * then add it to this query */

		for (j = 0; j < formatlen; j++)
			if (strchr(rawformat[j + 1], declaration[i]))
				query[i][n++] = format[j];

		/* Terminating 0 */
		query[i][n] = NULL;
		n = 0;
	}

	/* DEBUG
	for (i = 0; i < strlen(declaration); i++) {
		printf("For char %c\n", declaration[i]);

		for (j = 0; query[i][j]; j++) {
			for (n = 0; query[i][j][n] != -1; n++)
				printf("%d ", query[i][j][n]);
			printf("\n");
		}
	} */

	printf("Done !\n");

	printf("Constructing word tree... ");

	/* For each word, mask with each binary number
	 * from 0 to 2^(word length) to include wildcards
	 * and append to list */

	/* Playground for masking */
	instruction = malloc(sizeof(char) * 256);

	for (i = 0; i < wordlen; i++) {
		/* Root word */
		word = words[i];

		/* Length */
		n = strlen(word);

		/* Skip if possible */
		if (n >= 256 || !possibleLengths[n]) continue;

		/* Begin wildcard masking */
		for (j = 0; j < (1UL << n); j++) {
			strcpy(instruction, word);

			for (k = 0; k < n; k++)
				if (j & (1UL << k))
					instruction[k] = WILDCARD;

			/* Start at top */
			treenode = wtree;

			/* Build trie, use unsigned for array indexing */
			for (k = 0; k < n; k++) {
				if (treenode[(unsigned char) instruction[k]]) {
					treenode = (usf_dynarr **) treenode[(unsigned char) instruction[k]];
				} else {
					treenode = (usf_dynarr **) (treenode[(unsigned char) instruction[k]] = calloc(sizeof(usf_dynarr **), 256));
				}
			}

			/* First element is word ending */
			matches = (usf_dynarr *) treenode[0];

			if (matches == NULL)
				treenode[0] = (matches = usf_newda(0));

			declaration = malloc(sizeof(char) * (n + 1));
			strcpy(declaration, word);

			usf_daappend(matches, USFDATAP(declaration));
		}
	}

	free(instruction);

	/* DEBUG
	//Let's seach everything starting with air, 5 long !
	declaration = "air*";
	treenode = wtree;

	for (i = 0; i < strlen(declaration); i++) {
		if (treenode == NULL) {
			fprintf(stderr, "No such words\n");
			exit(1);
		}

		treenode = (usf_dynarr **) treenode[declaration[i]];
	}

	matches = treenode[0];
	if (matches == NULL) {
		fprintf(stderr, "No words like this\n");
		exit(1);
	}

	fprintf(stderr, "%lu words\n", matches -> size);

	for (i = 0; i < matches -> size; i++) {
		fprintf(stderr, "%s\n", matches -> array[i]);
	} */

	printf("Done !\n");

	printf("Starting process...\n");

	clock_gettime(CLOCK_REALTIME, &begin);

	/* Starting word length */
	n = strlen(rawformat[1]);

	/* Counter */
	j = 0;

#pragma omp parallel for
	for (i = 0; i < wordlen; i++) {
		/* Progress tracking */
		if (j % 1024 == 0)
			fprintf(stderr, "Completed %lu / %lu\n", j, wordlen);
		j++;

		/* Base word */
		word = words[i];

		/* Cannot fit */
		if (strlen(word) != n)
			continue;

		/* Start a process */
		start(word);
	}

	clock_gettime(CLOCK_REALTIME, &end);

	fprintf(stderr, "Finished after %f seconds, found %lu solutions.\n", (double) (end.tv_sec - begin.tv_sec) + (double) (end.tv_nsec - begin.tv_nsec) / 1000000000.0, solutions);

	/* Cleanup */

	/* Destroy format */
	for (i = 0; i < formatlen; i++)
		free(format[i]);
	free(format);

	/* Destroy query */
	for (i = 0; i < (unsigned) declen; i++)
		free(query[i]);

	/* Destroy word list */
	for (i = 0; i < wordlen; i++)
		free(words[i]);
	free(words);

	/* Destroy raw format */
	for (i = 0; i < rawlen; i++)
		free(rawformat[i]);
	free(rawformat);

	/* Close output stream when finished if need be */
	if (args > 3) fclose(outstream);

	return 0;
}

void start(char *baseword) {
	int i, *instruction;
	char *board;

	/* Allocate board +1 for printable \n char */
	board = malloc(sizeof(char) * (declen + 1));

	/* Fill with wildcards */
	for (i = 0; i < declen; i++)
		board[i] = WILDCARD;
	board[i] = '\n';

	for (i = 0, instruction = format[0]; *instruction != -1; i++, instruction++)
		board[*instruction] = baseword[i];

	/* Start digging */
	fit(1, board);

	/* Cleanup */
	free(board);
}

void fit(int step, char *board) {
	int modified[256], *form, *i, n, e;
	int **answer, *instruction;
	char **nominees, *attempt;
	unsigned char current;
	usf_dynarr *branches, **treenode;

	/* Instruction for this step */
	form = format[step];

	/* Keep track of modified characters */
	e = 0;

	/* Traverse trie */
	treenode = wtree;

	for (i = form; *i != -1; i++) {
		/* Current slot being accessed */
		n = *i;
		current = (unsigned char) board[n];

		/* Get the word defined by this format into candidate */
		if (current == WILDCARD)
			/* Will have to be reset later */
			modified[e++] = n;

		/* Cul-de-sac */
		if ((treenode = (usf_dynarr **) treenode[current]) == NULL) return;
	}

	/* Terminating -1 to stop resetting */
	modified[e] = -1;

	/* All valid branches */
	if ((branches = (usf_dynarr *) treenode[0]) == NULL)
		return; /* Also cul-de-sac */

	/* All possible next words */
	nominees = (char **) branches -> array;

	/* Test all nominees */
	for (n = 0; (unsigned) n < branches -> size; n++) {
		/* Current nominee */
		attempt = nominees[n];

		for (e = 0, i = form; *i != -1; i++, e++) {
			/* Place one character */
			board[*i] = attempt[e];

			/* Don't bother with the last one */
			if ((unsigned) step == formatlen - 1)
				continue;

			/* For every overlapping instruction (format) */
			for (answer = query[*i]; *answer; answer++) {
				/* Prepare to traverse trie */
				treenode = wtree;

				/* Check existence */
				for (instruction = *answer; *instruction != -1; instruction++) {
					if ((treenode = (usf_dynarr **) treenode[(unsigned char) board[*instruction]]) == NULL)
						goto cull;
				}

				if (treenode[0] == NULL) goto cull;
			}
		}

		/* Valid solution */
		if ((unsigned) step == formatlen - 1) {
			/* Include \n char to print board */
			fwrite(board, sizeof(char), declen + 1, outstream);
			solutions++;
		} else {/* Recurse one branch down */
			fit(step + 1, board);
		}
cull:
		/* Reset to avoid integrity malfunction */
		for (i = modified; *i != -1; i++)
			board[*i] = WILDCARD;
	}
}
