#include <string.h>
#include "usfdynarr.h"
#include "usfio.h"
#include <omp.h>
#include <time.h>
#include <unistd.h>

/* Last defined char is END, 27th is WILDCARD */
#define MAXCHARS 28
#define MAXWORD 20
#define UPDATE 15000

char **rawformat, **wordlist, *word;
uint64_t formatlen, words, wlen, blen, solutions;

FILE *outstream;
unsigned char **wtree, **format;

void fit(unsigned, unsigned char *);
void start(char *);

int main(int args, char *argv[]) {
	clock_t time;
	uint64_t i, j, k;
	int c;
	unsigned char **wnode, skip[MAXWORD] = {0};
	usf_dynarr *matchingwords;

	/* Syntax: topwords [FORMAT] [WORDLIST] */

	if (args < 3) {
		printf("Syntax: topwords [FORMAT] [WORDLIST] ?[OUTFILE]\n");
		exit(1);
	}

	if (args > 3) outstream = fopen(argv[3], "w");
	else outstream = stdout;

	printf("Starting TopologicalCrosswords V2.0 with format %s and wordlist %s.\n", argv[1], argv[2]);
	printf("NOTE: Wordlists must be in lowercase and must only contain ASCII a-z characters !\n");

	/* Read data */
	rawformat = usf_ftot(argv[1], "r", &formatlen);
	wordlist = usf_ftot(argv[2], "r", &words);

	blen = strlen(rawformat[0]) - 1; //Number of slots in board

	if (blen >= 256) {
		printf("TopologicalCrosswords: Error: format declaration exceeds byte limit.\n");
		exit(1);
	}

	printf("Building format...\n");

	/* Build absolute format */
	format = malloc(sizeof(unsigned char *) * (formatlen - 1));

	for (i = 1; i < formatlen; i++) {
		/* Len of one format line */
		k = strlen(rawformat[i]) - 1;
		skip[k] = 1;

		if (k >= 255) {
			printf("TopologicalCrosswords: ERROR: format exceeds byte limit.\n");
			exit(1);
		}

		format[i - 1] = malloc(k + 1); //Include ending char

		/* Transform format to offset in board */
		for (j = 0; j < k; j++)
			format[i - 1][j] = strchr(rawformat[0], rawformat[i][j]) - rawformat[0];

		format[i - 1][j] = 255;
	}

	printf("Building word tree...\n");

	/* Using the trie structure:
	 * each successive node is an array with
	 * 28 entries, representing chars a-z plus one wildcard and
	 * and leafnode. Using any other character is not supported,
	 * although this is doable provided the char is not { (wildcard)
	 */

	time = clock();
	wtree = calloc(sizeof(unsigned char *), MAXCHARS);

	for (i = 0; i < words; i++) {
		if (i % 5000 == 0)
			printf("%ld / %ld\n", i, words);

		/* For each word :
		 * Mask with all binary numbers from
		 * 0 to 2^n where n is len of word
		 * insert word into the usf_dynarr present
		 * when traversing the tree using the word's chars
		 * and substituting wildcards where there are binary 0s
		 */

		word = wordlist[i];
		wlen = strlen(word);

		word[--wlen] = '\0'; //Remove \n

		if (wlen >= MAXWORD || !skip[wlen]) continue;

		for (j = 0; j < 1UL << wlen; j++) {
			wnode = wtree;

			for (k = 0; k < wlen; k++) {
				/* Purple magic */
				c = (j & (1 << k) ? word[k] : '\173') - 'a';

				/* New subtree if not present */
				if (wnode[c] == NULL)
					wnode[c] = calloc(sizeof(unsigned char *), MAXCHARS);

				wnode = (unsigned char **) wnode[c];
			}

			matchingwords = (usf_dynarr *) wnode[MAXCHARS - 1];

			if (matchingwords == NULL)
				wnode[MAXCHARS - 1] = (unsigned char *) (matchingwords = usf_arrtodyn(NULL, 0));

			/* Append word */
			usf_daset(matchingwords, matchingwords -> size, USFDATAP(word));
		}
	}

	printf("Done ! (%f seconds)\n", ((double) clock() - time) / CLOCKS_PER_SEC);

	/* Visual confirmation for 0.2 seconds */
	usleep(200000);

	//Omit \n
	wlen = strlen(rawformat[1]) - 1;

	time = clock();
	j = 0;

	#pragma omp parallel for
	for (i = 0; i < words; i++) {
		j++;
		if (j % 100 == 0)
			printf("%lu / %lu\n", j, words);

		word = wordlist[i];
		if (strlen(word) != wlen)
			continue;

		start(word);
	}

	printf("Finished after %f CPU seconds, found %lu solutions.\n",
			((double) clock() - time) / CLOCKS_PER_SEC, solutions);

	return 0;
}

void start(char *word) {
	unsigned char *board, *form;
	uint64_t i, j;

	board = malloc(blen);

	/* Set all wildcard */
	for (i = 0; i < blen; i++)
		board[i] = 255;

	/* Set first word */
	for (j = 0, form = format[0]; *form != 255; form++, j++)
		board[*form] = word[j];

	fit(1, board);
	free(board);
}

void fit(unsigned step, unsigned char *board) {
	unsigned char c, *form, **wnode, *b;
	uint64_t i, j;
	usf_data *arr;
	usf_dynarr *matching;

	/* Get this step's format */
	form = format[step];

	/* Root note */
	wnode = wtree;

	/* Lookup matching words */
	for (; *form != 255; form++) {
		c = board[*form] - 'a';

		if (c >= MAXCHARS)
			/* Wildcard */
			wnode = (unsigned char **) wnode[26];
		else wnode = (unsigned char **) wnode[c];

		if (wnode == NULL)
			/* Invalid pattern; prune */
			return;
	}

	if ((matching = (usf_dynarr *) wnode[MAXCHARS - 1]) == NULL)
		/* No words */
		return;

	arr = matching -> array;

	/* Create new board and set all matching */
	b = malloc(blen);
	memcpy(b, board, blen);

	for (i = 0; i < matching -> size; i++) {
		/* Copy matching word */
		for (j = 0, form = format[step]; *form != 255; form++, j++)
			b[*form] = ((char *) arr[i].p)[j];

		/* Win condition */
		if (step == formatlen - 2) {
			fwrite(b, 1, blen, outstream);
			fputc('\n', outstream);
			solutions++;
		} else
			fit(step + 1, b);
	}

	/* Destroy intermediate board */
	free(b);
}
