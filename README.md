# TopologicalCrosswords
Program for bruteforcing custom crossword layouts

Important: usflib2 is required for compilation (available at https://github.com/lologar1/usflib2)

To use: simply execute and specify, in order, the FORMAT file, a WORDLIST, and an OUTPUT FILE (defaults to stdout)

WORDLISTs must contain one word per line, ending with a line feed (\n) only.

FORMAT files are specified as such:
The first line is a declaration of all tokens to be used in the arrangement. For example, a 4x4 word square would have 16 tokens (16 characters). Note that only ASCII chars are supported as tokens at the moment, limiting them to 255, although that is something I'm working on.
All subsequent lines declare a valid word (from the WORDLIST) to be formed by specific tokens. For instance, a 3x3 word square specification file could look like this, without the parenthesized stuff of course:

abcdefghi (nine tokens, for the nine characters)
abc (first row)
adg (first column)
def (second row)
beh (second column)
ghi (third row)
cfi (third column)

The program will try to fit the words one after the other, and does so "smartly" (will cull roots with no possible branches), but be aware that placing all rows then all columns will considerably slow the process.
Again, the behavior and algorithm are also something I'm working to improve at the moment.

Thanks, lologar1
