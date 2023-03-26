#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "trie.h"
#include "kseq.h"
#include "ds.h"
#include "dnautil.h"

KSEQ_INIT(int, read)

void look_up_ngg(size_t len, char *s, char* contig, Trie *trie);

Trie* build_ngg_trie_from_source(char* source, char* out_file_name);