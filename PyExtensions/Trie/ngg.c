#include "include/ngg.h"
#include "include/dnautil.h"
#include "ahocorasick/ahocorasick.h"
#include "include/convert.h"

Trie *build_ngg_trie_from_source(char *source, char* out_file_name)
{
    initNtCompTable();
    kseq_t *seq;
    FILE* fp;
    int n = 0, slen = 0, qlen = 0;
    if (!(fp = fopen(source, "r")))
        return NULL;
    seq = kseq_init(fileno(fp));
    struct ahocorasick aho;
    int id[8] = {0};

    aho_init(&aho);

    id[0] = aho_add_match_text(&aho, "AGG", 3);
    id[1] = aho_add_match_text(&aho, "TGG", 3);
    id[2] = aho_add_match_text(&aho, "CGG", 3);
    id[3] = aho_add_match_text(&aho, "GGG", 3);
    id[4] = aho_add_match_text(&aho, "CCA", 3);
    id[5] = aho_add_match_text(&aho, "CCT", 3);
    id[6] = aho_add_match_text(&aho, "CCC", 3);
    id[7] = aho_add_match_text(&aho, "CCG", 3);

    
    aho_create_trie(&aho);
    Trie *trie = Trie_new();
    printf("Out file name: %s\n", out_file_name);
    FILE* wfp = fopen(out_file_name, "w+");
    while (kseq_read(seq) >= 0)
    {
        printf("Processing Chromosome %s\n", seq->name.s);
        upper_str(seq->seq.s, strlen(seq->seq.s));
        findNGG(&aho, trie, seq->name.s, seq->seq.s, wfp, &callback_match_pos);
        // look_up_ngg(seq->seq.l, seq->seq.s, seq->name.s, trie);
        ++n;
    }
    aho_destroy(&aho);
    kseq_destroy(seq);
    fclose(fp); // FIXME!
    return trie;
}