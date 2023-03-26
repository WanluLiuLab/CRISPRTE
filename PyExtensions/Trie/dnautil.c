#include "include/dnautil.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * the struct Gloc behaves like a linked list
 */
struct Gloc
{
    char *pam;
    char *contig;
    char *up;
    char *dw;
    size_t start;
    size_t end;
    char strand;
    Gloc *next;
};

Gloc *Gloc_new_args(char *pam_, char *contig_, char *up_, char *dw_, size_t st, size_t en, char sd)
{
    Gloc *new_gloc = (Gloc *)malloc(sizeof(Gloc));
    Gloc_set(new_gloc, pam_, contig_, up_, dw_, st, en, sd);
    return new_gloc;
}

Gloc *Gloc_new(void)
{
    return (Gloc *)malloc(sizeof(Gloc));
}

void Gloc_del(Gloc *gloc)
{
    free(gloc->contig);
    free(gloc->dw);
    free(gloc->pam);
    free(gloc->up);
    free(gloc);
}

void Gloc_set(Gloc *gloc, char *pam_, char *contig_, char *up_, char *dw_, size_t st, size_t en, char sd)
{
    gloc->pam = pam_;
    gloc->contig = contig_;
    gloc->up = up_;
    gloc->dw = dw_;
    gloc->start = st;
    gloc->end = en;
    gloc->strand = sd;
    gloc->next = (void *)NULL;
}

void Gloc_map_start(Gloc *gloc, void *(*func)(void *))
{
    gloc->start = func(gloc->start);
}

void Gloc_map_end(Gloc *gloc, void *(*func)(void *))
{
    gloc->end = func(gloc->end);
}

void Gloc_map_strand(Gloc *gloc, void *(*func)(void *))
{
    gloc->strand = func(gloc->strand);
}

void Gloc_map_pam(Gloc *gloc, void *(*func)(void *))
{
    char *old = gloc->pam;
    gloc->pam = func(gloc->pam);
    free(old);
}

void Gloc_map_contig(Gloc *gloc, void *(*func)(void *))
{
    char *old = gloc->contig;
    gloc->contig = func(gloc->contig);
    free(old);
}

void Gloc_map_up(Gloc *gloc, void *(*func)(void *))
{
    char *old = gloc->up;
    gloc->up = func(gloc->up);
    free(old);
}

void Gloc_map_dw(Gloc *gloc, void *(*func)(void *))
{
    char *old = gloc->dw;
    gloc->dw = func(gloc->dw);
    free(old);
}

size_t Gloc_get_start(Gloc *gloc)
{
    return gloc->start;
}

size_t Gloc_get_end(Gloc *gloc)
{
    return gloc->end;
}

char *Gloc_get_pam(Gloc *gloc)
{
    return gloc->pam;
}

char *Gloc_get_contig(Gloc *gloc)
{
    return gloc->contig;
}

char *Gloc_get_up(Gloc *gloc)
{
    return gloc->up;
}

char *Gloc_get_dw(Gloc *gloc)
{
    return gloc->dw;
}

char Gloc_get_strand(Gloc *gloc)
{
    return gloc->strand;
}

Gloc *Gloc_get_next(Gloc *gloc)
{
    return gloc->next;
}

void Gloc_set_next(Gloc *gloc, Gloc *value)
{
    Gloc *cur = gloc;
    while (cur->next)
    {
        cur = cur->next;
    }
    cur->next = value;
}

int Gloc_length(Gloc *gloc)
{
    int i = 1;
    Gloc *cur = gloc;
    while (cur->next)
    {
        cur = cur->next;
        i += 1;
    }
    return i;
}

/* Another array to help us do complement of char */
static char ntCompTable[256];
static boolean inittedCompTable = false;

void initNtCompTable(void)
{
    ntCompTable[' '] = ' ';
    ntCompTable['-'] = '-';
    ntCompTable['='] = '=';
    ntCompTable['a'] = 't';
    ntCompTable['c'] = 'g';
    ntCompTable['g'] = 'c';
    ntCompTable['t'] = 'a';
    ntCompTable['u'] = 'a';
    ntCompTable['n'] = 'n';
    ntCompTable['-'] = '-';
    ntCompTable['.'] = '.';
    ntCompTable['A'] = 'T';
    ntCompTable['C'] = 'G';
    ntCompTable['G'] = 'C';
    ntCompTable['T'] = 'A';
    ntCompTable['U'] = 'A';
    ntCompTable['N'] = 'N';
    ntCompTable['R'] = 'Y';
    ntCompTable['Y'] = 'R';
    ntCompTable['M'] = 'K';
    ntCompTable['K'] = 'M';
    ntCompTable['S'] = 'S';
    ntCompTable['W'] = 'W';
    ntCompTable['V'] = 'B';
    ntCompTable['H'] = 'D';
    ntCompTable['D'] = 'H';
    ntCompTable['B'] = 'V';
    ntCompTable['X'] = 'N';
    ntCompTable['r'] = 'y';
    ntCompTable['y'] = 'r';
    ntCompTable['s'] = 's';
    ntCompTable['w'] = 'w';
    ntCompTable['m'] = 'k';
    ntCompTable['k'] = 'm';
    ntCompTable['v'] = 'b';
    ntCompTable['h'] = 'd';
    ntCompTable['d'] = 'h';
    ntCompTable['b'] = 'v';
    ntCompTable['x'] = 'n';
    ntCompTable['('] = ')';
    ntCompTable[')'] = '(';
    inittedCompTable = true;
}

static char ntMixedCaseChars[256];

void initNtMixedCaseChars(void)
{
    static boolean initted = false;

    if (!initted)
    {
        ntMixedCaseChars['a'] = 'a';
        ntMixedCaseChars['A'] = 'A';
        ntMixedCaseChars['c'] = 'c';
        ntMixedCaseChars['C'] = 'C';
        ntMixedCaseChars['g'] = 'g';
        ntMixedCaseChars['G'] = 'G';
        ntMixedCaseChars['t'] = 't';
        ntMixedCaseChars['T'] = 'T';
        ntMixedCaseChars['n'] = 'n';
        ntMixedCaseChars['N'] = 'N';
        ntMixedCaseChars['u'] = 'u';
        ntMixedCaseChars['U'] = 'U';
        ntMixedCaseChars['-'] = 'n';
        initted = true;
    }
}

/* Reverse the order of the bytes. */
void reverseBytes(char *bytes, long length)
{
    long halfLen = (length >> 1);
    char *end = bytes + length;
    char c;
    while (--halfLen >= 0)
    {
        c = *bytes;
        *bytes++ = *--end;
        *end = c;
    }
}

char *complement(char *dna, long length)
{
    int i;
    char *buf = (char *)malloc(sizeof(char) * length + 1);
    memset(buf, 0, length + 1);
    for (i = 0; i < length; i++)
    {
        buf[i] = ntCompTable[(int)dna[i]];
    }
    return buf;
}

char *reverseComplement(char *dna, long length)
{
    reverseBytes(dna, length);
    char *cpld = complement(dna, length);
    free(dna);
    return cpld;
}

/* A little array to help us decide if a character is a 
 * nucleotide, and if so convert it to lower case. */
char ntChars[256];
static boolean initted = false;
static void initNtChars()
{

    if (!initted)
    {
        memset(ntChars, 0, sizeof(ntChars));
        ntChars['a'] = ntChars['A'] = 'a';
        ntChars['c'] = ntChars['C'] = 'c';
        ntChars['g'] = ntChars['G'] = 'g';
        ntChars['t'] = ntChars['T'] = 't';
        ntChars['u'] = ntChars['U'] = 'u';
        initted = true;
    }
}

boolean isValidDna(char *poly, int size)
{
    if (!initted)
        initNtChars();
    int i;
    int dnaCount = 0;

    for (i = 0; i < size; ++i)
    {
        if (!ntChars[(int)poly[i]])
            return false;
    }
    return true;
}

boolean isValidAmbiguityCodes(char *poly, int size)
{
    if (!initted)
        initNtChars();
    int i;
    int dnaCount = 0;

    for (i = 0; i < size; ++i)
    {
        if (!ntCompTable[(int)poly[i]])
            return false;
    }
    return true;
}

struct tmp
{
    Trie *trie;
    char *seq;
    char *contig;
    FILE* fp;
    size_t len;
};

int strngcmp(const char *__s1, const char *__s2, size_t __n)
{
    printf("%s,%s!!!!\n", __s1, __s2);
    unsigned char c1;
    unsigned char c2;
    for (size_t i = 0; i < __n; i++)
    {
        c1 = (unsigned char)*(__s1 + i);
        c2 = (unsigned char)*(__s2 + i);
        if (c1 == '*')
            continue;
        if (c1 != '\0' || c1 != c2)
            return c1 - c2;
    }
    return 0;
}

void upper_str(char *__s, size_t __size)
{
    for (size_t i = 0; i < __size; i++)
    {
        __s[i] = toupper(__s[i]);
    }
}

void slice_str(const char *__s, char *buf, size_t start, size_t end)
{
    size_t j = 0;
    for (size_t i = start; i <= end; ++i)
    {
        buf[j++] = __s[i];
    }
    buf[j] = 0;
}

char *slice_str_cpy(const char *__s, size_t start, size_t end)
{
    char *buffer = (char *)malloc(sizeof(char) * (end - start + 2));
    memset(buffer, 0, end - start + 2);
    size_t j = 0;
    for (size_t i = start; i <= end; i++)
    {
        buffer[j] = __s[i];
        j++;
    }
    return buffer;
}

char *cat_str(char *__s1, char *__s2)
{
    char *result = malloc(strlen(__s1) + strlen(__s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, __s1);
    strcat(result, __s2);
    return result;
}

char *nggcpy(char *__s)
{
    char *pam_ = (char *)malloc(4);
    memset(pam_, 0, 4);
    memcpy(pam_, __s, 3);
    upper_str(pam_, 3);
    return pam_;
}

char *strcpyn(char *__s, size_t __sz)
{
    char *new = (char *)calloc(0, __sz + 1);
    memcpy(new, __s, __sz);
    return new;
}

int _isValidGRNA_helper(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return 4;
    }
}

boolean isValidGRNA(char *seq, int thres)
{
    int count[4] = {0, 0, 0, 0};
    for (int i = 0; i < strlen(seq); ++i)
    {
        count[_isValidGRNA_helper(seq[i])]++;
    }
    for (int i = 0; i < 4; ++i)
    {

        if (count[i] >= thres)
            return false;
    }
    return true;
}

void callback_match_pos(void *arg, struct aho_match_t *m)
{
    
    size_t pos = m->pos;
    struct tmp *t = (struct tmp *)arg;

    char *s = strcpyn(t->seq + m->pos, 3);
    free(s);


    // if (t->seq[pos] == 'G' && pos > GRNA_L + 8 && pos < t->len - 9)
    if (m->id < 4 && pos > GRNA_L + 7 && pos < t->len - 10)
    {
        char *grna = slice_str_cpy(t->seq, pos - GRNA_L, pos - 1);

        if (isValidDna(grna, strlen(grna)) && isValidGRNA(grna, 17))
        {
            char *s = nggcpy(&t->seq[pos]);
            free(s);
            Gloc *g = Gloc_new_args(nggcpy(&t->seq[pos]), strcpyn(t->contig, strlen(t->contig)), slice_str_cpy(t->seq, pos - GRNA_L - 6, pos - GRNA_L - 1), slice_str_cpy(t->seq, pos + 3, pos + 8), pos - GRNA_L, pos, '+');
            Trie_set(t->trie, grna, g, t->fp, write_callback);
            Gloc_del(g);
        }
        else
            free(grna);
    }
    // if (t->seq[pos] == 'C' && pos > 7 && pos < t->len - GRNA_L - 9) {
    if (m->id >= 4 && pos > 6 && pos < t->len - GRNA_L - 10) {
        char *rcd = reverseComplement(slice_str_cpy(t->seq, pos + 3, pos + 2 + GRNA_L), GRNA_L);
        if (isValidDna(rcd, strlen(rcd)) && isValidGRNA(rcd, 17))
        {
            char *s = nggcpy(&t->seq[pos]);
            char *ngg = reverseComplement(s, 3);
           
            char *up = reverseComplement(slice_str_cpy(t->seq, pos + GRNA_L + 3, pos + GRNA_L + 8), 6);

            char *down = reverseComplement(slice_str_cpy(t->seq, pos - 6, pos - 1), 6);

            Gloc *g = Gloc_new_args(ngg, strcpyn(t->contig, strlen(t->contig)), up, down, pos + 3, pos + GRNA_L + 3, '-');

            Trie_set(t->trie, rcd, g, t->fp, write_callback);
            Gloc_del(g);
        }
        else
            free(rcd);
    }
}

void findNGG(struct ahocorasick *aho, Trie *trie, char *contig, char *target, FILE* fp, void(callback)(Trie *trie, struct aho_match_t *m))
{
    struct tmp *arg = (struct tmp *)malloc(sizeof(struct tmp));
    arg->seq = target;
    arg->trie = trie;
    arg->contig = contig;
    arg->len = strlen(target);
    arg->fp = fp;
    aho_register_match_callback(aho, callback, (void *)arg);
    aho_findtext(aho, target, strlen(target));
}


static int write_callback(FILE* fp, Gloc* gloc, int uid)
{
    char *insert = (char*) malloc(60);
    memset(insert,0,60);
    int write_bytes = sprintf(insert, "%s\t%d\t%d\t%s\t%d\t%c\t%s\t%s\n", gloc->contig, gloc->start, gloc->end, gloc->pam, uid, gloc->strand, gloc->up, gloc->dw);
    if (fwrite(insert ,1, strlen(insert), fp) != write_bytes){
        fprintf(stderr, "Write Failure");
    }

    free(insert);
}
