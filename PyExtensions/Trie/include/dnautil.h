/* Some stuff that you'll likely need in any program that works with
 * char.  Includes stuff for amino acids as well. 
 *
 * Assumes that char is stored as a character.
 * The char it generates will include the bases 
 * as lowercase tcag.  It will generally accept
 * uppercase as well, and also 'n' or 'N' or '-'
 * for unknown bases. 
 *
 * Amino acids are stored as single character upper case. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial.
 * 
 * ---------------------------------------------------
 * The file has been minified and modified by Ziwei Xue, 
 * introducing struct Gloc to store genomic location data.
 * Last Modified: 04MAR2021
 * Contact: xueziweisz@gmail.com
 */
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "../ahocorasick/ahocorasick.h"
#include "trie.h"

#define GRNA_L 20
typedef bool boolean;
typedef unsigned long long bits64;
typedef unsigned bits32;
typedef unsigned short bits16;
typedef unsigned char UBYTE;

#ifndef charUTIL_H
#define charUTIL_H
#endif

/* Generic string compare function */
int strngcmp(const char *__s1, const char *__s2, size_t __n);

/* make length of len */
void upper_str(char *__s, size_t __size);

/* slice of a string */
void slice_str(const char *__s, char *buf, size_t start, size_t end);

/* slice of a string and copy it */
char *slice_str_cpy(const char *__s, size_t start, size_t end);

/* concat two string and copy them */
char* cat_str(char *__s1, char* __s2);

/* copy an ngg sequence */
char *nggcpy(char *__s);

/* copy n elements of a sequence */
char *strcpyn(char* __s, size_t __sz);

void callback_match_pos(void *arg, struct aho_match_t *m);

void findNGG(struct ahocorasick *aho, Trie *trie, char *contig, char *target, FILE* fp, void(callback)(Trie *trie, struct aho_match_t *m));

void dnaUtilOpen(void); /* Good idea to call this before using any arrays here */

/* Numerical values for bases. */
#define T_BASE_VAL 0
#define U_BASE_VAL 0
#define C_BASE_VAL 1
#define A_BASE_VAL 2
#define G_BASE_VAL 3
#define N_BASE_VAL 4 /* Used in 1/2 byte representation. */

/* Initializing a new Gloc struct without any field */
Gloc *Gloc_new(void);

/* Free a Gloc struct */
void Gloc_del(Gloc *gloc);

/* Initializing a new with fields */
Gloc *Gloc_new_args(char* pam_, char *contig_, char *up_, char *dw_, size_t st, size_t en, char sd);

/* Setting the Gloc struct with fields*/

void Gloc_set(Gloc *gloc, char* pam_, char *contig_, char *up_, char *dw_, size_t st, size_t en, char sd);

/**
 * Map the start pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_start(Gloc *gloc, void *(*func)(void *));

/**
 * Map the end pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_end(Gloc *gloc, void *(*func)(void *));

/**
 * Map the strand pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_strand(Gloc *gloc, void *(*func)(void *));

/**
 * Map the pam char pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_pam(Gloc *gloc, void *(*func)(void *));

/**
 * Map the contig char pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_contig(Gloc *gloc, void *(*func)(void *));


/**
 * Map the up char pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_up(Gloc *gloc, void *(*func)(void *));

/**
 * Map the dw char pointer and make it points to a new 
 * address. A function pointer is passed as arguments
 * for conveniency in mapping to Python object.
 * */
void Gloc_map_dw(Gloc *gloc, void *(*func)(void *));


/* Get the start of the Gloc struct. */
size_t Gloc_get_start(Gloc *gloc);

/* Get the end position of the Gloc struct. */
size_t Gloc_get_end(Gloc *gloc);

/* Get the pam sequence of the Gloc struct. */
char* Gloc_get_pam(Gloc *gloc);

/* Get the contig of the Gloc struct. */
char *Gloc_get_contig(Gloc *gloc);

/* Get the upstream sequence of the Gloc struct. */
char *Gloc_get_up(Gloc *gloc);

/* Get the downstream of the Gloc struct. */
char *Gloc_get_dw(Gloc *gloc);

/* Get the strand of the Gloc struct. */
char Gloc_get_strand(Gloc *gloc);

/* Get the next element of the Gloc struct. */
Gloc* Gloc_get_next(Gloc *gloc);

/* Set the next element of the Gloc struct. */
void Gloc_set_next(Gloc *gloc, Gloc *value);

/* Initialize the nucleotides compare characters table. */
void initNtCompTable(void);

/* Initialize the mixed case characters table. */
void initNtMixedCaseChars(void);

/* Reverse the order of the bytes. */
void reverseBytes(char *bytes, long length);

/* Complement char (not reverse). */
char *complement(char *dna, long length);

/* Reverse complement char. */
char *reverseComplement(char *dna, long length);

/** checks whether the gRNA sequences is a 
 * valid sequence. If characters rather than
 * A,T,C,G appear in the sequence, return false
**/ 
boolean isValidDna(char *poly, int size);

/** checks whether the DNA sequence is a valid 
 * DNA sequence. If characters rather than IAPUC Ambiguty codes
 * appears in the sequence, return false
**/
boolean isValidAmbiguityCodes(char *poly, int size);

/** pre-filtering the grna being inserted to the trie tree.
**/
boolean isValidGRNA(char* seq, int thres);

/** write callback function while building the trie.
**/ 
static int write_callback(FILE* fp, Gloc* gloc, int uid);