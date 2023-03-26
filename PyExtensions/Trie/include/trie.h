#include <stdio.h>
/* The trie.h header file is from BioPython source code.
 * 
 * ---------------------------------------------------
 * The file has been minified and modified by Ziwei Xue, 
 * making it available to store GRNA sequence and introducing
 * new field to store genomic location.
 * Last Modified: 05MAR2021
 * Contact: xueziweisz@gmail.com
 */

/* forward declaration */
typedef struct Gloc Gloc;
typedef struct Trie Trie;			  /* forward declaration */
typedef struct Transition Transition; /* forward declaration */
/* Trie_new
 * --------
 * Create a new trie.  Return a Trie structure, which is an abstract
 * data structure.  The client should not have to know about the
 * details of this structure.  When finished, each Trie should be
 * freed with Trie_del.
 */
Trie *Trie_new(void);

/* Trie_del
 * --------
 * Free a Trie data structure.
 */
void Trie_del(Trie *trie);

/* Trie_set
 * --------
 * Set a string in the Trie to some value.  Returns a 0 if the
 * function succeeded.
 */
int Trie_set(Trie* trie, const char *key, const void *value, FILE* fp, int (callback) (FILE* fp, Gloc* gloc, int uid));


int Trie_set_pyval(Trie *trie, const char *key, const void *value);

/* Trie_get
 * --------
 * Lookup whether a key exists in the Trie.  Returns the value that
 * was previous set in the Trie, or NULL if it doesn't exist.
 */
void *Trie_get(const Trie *trie, const char *key);

/* Trie_get_pyval
 * --------
 * Lookup whether a key exists in the Trie.  Returns the value that
 * was previous set in the Trie, or NULL if it doesn't exist.
 */
void *Trie_get_pyval(const Trie *trie, const char *key);

/* Trie_get_data
 * --------
 * Lookup whether a key exists in the Trie.  Returns the value that
 * was previous set in the Trie, or NULL if it doesn't exist.
 */
void *Trie_get_data(const Trie *trie);

/* Trie_get_approximate
 * --------------------
 * Lookup whether a key exists in the Trie, allowing for mismatches to
 * the dictionary (Levenshtein distance from key to strings in dictionary).
 * Passes back values using a callback function.
 */
void Trie_get_approximate(const Trie *trie, const char *key, const int k,
						  void (*callback)(const char *key,
										   const void *value,
										   const int mismatches,
										   void *data),
						  void *data);

/* add code to Biopython original code
   in order to implement new function get_approximate_hamming;
   copy all declarations corresponding to get_approximate */

/* Trie_get_approximate_hamming
 * ----------------------------
 * Lookup whether a key exists in the Trie, allowing for mismatches to
 * the dictionary (Hamming distance from key to strings in dictionary).
 * Passes back values using a callback function.
 */
void Trie_get_approximate_hamming(const Trie *trie, const char *key, const int k,
								  void (*callback)(const char *key,
												   const void *value,
												   const int mismatches,
												   void *data),
								  void *data);

/* end of added code for get_approximate_hamming */

/* Trie_len
 * --------
 * Return the number of strings in the trie.
 */
int Trie_len(const Trie *trie);

/* Trie_has_key
 * ------------
 * Return whether a key exists in the trie.
 */
int Trie_has_key(const Trie *trie, const char *key);

/* Trie_has_prefix
 * ---------------
 * Return whether a string is a prefix of a key in the trie.
 */
int Trie_has_prefix(const Trie *trie, const char *prefix);

/* Trie_with_prefix
 * ----------------
 * Iterate over all the keys in the trie that start with a prefix.
 */
void Trie_with_prefix(const Trie *trie, const char *prefix,
					  void (*callback)(const char *key,
									   const void *value,
									   void *data),
					  void *data);

/* Trie_iterate
 * ------------
 * Iterate through everything stored in the trie.  callback is a
 * function that gets called for each thing in the trie.  It is called
 * in arbitrary order.  data is a pointer to some arbitrary data and
 * gets passed unchanged to the callback.
 */
void Trie_iterate(const Trie *trie,
				  void (*callback)(const char *key,
								   const void *value,
								   void *data),
				  void *data);

/* Trie_serialize
 * --------------
 * Serialize a tree into a stream of bytes.  This function takes a
 * callback 'write' that should take a pointer to data and the length
 * of the data in bytes.  This will be called repeatedly until the
 * whole Trie is serialized.  When it is done, this function will call
 * 'write' with a length of 0.  Since the values are handled by the
 * client, this function also takes a callback function 'write_value'
 * so that the client can serialize their own values.
 *
 * This function is platform-dependent, so byte streams created on one
 * machine may not necessarily port to another.
 */
int Trie_serialize(const Trie *trie,
				   int (*write)(const void *towrite, const int length,
								void *data),
				   int (*write_value)(const void *value, void *data),
				   void *data);

/* Trie_deserialize
 * ----------------
 * Deserialize a tree that was previously serialized with
 * Trie_serialize.  This function takes a callback 'read' that should
 * read 'length' bytes and save it to 'wasread'.  'read_value' should
 * read a value and return a pointer to it.  'data' is a pointer that
 * will be passed unchanged to 'read' and 'read_value'.
 */
Trie *Trie_deserialize(int (*read)(void *wasread, const int length, void *data),
					   void *(*read_value)(void *data),
					   void *data);

/** Trie_traverse_trie
 * ----------------
 * traverse a trie and save the pointer to 
 * the tries to an array.
**/
Trie **Trie_traverse_trie(const Trie *trie);

/** Trie_traverse_transition
 * ----------------
 * traverse a trie and save the keys to 
 * an array.
 */
char **Trie_traverse_transition(const Trie *trie);

/** Trie_map_trie_data
 * ----------------
 * traverse a trie and map the trie values by 
 * a function pointer.
 */
int Trie_map_trie_data(Trie *trie, int (*func)(void *));

/** Trie_map_trie_trie
 * ----------------
 * traverse a trie and map the trie by 
 * a function pointer.
 */
int Trie_map_trie(Trie *trie, int (*func)(void *));

/** Trie_merge
 * ----------------
 * Merge two tries and returns a new Trie.
 */
Trie *Trie_merge(const Trie *trie_l, const Trie *trie_r);

/** Trie_size
 * ----------------
 * returns the total size (length of
 * the key)
 */
int Trie_size(const Trie *trie);

/** returns the exact trie given the prefix
 */ 
Trie *Trie_get_trie(const Trie *trie, const char *key);

/** implementation of the lazy iterator of the
 * trie data structure. Required for making 
 * the python trie itself to be a generator object
**/ 
void *Trie_lazy_iter(Trie *trie);