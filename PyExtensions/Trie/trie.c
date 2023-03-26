#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include "include/trie.h"
#include "include/triemodule.h"
#include "include/dnautil.h"

/** duplicate is originally in this file.
 */
static char *duplicate(const char *__s)
{
    /* Don't use strdup, as it's not ANSI C. */
    char *t = malloc(strlen(__s) + 1);
    memset(t, 0, strlen(__s)+1);
    if (!t)
    // error in allocating memory
        return NULL;
    // better not to use strcpy
    // strcpy(t, s);
    memcpy(t, __s, strlen(__s));
    return t;
}

/* Transition holds information about the transitions leading from
 * one Trie to another.  The trie structure here is different from
 * typical ones, because the transitions between nodes can contain
 * strings of arbitrary length, not just single characters.  Suffix is
 * the string that is matched from one node to the next.
 */
struct Transition
{
    char *suffix;
    Trie *next;
    Trie *prev;
};

/* Trie is a recursive data structure.  A Trie contains zero or more
 * Transitions that lead to more Tries.  The transitions are stored
 * in alphabetical order of the suffix member of the data structure.
 * Trie also contains a pointer called value where the user can store
 * arbitrary data.  If value is NULL, then no data is stored here.
 */
struct Trie
{
    Transition *transitions;
    Transition *parent_transition;
    unsigned char num_transitions;
    size_t size;
    void *value; /* never set this field in python! */
    int setted;
    void *py_value;
    unsigned long long uid;
    Trie *cur_trie; // for iter
    char* cur_suffix; // for iter
};

#define MAX_KEY_LENGTH (200)
static char KEY[MAX_KEY_LENGTH];

Trie *Trie_new(void)
{
    Trie *trie;

    if (!(trie = malloc(sizeof(struct Trie))))
        return NULL;
    trie->transitions = NULL;
    trie->parent_transition = NULL;
    trie->num_transitions = 0;
    trie->size = 0;
    trie->setted = 0;
    trie->value = NULL;
    trie->py_value = NULL;
    trie->cur_trie = NULL;
    return trie;
}


void *Trie_get_data(const Trie* trie)
{
    return trie->value;
}

Trie *Trie_get_parent_transition(const Trie *trie)
{
    return trie->parent_transition;
}

Trie *Trie_get_parent_trie(const Trie *trie)
{
    return trie->parent_transition->prev;
}

Trie *Trie_get_trie(const Trie *trie, const char *key)
{
    int first, last, mid;

    if (!key[0])
    {
        return trie;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        char *suffix;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        /* If suffix is a substring of key, then get the value from
	   the next trie.
	*/
        c = strncmp(key, suffix, strlen(suffix));
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
            return Trie_get_trie(transition->next, key + strlen(suffix));
    }
    return NULL;
}


void *Trie_get_pyval(const Trie* trie, const char* key)
{
    int first, last, mid;

    if (!key[0])
    {
        return trie->py_value;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        char *suffix;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        /* If suffix is a substring of key, then get the value from
	   the next trie.
	*/
        c = strncmp(key, suffix, strlen(suffix));
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
            return Trie_get_pyval(transition->next, key + strlen(suffix));
    }
    return NULL;
}

unsigned long long Trie_get_uid(const Trie* trie, const char* key)
{
    int first, last, mid;

    if (!key[0])
    {
        return trie->uid;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        char *suffix;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        /* If suffix is a substring of key, then get the value from
	   the next trie.
	*/
        c = strncmp(key, suffix, strlen(suffix));
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
            return Trie_get_uid(transition->next, key + strlen(suffix));
    }
    return NULL;
}

int _trie_set_helper(Trie *head, Trie *trie, const char *key, const void *value, FILE* fp, int (callback) (FILE* fp, Gloc* gloc, int uid))
{
    int i;
    Transition *transition = NULL;
    const char *suffix = NULL;
    int retval = 0;
    int first, last, mid;
    /**
     * if the key is a empty array
     */
    if (!key[0])
    {
        if (trie->setted && value) callback(fp, value, trie->uid);
        else {
            if (value) {
                trie->value = PyLong_FromLong( trie->uid );
                callback(fp, value, trie->uid);
                trie->setted = 1;
            }
        }
        // trie->parent_transition = parent_transition;
        return 0;
    }

    /* Insert the key in alphabetical order.  Do a binary search to
       find the proper place. */
    first = 0;
    last = trie->num_transitions - 1;
    i = -1;
    while (first <= last)
    {
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        if (key[0] < suffix[0])
            last = mid - 1;
        else if (key[0] > suffix[0])
            first = mid + 1;
        else
        {
            i = mid;
            break;
        }
    }
    /* If no place was found for it, then the indexes will be in the
       order last,first.  Place it at index first. */
    if (i == -1)
        i = first;
    /* If nothing matches, then insert a new trie here. */
    if ((i >= trie->num_transitions) || (key[0] != suffix[0]))
    {
        char *new_suffix = NULL;
        Trie *newtrie = NULL;
        Transition *new_transitions = NULL;
        /* Create some variables for the new transition.  I'm going to
	       allocate these firsst so that if I can detect memory errors
	       before I mess up the data structure of the transitions.
	    */
        if (!(new_suffix = duplicate(key)))
            goto insert_memerror;
        if (!(newtrie = Trie_new()))
            goto insert_memerror;

        /* Create some space for the next transition.  Allocate some
	       memory and shift the old transitions over to make room for
	      this one.
	    */
        if (!(new_transitions = malloc(sizeof(Transition) *
                                       (trie->num_transitions + 1))))
            goto insert_memerror;
        memcpy(new_transitions, trie->transitions, sizeof(Transition) * i);
        memcpy(&new_transitions[i + 1], &trie->transitions[i], sizeof(Transition) * (trie->num_transitions - i));
        free(trie->transitions);
        trie->transitions = new_transitions;
        for (int j =0; j < trie->num_transitions;++j)
        {
            if (i != j) {
                // trie->transitions[j].next->parent_transition = &trie->transitions[j];
                // trie->transitions[j].prev = trie;
            }
        }
        
        new_transitions = NULL;
        trie->num_transitions += 1;
        /* Initialize the new transition. */
        transition = &trie->transitions[i];
        // transition->prev = trie;
        transition->suffix = new_suffix;
        transition->next = newtrie;
        // newtrie->parent_transition = transition;
        if (transition->next->setted && value) callback(fp, value, transition->next->uid);
        else {
            if (value) {
                transition->next->setted = 1;
                transition->next->uid = head->size + 1;
                transition->next->value = PyLong_FromLong( transition->next->uid );
                callback(fp, value, transition->next->uid);
            }
        }

        if (0)
        {
        insert_memerror:
            if (new_transitions)
                free(new_transitions);
            if (newtrie)
                free(newtrie);
            if (new_suffix)
                free(new_suffix);
            return 1;
        }
        head->size += 1;
    }
    /* There are three cases where the key and suffix share some
       letters.
       1.  suffix is proper substring of key.
       2.  key is proper substring of suffix.
       3.  neither is proper substring of other.

       For cases 2 and 3, I need to first split up the transition
       based on the number of characters shared.  Then, I can insert
       the rest of the key into the next trie.
    */
    else
    {
        /* Count the number of characters shared between key
	   and suffix. */
        int chars_shared = 0;
        while (key[chars_shared] && key[chars_shared] == suffix[chars_shared])
            chars_shared++;
        /* Case 2 or 3, split this sucker! */
        if (chars_shared < strlen(suffix))
        {
            Trie *newtrie = NULL;
            char *new_suffix1 = NULL, *new_suffix2 = NULL;
            if (!(new_suffix1 = malloc(chars_shared + 1)))
                goto split_memerror;
            strncpy(new_suffix1, key, chars_shared);
            new_suffix1[chars_shared] = 0;
            if (!(new_suffix2 = duplicate(suffix + chars_shared)))
                goto split_memerror;
            if (!(newtrie = Trie_new()))
                goto split_memerror;
            if (!(newtrie->transitions = malloc(sizeof(Transition))))
                goto split_memerror;
            newtrie->num_transitions = 1;
            newtrie->transitions[0].next = transition->next;
            // newtrie->transitions[0].next->parent_transition = &newtrie->transitions[0];
            newtrie->transitions[0].suffix = new_suffix2;
            // newtrie->transitions[0].prev = newtrie;
            
            free(transition->suffix);
            transition->suffix = new_suffix1;
            transition->next = newtrie;
            newtrie->parent_transition = transition;
            if (0)
            {
            split_memerror:
                if (newtrie && newtrie->transitions)
                    free(newtrie->transitions);
                if (newtrie)
                    free(newtrie);
                if (new_suffix2)
                    free(new_suffix2);
                if (new_suffix1)
                    free(new_suffix1);
                return 1;
            }
        }
        retval = _trie_set_helper(head, transition->next, key + chars_shared, value, fp, callback);
    }

    return retval;
}

int Trie_set(Trie* trie, const char *key, const void *value, FILE* fp, int (callback) (FILE* fp, Gloc* gloc, int uid))
{
    return _trie_set_helper(trie, trie, key, value, fp, callback);
}

int _trie_set_pyval_helper(Trie *trie, const char *key, const void *value)
{
    Trie* t;
    if (!(t = Trie_get_trie(trie, key))) return 0;
    t->py_value = value;
    return 1;
}

int Trie_set_pyval(Trie* trie, const char *key, const void *value)
{
    return _trie_set_pyval_helper(trie, key, value);
}

void Trie_del(Trie *trie)
{
    int i;
    if (!trie)
        return;
    for (i = 0; i < trie->num_transitions; i++)
    {
        Transition *transition = &trie->transitions[i];
        if (transition->suffix)
            free(transition->suffix);
        Trie_del(transition->next);
    }
    free(trie);
}

void *Trie_get(const Trie *trie, const char *key)
{
    int first, last, mid;

    if (!key[0])
    {
        return trie->value;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        char *suffix;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        /* If suffix is a substring of key, then get the value from
	   the next trie.
	*/
        c = strncmp(key, suffix, strlen(suffix));
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
            return Trie_get(transition->next, key + strlen(suffix));
    }
    return NULL;
}

/* DEBUG
static void
_print_trie_helper(const Trie* trie, int num_indent) {
    int i, j;
    char *message = "no value";
    Transition *t;

    if(trie->value != NULL)
        message = "has value";
    for(j=0; j<num_indent; j++)
        printf(" ");
    printf("%d transitions, %s.\n", trie->num_transitions, message);
    for(i=0; i<trie->num_transitions; i++) {
        t = &trie->transitions[i];
        for(j=0; j<num_indent; j++)
            printf(" ");
        printf("suffix %s\n", t->suffix);
        _print_trie_helper(t->next, num_indent+2);
    }
}

static void
_print_trie(const Trie* trie) {
    _print_trie_helper(trie, 0);
}
*/

/* Mutually recursive, so need to make a forward declaration. */
static void
_get_approximate_trie(const Trie *trie, const char *key, const int k,
                      void (*callback)(const char *key,
                                       const void *value,
                                       const int mismatches,
                                       void *data),
                      void *data,
                      const int mismatches,
                      char *current_key, const int max_key);

static void
_get_approximate_transition(const char *key,
                            const int k,
                            const Transition *transition,
                            const char *suffix,
                            void (*callback)(const char *key,
                                             const void *value,
                                             const int mismatches,
                                             void *data),
                            void *data,
                            const int mismatches,
                            char *current_key, const int max_key)
{
    int i;
    int prev_keylen = strlen(current_key);

    /* Short circuit optimization.  If there's too many characters to
       possibly be a match, then don't even try to match things. */
    if ((int)(strlen(suffix) - strlen(key)) > k)
        return;

    /* Match as many characters as possible. */
    i = 0;
    while (suffix[i] && (key[i] == suffix[i]))
    {
        i++;
    }
    /* Check to make sure the key is not too long.  BUG: If it is,
       fails silently. */
    if ((prev_keylen + i) >= max_key)
        return;
    strncat(current_key, suffix, i);

    /* If all the letters in the suffix matched, then move to the
       next trie. */
    if (!suffix[i])
    {
        _get_approximate_trie(transition->next, &key[i], k, callback, data,
                              mismatches, current_key, max_key);
    }
    /* Otherwise, try out different kinds of mismatches. */
    else if (k)
    {
        int new_keylen = prev_keylen + i;

        /* Letter replacement, skip the next letter in both the key and
	   suffix. */
        if ((new_keylen + 1 < max_key) && key[i] && suffix[i])
        {
            current_key[new_keylen] = suffix[i];
            current_key[new_keylen + 1] = 0;
            _get_approximate_transition(&key[i + 1], k - 1,
                                        transition, &suffix[i + 1],
                                        callback, data,
                                        mismatches + 1, current_key, max_key);
            current_key[new_keylen] = 0;
        }

        /* Insertion in key, skip the next letter in the key. */
        if (key[i])
        {
            _get_approximate_transition(&key[i + 1], k - 1,
                                        transition, &suffix[i],
                                        callback, data,
                                        mismatches + 1, current_key, max_key);
        }

        /* Deletion from key, skip the next letter in the suffix. */
        if ((new_keylen + 1 < max_key) && suffix[i])
        {
            current_key[new_keylen] = suffix[i];
            current_key[new_keylen + 1] = 0;
            _get_approximate_transition(&key[i], k - 1,
                                        transition, &suffix[i + 1],
                                        callback, data,
                                        mismatches + 1, current_key, max_key);
            current_key[new_keylen] = 0;
        }
    }
    current_key[prev_keylen] = 0;
}

static void
_get_approximate_trie(const Trie *trie, const char *key, const int k,
                      void (*callback)(const char *key,
                                       const void *value,
                                       const int mismatches,
                                       void *data),
                      void *data,
                      const int mismatches,
                      char *current_key, const int max_key)
{
    int i;

    /* If there's no more key to match, then I'm done. */
    if (!key[0])
    {
        if (trie->value)
            (*callback)(current_key, trie->value, mismatches, data);
    }
    /* If there are no more mismatches allowed, then fall back to the
       faster Trie_get. */
    else if (!k)
    {
        void *value = Trie_get(trie, key);
        if (value)
        {
            int l = strlen(current_key);
            /* Make sure I have enough space for the full key. */
            if (l + strlen(key) < max_key)
            {
                strcat(current_key, key);
                (*callback)(current_key, value, mismatches, data);
                current_key[l] = 0;
            }
            /* BUG: Ran out of space for the key.  This fails
	       silently, but should signal an error. */
        }
    }
    /* If there are no more transitions, then all the characters left
       in the key are mismatches. */
    else if (!trie->num_transitions)
    {
        if (trie->value && (strlen(key) <= k))
        {
            (*callback)(current_key, trie->value,
                        mismatches + strlen(key), data);
        }
    }
    /* Otherwise, try to match each of the transitions. */
    else
    {
        for (i = 0; i < trie->num_transitions; i++)
        {
            Transition *transition = &trie->transitions[i];
            const char *suffix = transition->suffix;
            _get_approximate_transition(key, k, transition, suffix,
                                        callback, data,
                                        mismatches, current_key, max_key);
        }
    }
}

void Trie_get_approximate(const Trie *trie, const char *key, const int k,
                          void (*callback)(const char *key,
                                           const void *value,
                                           const int mismatches,
                                           void *data),
                          void *data)
{
    KEY[0] = 0;
    _get_approximate_trie(trie, key, k, callback, data, 0, KEY, MAX_KEY_LENGTH);
}

/* add code to Biopython original code
   in order to implement new function get_approximate_hamming;
   copy all code corresponding to get_approximate,
   and modify only _get_approximate_transition() where
   two options for traversing the trie are removed,
   making it to correspond to Hamming distance rather than
   edit (Levenshtein) distance; also this version becomes much faster */

/* Mutually recursive, so need to make a forward declaration. */
static void
_get_approximate_hamming_trie(const Trie *trie, const char *key, const int k,
                              void (*callback)(const char *key,
                                               const void *value,
                                               const int mismatches,
                                               void *data),
                              void *data,
                              const int mismatches,
                              char *current_key, const int max_key);

static void
_get_approximate_hamming_transition(const char *key,
                                    const int k,
                                    const Transition *transition,
                                    const char *suffix,
                                    void (*callback)(const char *key,
                                                     const void *value,
                                                     const int mismatches,
                                                     void *data),
                                    void *data,
                                    const int mismatches,
                                    char *current_key, const int max_key)
{
    int i;
    int prev_keylen = strlen(current_key);

    /* Short circuit optimization.  If there's too many characters to
       possibly be a match, then don't even try to match things. */
    if ((int)(strlen(suffix) - strlen(key)) > k)
        return;

    /* Match as many characters as possible. */
    i = 0;
    while (suffix[i] && (key[i] == suffix[i]))
    {
        i++;
    }
    /* Check to make sure the key is not too long.  BUG: If it is,
       fails silently. */
    if ((prev_keylen + i) >= max_key)
        return;
    strncat(current_key, suffix, i);

    /* If all the letters in the suffix matched, then move to the
       next trie. */
    if (!suffix[i])
    {
        _get_approximate_hamming_trie(transition->next, &key[i], k, callback, data,
                                      mismatches, current_key, max_key);
    }
    /* Otherwise, try out different kinds of mismatches. */
    else if (k)
    {
        int new_keylen = prev_keylen + i;

        /* Letter replacement, skip the next letter in both the key and
       suffix. */
        if ((new_keylen + 1 < max_key) && key[i] && suffix[i])
        {
            current_key[new_keylen] = suffix[i];
            current_key[new_keylen + 1] = 0;
            _get_approximate_hamming_transition(&key[i + 1], k - 1,
                                                transition, &suffix[i + 1],
                                                callback, data,
                                                mismatches + 1, current_key, max_key);
            current_key[new_keylen] = 0;
        }

        /* Insertion in key, skip the next letter in the key. * /
    if(key[i]) {
        _get_approximate_hamming_transition(&key[i+1], k-1,
                    transition, &suffix[i],
                    callback, data,
                    mismatches+1, current_key, max_key);
    }
    */
        /* Deletion from key, skip the next letter in the suffix. * /
    if((new_keylen+1 < max_key) && suffix[i]) {
        current_key[new_keylen] = suffix[i];
        current_key[new_keylen+1] = 0;
        _get_approximate_hamming_transition(&key[i], k-1,
                    transition, &suffix[i+1],
                    callback, data,
                    mismatches+1, current_key, max_key);
        current_key[new_keylen] = 0;
    }
*/
    }
    current_key[prev_keylen] = 0;
}

static void
_get_approximate_hamming_trie(const Trie *trie, const char *key, const int k,
                              void (*callback)(const char *key,
                                               const void *value,
                                               const int mismatches,
                                               void *data),
                              void *data,
                              const int mismatches,
                              char *current_key, const int max_key)
{
    int i;

    /* If there's no more key to match, then I'm done. */
    if (!key[0])
    {
        if (trie->value)
            (*callback)(current_key, trie->value, mismatches, data);
    }
    /* If there are no more mismatches allowed, then fall back to the
       faster Trie_get. */
    else if (!k)
    {
        void *value = Trie_get(trie, key);
        if (value)
        {
            int l = strlen(current_key);
            /* Make sure I have enough space for the full key. */
            if (l + strlen(key) < max_key)
            {
                strcat(current_key, key);
                (*callback)(current_key, value, mismatches, data);
                current_key[l] = 0;
            }
            /* BUG: Ran out of space for the key.  This fails
           silently, but should signal an error. */
        }
    }
    /* If there are no more transitions, then all the characters left
       in the key are mismatches. */
    else if (!trie->num_transitions)
    {
        if (trie->value && (strlen(key) <= k))
        {
            (*callback)(current_key, trie->value,
                        mismatches + strlen(key), data);
        }
    }
    /* Otherwise, try to match each of the transitions. */
    else
    {
        for (i = 0; i < trie->num_transitions; i++)
        {
            Transition *transition = &trie->transitions[i];
            const char *suffix = transition->suffix;
            _get_approximate_hamming_transition(key, k, transition, suffix,
                                                callback, data,
                                                mismatches, current_key, max_key);
        }
    }
}

void Trie_get_approximate_hamming(const Trie *trie, const char *key, const int k,
                                  void (*callback)(const char *key,
                                                   const void *value,
                                                   const int mismatches,
                                                   void *data),
                                  void *data)
{
    KEY[0] = 0;
    _get_approximate_hamming_trie(trie, key, k, callback, data, 0, KEY, MAX_KEY_LENGTH);
}

/* end of added code for get_approximate_hamming */

int Trie_len(const Trie *trie)
{
    int i;
    int length = 0;

    if (!trie)
        return 0;
    if (trie->value)
        length += 1;
    for (i = 0; i < trie->num_transitions; i++)
    {
        length += Trie_len(trie->transitions[i].next);
    }
    return length;
}

int Trie_has_key(const Trie *trie, const char *key)
{
    return Trie_get(trie, key) != NULL;
}

int Trie_has_prefix(const Trie *trie, const char *prefix)
{
    int first, last, mid;

    if (!prefix[0])
    {
        return 1;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        char *suffix;
        int suffixlen, prefixlen, minlen;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        suffixlen = strlen(suffix);
        prefixlen = strlen(prefix);
        minlen = (suffixlen < prefixlen) ? suffixlen : prefixlen;
        c = strncmp(prefix, suffix, minlen);
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
            return Trie_has_prefix(transition->next, prefix + minlen);
    }
    return 0;
}

static void
_iterate_helper(const Trie *trie,
                void (*callback)(const char *key,
                                 const void *value,
                                 void *data),
                void *data,
                char *current_key, const int max_key)
{
    int i;
    if (trie->value)
        (*callback)(current_key, trie->value, data);
    for (i = 0; i < trie->num_transitions; i++)
    {
        Transition *transition = &trie->transitions[i];
        const char *suffix = transition->suffix;
        int keylen = strlen(current_key);

        if (keylen + strlen(suffix) >= max_key)
        {
            /* BUG: This will fail silently.  It should raise some
	       sort of error. */
            continue;
        }
        strcat(current_key, suffix);
        _iterate_helper(transition->next, callback, data,
                        current_key, max_key);
        current_key[keylen] = 0;
    }
}

void Trie_iterate(const Trie *trie,
                  void (*callback)(const char *key,
                                   const void *value,
                                   void *data),
                  void *data)
{
    KEY[0] = 0;
    _iterate_helper(trie, callback, data, KEY, MAX_KEY_LENGTH);
}

static void
_with_prefix_helper(const Trie *trie, const char *prefix,
                    void (*callback)(const char *key,
                                     const void *value,
                                     void *data),
                    void *data,
                    char *current_key, const int max_key)
{
    int first, last, mid;

    if (!prefix[0])
    {
        _iterate_helper(trie, callback, data, current_key, max_key);
        return;
    }

    /* The transitions are stored in alphabetical order.  Do a binary
     * search to find the proper one.
     */
    first = 0;
    last = trie->num_transitions - 1;
    while (first <= last)
    {
        Transition *transition;
        const char *suffix;
        int suffixlen, prefixlen, minlen;
        int c;
        mid = (first + last) / 2;
        transition = &trie->transitions[mid];
        suffix = transition->suffix;
        suffixlen = strlen(suffix);
        prefixlen = strlen(prefix);
        minlen = (suffixlen < prefixlen) ? suffixlen : prefixlen;
        c = strncmp(prefix, suffix, minlen);
        if (c < 0)
            last = mid - 1;
        else if (c > 0)
            first = mid + 1;
        else
        {
            /* Three cases here.
             * 1.  Suffix and prefix are the same.
             *     Add suffix to current_key, advance prefix to the
             *     end.  Continue recursively.  Since there is no more
             *     prefix, every sub-trie will be considered as having
             *     this prefix.
             * 2.  Suffix shorter than prefix.
             *     suffix  A
             *     prefix  AN
             *     Add suffix (A) to current_key.  Advance prefix by
             *     1.  Continue recursively to match rest of prefix.
             * 3.  Suffix longer than prefix.
             *     suffix  AN
             *     prefix  A
             *     Add suffix (AN) to current_key.  Advance prefix to
             *     the end.  Continue recursively.  Since there is no
             *     more prefix, every sub-trie will be considered as
             *     having this prefix.
             */
            int keylen = strlen(current_key);
            if (keylen + suffixlen >= max_key)
            {
                /* BUG: This will fail silently.  It should raise some
		   sort of error. */
                break;
            }
            strncat(current_key, suffix, suffixlen);
            _with_prefix_helper(transition->next, prefix + minlen,
                                callback, data, current_key, max_key);
            current_key[keylen] = 0; /* reset current_key */
            break;
        }
    }
}

void Trie_with_prefix(const Trie *trie, const char *prefix,
                      void (*callback)(const char *key,
                                       const void *value,
                                       void *data),
                      void *data)
{
    /*_print_trie(trie);*/
    KEY[0] = 0;
    _with_prefix_helper(trie, prefix, callback, data, KEY, MAX_KEY_LENGTH);
}

/* Need to declare _serialize_transition here so it can be called from
   _serialize_trie. */
static int _serialize_transition(const Transition *transition,
                                 int (*write)(const void *towrite, const int length,
                                              void *data),
                                 int (*write_value)(const void *value, void *data),
                                 void *data);

/* This library also provides code for flattening tries so that they
 * can be saved and read back in later.  The format of a serialized
 * trie is:
 * TYPE        NBYTES    DESCRIPTION
 * byte        1         Whether or not there is a value
 * variable    variable  If there is a value, let the client store it.
 * byte        1         Whether or not there is a py value
 * variable    variable  If there is a pyvalue, let the client store it.
 * byte        1         Whether or not there is a uid
 * variable    variable  If there is a uid, let the client store it.
 * byte        1         Number of transitions for this Trie.
 * transition  variable
 *   int       4         Number of characters in the suffix.
 *   suffix    variable  the suffix for this transition
 *   byte      1         Whether or not there is a trie
 *   trie      variable  Recursively points to another trie.
 *
 * The number of bytes and the endian may vary from platform to
 * platform.
 */

static int _serialize_trie(const Trie *trie,
                           int (*write)(const void *towrite, const int length,
                                        void *data),
                           int (*write_value)(const void *value, void *data),
                           void *data)
{
    int i;
    unsigned char has_value;
    has_value = (trie->value != NULL);
    if (!(*write)(&has_value, sizeof(has_value), data))
        return 0;
    
    if (has_value)
    {
        if (!(*write_value)(trie->value, data))
            return 0;
    }
    has_value = (trie->py_value != NULL);
    if (!(*write)(&has_value, sizeof(has_value), data))
        return 0;

    if (has_value)  
    {
        if (!(*write_value)(trie->py_value, data))
            return 0;
    } 

    has_value = (trie->uid != NULL);
    if (!(*write)(&has_value, sizeof(has_value), data))
        return 0;

    if (has_value)  
    {
        if (!(*write_value)(PyLong_FromLong(trie->uid), data))
            return 0;
    } 


    if (!(*write)(&trie->num_transitions, sizeof(trie->num_transitions), data))
        return 0;
    for (i = 0; i < trie->num_transitions; i++)
    {
        if (!_serialize_transition(&trie->transitions[i],
                                   write, write_value, data))
            return 0;
    }

    return 1;
}

static int _serialize_transition(const Transition *transition,
                                 int (*write)(const void *towrite, const int length,
                                              void *data),
                                 int (*write_value)(const void *value, void *data),
                                 void *data)
{
    int suffixlen;
    unsigned char has_trie;

    suffixlen = strlen(transition->suffix);
    if (!(*write)(&suffixlen, sizeof(suffixlen), data))
        return 0;
    if (!(*write)(transition->suffix, suffixlen, data))
        return 0;

    has_trie = (transition->next != NULL);
    if (!(*write)(&has_trie, sizeof(has_trie), data))
        return 0;
    if (has_trie)
    {
        if (!_serialize_trie(transition->next, write, write_value, data))
            return 0;
    }
    return 1;
}

int Trie_serialize(const Trie *trie,
                   int (*write)(const void *towrite, const int length,
                                void *data),
                   int (*write_value)(const void *value, void *data),
                   void *data)
{
    int success = _serialize_trie(trie, write, write_value, data);
    (*write)(NULL, 0, data);
    return success;
}

static int _deserialize_transition(Transition *transition,
                                   int (*read)(void *wasread, const int length,
                                               void *data),
                                   void *(*read_value)(void *data),
                                   void *data);

static int _deserialize_trie(Trie *trie,
                             int (*read)(void *wasread, const int length, void *data),
                             void *(*read_value)(void *data),
                             void *data)
{
    int i;
    unsigned char has_value;
    if (!(*read)(&has_value, sizeof(has_value), data))
        goto _deserialize_trie_error;
    if (has_value != 0 && has_value != 1)
        goto _deserialize_trie_error;
    if (has_value)
    {
        if (!(trie->value = (*read_value)(data)))
            goto _deserialize_trie_error;
    }

    if (!(*read)(&has_value, sizeof(has_value), data))
        goto _deserialize_trie_error;
    if (has_value != 0 && has_value != 1)
        goto _deserialize_trie_error;
    if (has_value)
    {
        if (!(trie->py_value = (*read_value)(data)))
            goto _deserialize_trie_error;
    }

    if (!(*read)(&has_value, sizeof(has_value), data))
        goto _deserialize_trie_error;
    if (has_value != 0 && has_value != 1)
        goto _deserialize_trie_error;
    if (has_value)
    {
        if (!(trie->uid = PyLong_AsLong ((*read_value)(data))))
            goto _deserialize_trie_error;
    }

    if (!(*read)(&trie->num_transitions, sizeof(trie->num_transitions), data))
        goto _deserialize_trie_error;
    if (!(trie->transitions =
              malloc(trie->num_transitions * sizeof(Transition))))
        goto _deserialize_trie_error;
    for (i = 0; i < trie->num_transitions; i++)
    {
        trie->transitions[i].suffix = NULL;
        trie->transitions[i].next = NULL;
    }
    for (i = 0; i < trie->num_transitions; i++)
    {
        if (!_deserialize_transition(&trie->transitions[i],
                                     read, read_value, data))
            goto _deserialize_trie_error;
    }
    return 1;

_deserialize_trie_error:
    trie->num_transitions = 0;
    if (trie->transitions)
    {
        free(trie->transitions);
        trie->transitions = NULL;
    }
    trie->value = NULL;
    return 0;
}

static int _deserialize_transition(Transition *transition,
                                   int (*read)(void *wasread, const int length,
                                               void *data),
                                   void *(*read_value)(void *data),
                                   void *data)
{
    int suffixlen;
    unsigned char has_trie;

    if (!(*read)(&suffixlen, sizeof(suffixlen), data))
        goto _deserialize_transition_error;
    if (suffixlen < 0 || suffixlen >= MAX_KEY_LENGTH)
    {
        printf("MAX_KEY_LENGTH too short [%d:%d]\n",
               MAX_KEY_LENGTH, suffixlen);
        goto _deserialize_transition_error;
    }
    if (!(*read)(KEY, suffixlen, data))
        goto _deserialize_transition_error;
    KEY[suffixlen] = 0;
    if (!(transition->suffix = duplicate(KEY)))
        goto _deserialize_transition_error;
    if (!(*read)(&has_trie, sizeof(has_trie), data))
        goto _deserialize_transition_error;
    if (has_trie != 0 && has_trie != 1)
        goto _deserialize_transition_error;
    if (has_trie)
    {
        transition->next = Trie_new();
        if (!_deserialize_trie(transition->next, read, read_value, data))
            goto _deserialize_transition_error;
    }
    return 1;

_deserialize_transition_error:
    if (transition->suffix)
    {
        free(transition->suffix);
        transition->suffix = NULL;
    }
    if (transition->next)
    {
        Trie_del(transition->next);
        transition->next = NULL;
    }
    return 0;
}

Trie *Trie_deserialize(int (*read)(void *wasread, const int length, void *data),
                       void *(*read_value)(void *data),
                       void *data)
{
    Trie *trie = Trie_new();
    if (!_deserialize_trie(trie, read, read_value, data))
    {
        Trie_del(trie);
        return NULL;
    }
    return trie;
}


void _Trie_map_trie_data_helper(Trie* trie, void*(*func)(void*))
{
    if (trie->num_transitions == 0)
    {
        trie->value = func(trie->value);
    }
    else
    {
        for (int i = 0; i < trie->num_transitions; i++)
        {
            _Trie_map_trie_data_helper(trie->transitions[i].next, func);
        }
    }
}

int Trie_map_trie_data(Trie *trie, int(*func)(void*))
{
    _Trie_map_trie_data_helper(trie, func);
    return 1;
}

void _Trie_map_trie_helper(Trie* trie, void*(*func)(void*))
{
    if (trie->num_transitions == 0)
    {
        func(trie);
    }
    else
    {
        for (int i = 0; i < trie->num_transitions; i++)
        {
            _Trie_map_trie_helper(trie->transitions[i].next, func);
        }
    }
}

int Trie_map_trie(Trie *trie, int(*func)(void*))
{
    _Trie_map_trie_helper(trie, func);
    return 1;
}

void _Trie_map_parent_helper(Trie* trie, Transition *parent)
{
    trie->parent_transition = parent;
    if (trie->num_transitions != 0)
    {
        for (int i = 0; i < trie->num_transitions; i++)
        {
            trie->transitions[i].prev = trie;
            _Trie_map_parent_helper(trie->transitions[i].next, &trie->transitions[i]);
        }
    }
}

int Trie_map_parent(Trie* trie)
{
    _Trie_map_parent_helper(trie, NULL);
    return 1;
}

int Trie_size(const Trie *trie)
{
    return trie->size;
}

void* Trie_lazy_iter_next(Trie *trie)
{
    Trie* cur_trie = trie->cur_trie;
    char* cur_suffix = trie->cur_suffix;
    Transition* next_transition;    
    Transition* parent_transition = cur_trie->parent_transition;
    int i = 0;
    while (i < parent_transition->prev->num_transitions - 1) {
        if (parent_transition == &parent_transition->prev->transitions[i]) {
            next_transition = &parent_transition->prev->transitions[i+1];
            memset(cur_suffix + GRNA_L - strlen(parent_transition->suffix), 0, strlen(parent_transition->suffix));
            if (next_transition->next->setted) {
                trie->cur_trie = next_transition->next;
                memcpy(cur_suffix + GRNA_L - strlen(parent_transition->suffix), next_transition->suffix, strlen(next_transition->suffix));
                trie->cur_suffix = cur_suffix;
                PyTupleObject *py_tuple = PyTuple_New(3);
                PyTuple_SetItem(py_tuple, 0, PyLong_FromLong(trie->cur_trie->uid));
                PyTuple_SetItem(py_tuple, 1, PyUnicode_FromString(strcpyn(trie->cur_suffix, GRNA_L)));
                PyTuple_SetItem(py_tuple, 2, cur_trie->value);
                return py_tuple;
            } else {
                goto _traverse_down;
            }
        } else {
            i++;
        }
    }

    int cur_suffix_len = strlen(parent_transition->suffix);
    if (!(parent_transition = parent_transition->prev->parent_transition)) 
    {
        trie->cur_trie=NULL; 
        free(trie->cur_suffix);
        trie->cur_suffix = NULL;
        return NULL; // exhausted
    }
    memset(cur_suffix + GRNA_L - cur_suffix_len - strlen(parent_transition->suffix), 0, cur_suffix_len + strlen(parent_transition->suffix));

    i = 0;

    while (i < parent_transition->prev->num_transitions - 1) {
        if (parent_transition == &parent_transition->prev->transitions[i]) {
            next_transition = &parent_transition->prev->transitions[i+1];
            break;
        } else {
            i++;
            if (i == parent_transition->prev->num_transitions - 1) {
                if (!(parent_transition = parent_transition->prev->parent_transition)) 
                {
                    trie->cur_trie=NULL; 
                    free(trie->cur_suffix);
                    trie->cur_suffix = NULL;
                    return NULL; // exhausted
                }
                memset(cur_suffix + strlen(cur_suffix) - strlen(parent_transition->suffix), 0, strlen(parent_transition->suffix));
                i = 0;
            }
        }
    }
    goto _traverse_down;
    
_traverse_down:
    memcpy(cur_suffix+strlen(cur_suffix), next_transition->suffix, strlen(next_transition->suffix));
    while (!next_transition->next->setted)
    {
        next_transition = &next_transition->next->transitions[0];
        memcpy(cur_suffix+strlen(cur_suffix), next_transition->suffix, strlen(next_transition->suffix));
    }
    trie->cur_trie = next_transition->next;
    trie->cur_suffix = cur_suffix;
    PyTupleObject *py_tuple = PyTuple_New(3);
    PyTuple_SetItem(py_tuple, 0, PyLong_FromLong(trie->cur_trie->uid));
    PyTuple_SetItem(py_tuple, 1, PyUnicode_FromString(strcpyn(trie->cur_suffix, GRNA_L)));
    PyTuple_SetItem(py_tuple, 2, cur_trie->value);
    return py_tuple;
}

int Trie_lazy_iter_init(Trie *trie)
{
    trie->cur_trie = trie; // points to self initialization of the lazy iterator
    return 1;
}

void* Trie_lazy_iter(Trie* trie)
{
    Trie_map_parent(trie);
    char* cur_suffix = (char*) malloc (GRNA_L+1);
    memset(cur_suffix, 0, GRNA_L+1);
    if (!trie->cur_trie) {
        Trie_lazy_iter_init(trie);
        Transition *next_transition = trie->transitions;
        memcpy(cur_suffix, next_transition->suffix, strlen(next_transition->suffix));
        while (!next_transition->next->setted)
        {
            next_transition = &next_transition->next->transitions[0];  
            memcpy(cur_suffix+strlen(cur_suffix), next_transition->suffix, strlen(next_transition->suffix));
        } 
        Trie* cur_trie = next_transition->next;
        trie->cur_trie = cur_trie;
        trie->cur_suffix = cur_suffix;
        PyTupleObject *py_tuple = PyTuple_New(3);
        PyTuple_SetItem(py_tuple, 0, PyLong_FromLong(trie->cur_trie->uid));
        PyTuple_SetItem(py_tuple, 1, PyUnicode_FromString(strcpyn(trie->cur_suffix, GRNA_L)));
        PyTuple_SetItem(py_tuple, 2, cur_trie->value);
        return py_tuple;
    } else {
        Trie_lazy_iter_next(trie);
    }
}