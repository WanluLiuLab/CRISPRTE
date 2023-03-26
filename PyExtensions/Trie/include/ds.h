#include <stdlib.h>

#define VARLST_BS_L 16
#define VARLST_SEQ_L 23

typedef struct varlst
{
    char** _lst;
    size_t _sz;
    size_t _cur;
    size_t _bytes;
} Varlst;

Varlst* varlst_init();

int varlst_is_empty(Varlst* lst);

int varlst_insert(Varlst* lst, char* s);

int varlst_expand(Varlst* lst, size_t new_sz);

