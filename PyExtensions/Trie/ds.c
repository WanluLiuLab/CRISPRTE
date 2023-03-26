#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "include/ds.h"

struct Varlst
{
    char** _lst;
    size_t _sz;
    size_t _cur;
    size_t _bytes;
};

Varlst* varlst_init()
{   
    Varlst* new_lst;
    if (!(new_lst = malloc(sizeof(Varlst)))) return (Varlst*) NULL;
    new_lst->_sz = VARLST_BS_L;
    new_lst->_cur = 0;
    char** l;
    if (!(l =  malloc(sizeof(char*)*VARLST_BS_L))) return (Varlst*) NULL;
    new_lst->_lst = l;
    return new_lst;
}


int varlst_is_empty(Varlst* lst)
{
    return lst->_cur == 0;
}

int varlst_insert(Varlst* lst, char* s)
{
    if (!(lst->_lst[lst->_cur] = (char*)malloc(sizeof(char)*strlen(s)+1))) return 0;
    lst->_bytes += sizeof(char)*strlen(s)+1;
    memset(lst->_lst[lst->_cur], 0, strlen(s));
    memcpy(lst->_lst[lst->_cur], s, strlen(s)+1);
    lst->_cur++;
    if (lst->_cur >= lst->_sz) 
    { if (!varlst_expand(lst, lst->_sz * 2)) return 0;}
    return 1;
}

int varlst_expand(Varlst* lst, size_t new_sz)
{
    char** new_lst;
    if (!(new_lst = malloc(lst->_bytes))) return 0;
    memcpy(new_lst, lst->_lst, lst->_bytes);
    free(lst->_lst);
    lst->_lst = new_lst;
    lst->_sz = new_sz;
    return 1;
}

int varlst_free(Varlst* lst)
{
    free(lst->_lst);
    free(lst);
    return 0;
}