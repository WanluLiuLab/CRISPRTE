#include <Python.h>
#include <marshal.h>

/* forward declaration */
typedef struct Gloc Gloc;

/* convert Gloc obeject to python object. Deprecated soon. */
PyObject* _gloc_c2pytuple_helper(Gloc *mp);

/* convert Gloc obeject to python object. Deprecated soon. */
PyObject* gloc_c2pytuple(Gloc *mp);