#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>

#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <map>

#include <Python.h>
#include <sstream>
#include <unordered_map>          

#include <time.h>

using namespace std;

typedef long long int64; typedef unsigned long long uint64;

uint64 getTimeMs64()
{
  #ifdef _WIN32
  /* Windows */
  FILETIME ft;
  LARGE_INTEGER li;

  /* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
   * to a LARGE_INTEGER structure. */
  GetSystemTimeAsFileTime(&ft);
  li.LowPart = ft.dwLowDateTime;
  li.HighPart = ft.dwHighDateTime;

  uint64 ret = li.QuadPart;
  ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
  ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

  return ret;
  #else
  /* Linux */
  struct timeval tv;

  gettimeofday(&tv, NULL);

  uint64 ret = tv.tv_usec;
  /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
  ret /= 1000;

  /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
  ret += (tv.tv_sec * 1000);

  return ret;
  #endif
}


static PyObject *parse_file(PyObject *module, PyObject *args);
PyObject *read_file(unordered_map<string, bool> accession_list, const char * file);
unordered_map<string, bool> listtoSet(PyObject* incoming);
PyObject* vectorToList_Str(const vector<string> &data);

static PyMethodDef module_methods[] = {
        {"parse_file",
        parse_file,
        METH_VARARGS,
        NULL},
        {NULL, NULL, 0, NULL}
};

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#if PY_MAJOR_VERSION >= 3

    static int module_traverse(PyObject *m, visitproc visit, void *arg) {
        Py_VISIT(GETSTATE(m)->error);
        return 0;
    }

    static int module_clear(PyObject *m) {
        Py_CLEAR(GETSTATE(m)->error);
        return 0;
    }

    static struct PyModuleDef module_def = {
        PyModuleDef_HEAD_INIT,
        "_parser_factory",
        "This module parses tax ids from fileP",
        sizeof(struct module_state),
        module_methods,    /* m_methods */
        NULL,                /* m_reload */
        module_traverse,                /* m_traverse */
        module_clear,                /* m_clear */
        NULL,                /* m_free */
    };

#define INITERROR return NULL

PyMODINIT_FUNC PyInit__parser_factory(void)

#else
#define INITERROR return

PyMODINIT_FUNC init_parser_factory(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *m = PyModule_Create(&module_def);
#else
    static char module_docstring[] =
        "This module parses tax ids.";
    PyObject *m = Py_InitModule3("_parser_factory", module_methods, module_docstring);
#endif

    if (m == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(m);

    st->error = PyErr_NewException("_parser_factory.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(m);
        INITERROR;
    }

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}

