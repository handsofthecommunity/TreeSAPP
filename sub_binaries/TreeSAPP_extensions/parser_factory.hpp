#include <Python.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_set>

using namespace std;

static PyObject *parse_file(PyObject *module, PyObject *args);
PyObject *read_file(unordered_set<string> accession_list, const char * file);
unordered_set<string> listtoSet(PyObject* incoming);
PyObject* vectorToList_Str(const vector<string> &data);

static PyMethodDef parsers_methods[] = {
	{"parse_file", (PyCFunction)parse_file, METH_VARARGS, NULL },
	{ NULL, NULL, 0, NULL }
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
  "parser_factory",
  "File parser",
  0,
  parsers_methods,
  NULL,
  module_traverse,
  module_clear,
  NULL,
};

#define INITERROR return NULL

PyMODINIT_FUNC PyInit_parser_factory(void){
     return PyModule_Create(&module_def);
}
#endif
