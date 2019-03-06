#pragma once
#include <Python.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_set>

using namespace std;

static PyObject *read_file(unordered_set<string> accession_list, const char * file);
unordered_set<string> listtoSet(PyObject* incoming);
PyObject* vectorToList_Str(const vector<string> &data);
static PyObject *parse_file(PyObject *module, PyObject *args);

static PyMethodDef parsers_methods[] = {
	{"parsefile", (PyCFunction)parse_file, METH_VARARGS, nullptr },
	{ NULL, NULL, 0, NULL }
};

static PyModuleDef parsers_module = {
	PyModuleDef_HEAD_INIT,
	"parsers",
	"File parser",
	0,
	parsers_methods
};

PyMODINIT_FUNC PyInit_parsers() {
	return PyModule_Create(&parsers_module);
}

