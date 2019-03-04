#pragma once
#include <Python.h>

#include <Windows.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

static PyObject *read_file(vector<string> accession_list, char * file);
vector<string> listToVector_Str(PyObject* incoming);
PyObject* vectorToList_Str(const vector<string> &data);
static PyObject *parse_file(PyObject *module, PyObject *args);
PyObject* tanh_impl(PyObject *, PyObject* o);

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

