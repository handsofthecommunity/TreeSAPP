#include "parser_factory.hpp"


PyObject *read_file(unordered_map<string, bool> accession_list, const char * file) {
  /*                                                                                                         
   * Function for reading file with accession and tax id data                                                
   * to a list of accession ids, tax ids, and versions                                                       
   */
  try {
    string line;
    ifstream acc_file(file);
    string tmp;
    vector<string> data;

    vector<string> accession_ids = {};
    vector<string> version = {};
    vector<string> taxid = {};

    vector<vector<string>> allLists = {};

    int count = 0;
    int check = 0;
    int size = accession_list.size();
    
    clock_t start = clock();
    clock_t middle;
    clock_t end;

    while (getline(acc_file, line)) {
      std::istringstream iss(line);
      while (getline(iss, tmp, '\t')) {
      	data.push_back(tmp);
      	if (accession_list.find(data[0]) == accession_list.end() || accession_list[data[0]]) {
      	  check = 1;
      	  break;
      	}
      }
      
      if (accession_list.find(data[0]) != accession_list.end() && check == 0) {
      	count++;
      	accession_list[data[0]] = true;
      	accession_ids.push_back(data[0]);
      	version.push_back(data[1]);
      	taxid.push_back(data[2]);
      }
      
      if (count >= size) {
      	break;
      }
      check = 0;
      data.clear();
    }

    cout << ((float_t) clock() - start)/CLOCKS_PER_SEC << endl;
    acc_file.close();
    allLists.push_back(accession_ids);
    allLists.push_back(version);
    allLists.push_back(taxid);

    PyObject* retList = PyList_New(3);
    if (!retList) throw logic_error("Unable to allocate enough memory for output...");
    for (unsigned int i = 0; i < 3; i++) {
      PyObject *list = vectorToList_Str(allLists[i]);
      if (!list) {
	Py_DECREF(retList);
	throw logic_error("Unable to allocate memory...");
      }
      PyList_SetItem(retList, i, list);
    }

    return retList;
  }
  catch (std::exception const& e) {
    cout << "Error reading accession file..." << e.what() << endl;
  }

  return NULL;
}

unordered_map<string, bool> listtoSet(PyObject* incoming) {
  unordered_map<string, bool> accession_list;
  int numLines = PyList_Size(incoming);
  for (Py_ssize_t i = 0; i < numLines; i++) {
    PyObject *strObj = PyList_GetItem(incoming, i);
    PyObject *str = PyUnicode_AsEncodedString(strObj, "utf-8", "~E~");
    string accession_id = PyBytes_AsString(str);
    accession_list[accession_id] = false;
  }
  return accession_list;
}

PyObject* vectorToList_Str(const vector<string> &data) {
  PyObject* listObj = PyList_New(data.size());
  if (!listObj) throw logic_error("Unable to allocate enough memory for output...");
  for (unsigned int i = 0; i < data.size(); i++) {
    PyObject *str = PyUnicode_FromString((const char *) data[i].c_str());
    if (!str) {
      Py_DECREF(listObj);
      throw logic_error("Unable to allocate memory...");
    }
    PyList_SetItem(listObj, i, str);
  }
  return listObj;
}

static PyObject* parse_file(PyObject *module, PyObject* args) {
  PyObject *accList;
  const char * fileName;
  if (!PyArg_ParseTuple(args, "sO!", &fileName, &PyList_Type, &accList)) {
    fprintf(stderr, "ERROR: Argument types are incorrect...\n");
    return NULL;
  }
  unordered_map<string, bool> accession_list = listtoSet(accList);

  PyObject *retLists = read_file(accession_list, fileName);

  return retLists;
}

