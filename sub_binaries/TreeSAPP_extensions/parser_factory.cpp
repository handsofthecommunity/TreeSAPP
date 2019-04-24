#include "parser_factory.hpp"

vector<StringRef> split3(string const& str, bool check, char delimiter = '\t'){
  vector<StringRef> result;
  enum State {inSpace, inToken};

  State state = inSpace;
  char const* pTokenBegin = 0;
  for (auto it = str.begin(); it != str.end(); ++it){
    State const newState = (*it == delimiter? inSpace : inToken);

    if (newState != state){
      switch (newState){
      case inSpace:
	result.push_back(StringRef(pTokenBegin, &*it - pTokenBegin));
	if (check){
	  return result;
	} 
	break;
      case inToken:
	pTokenBegin = &*it;
      }
    }
    state = newState;
    
  }
  if (state == inToken){
    result.push_back(StringRef(pTokenBegin, &*str.end() - pTokenBegin));
  }
  return result;
}

PyObject *read_file(unordered_map<string, bool> accession_list, const char * file, int listSize) {
  /*                                                                                                         
   * Function for reading file with accession and tax id data                                                
   * to a list of accession ids, tax ids, and versions                                                       
   */
  try {
    std::ios_base::sync_with_stdio(false);
    cin.sync_with_stdio(false);
    // disable async IO
    string line;
    ifstream acc_file(file);
    vector<StringRef> data;
    data.reserve(3);
    vector<string> accession_ids = {};
    accession_ids.reserve(listSize);
    vector<string> versions = {};
    versions.reserve(listSize);
    vector<string> taxid = {};
    taxid.reserve(listSize);

    vector<vector<string>> allLists = {};
    allLists.reserve(3);
    
    int count = 0;
    int check = 0;
    int size = listSize;

    while (getline(acc_file, line)) {
      std::istringstream iss(line);
    
      data = split3(line, true);

      string accession_id = string(data[0].begin(), data[0].end());

      if (accession_list.find(accession_id) == accession_list.end() || accession_list[accession_id]){
	check = 1;
      }

      if (check == 0) {
	data = split3(line, false);
	string version = string(data[1].begin(), data[1].end());
	string tax_id = string(data[2].begin(), data[2].end());

      	count++;
      	accession_list[accession_id] = true;
      	accession_ids.push_back(accession_id);
      	versions.push_back(version);
      	taxid.push_back(tax_id);
      }
      
      if (count >= size) {
      	break;
      }
      check = 0;
      data.clear();
      
    }

    acc_file.close();

    allLists.push_back(accession_ids);
    allLists.push_back(versions);
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

  PyObject *retLists = read_file(accession_list, fileName, accession_list.size());

  return retLists;
}

