#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

static PyObject *read_the_reference_tree(PyObject *self, PyObject *args);
static PyObject *assign_parents_and_children(PyObject *self, PyObject *args);
char *get_node_relationships(char *tree_string);
char *split_tree_string(char *tree_string);

static char module_docstring[] =
        "This module provides an interface for parsing Newick formatted trees using C from within MLTreeMap";
static char read_the_reference_tree_docstring[] =
        "Reads the labelled_tree_file and reformats it for downstream interpretation";
static char assign_parents_and_children_docstring[] = 
        "Stores the input tree as a binary search tree before recursively finding the children and parent of each node";

static PyMethodDef module_methods[] = {
        {"_read_the_reference_tree",
        read_the_reference_tree,
        METH_VARARGS,
        read_the_reference_tree_docstring},
        {"_assign_parents_and_children",
        assign_parents_and_children,
        METH_VARARGS,
        assign_parents_and_children_docstring},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_tree_parser(void) {
    PyObject *m = Py_InitModule3("_tree_parser", module_methods, module_docstring);
    if (m == NULL)
        return;
}

/**
  A tree node
  */
struct Node {
    int key;
    int left;
    int right;
    Node* next;
};
char const *float_chars = "0.123456789";
char const *real_number_chars = "-0123456789";

/*
 Returns 1 is character sub is a substring of super, 0 otherwise
 */
int is_char_substr(char sub, const char * super) {
    int x = 0;
    while(super[x]) {
        if (super[x] == sub)
            return 1;
        x++;
    }
    return 0;
}

/*
 Returns the length of the char_array as int
 */
int get_char_array_length(char * char_array) {
    int x = 0;
    while (char_array[x])
        x++;
    return x;
}

int append_char_array(int start, char * source, char *&dest) {
    int i = 0;
    while (source[i])
        dest[start++] = source[i++];
    return start;
}

char reverse_char_array(char * char_array, char *&flipped, int first, int last) {
    if (last == -1)
        return '\0';
    flipped[first] = reverse_char_array(char_array, flipped, first+1, last-1);
    return char_array[last];
}

// TODO: Fix memory leaks
static PyObject *read_the_reference_tree(PyObject *self, PyObject *args) {
    char* reference_tree_file;
    if (!PyArg_ParseTuple(args, "s", &reference_tree_file)) {
        return NULL;
    }

    FILE *reference_tree = fopen(reference_tree_file, "r");
    if (reference_tree == 0) {
        printf ("The reference tree file %s could not be opened for reading!\n", reference_tree_file);
        exit(0);
    }

    int _MAX = 1000;
    int tree_len = 0;
    int count = 2;
    int x = 0;
    char *tree_string = (char *) malloc ( _MAX * sizeof(char));
    char count_char[10];
    char c = fgetc(reference_tree);
    while ( c != EOF ) {
        if (tree_len <= _MAX && tree_len >= _MAX - 100) {
            _MAX = _MAX + 1000;
            tree_string = (char *) realloc(tree_string, _MAX * sizeof(char));
        }
        if (c == ')') {
            tree_string[tree_len] = c;
            tree_len++;
            tree_string[tree_len] = '-';
            tree_len++;

            sprintf(count_char, "%d", count);
            count++;
            x = 0;
            while (count_char[x]) {
                c = count_char[x];
                tree_string[tree_len] = c;
                tree_len++;
                x++;
            }
            c = fgetc(reference_tree);
        }
        if (c == ':') {
            c = fgetc(reference_tree);
            // while c is a substring of float_chars, continue reading in characters
            while (is_char_substr(c, float_chars) == 1)
                c = fgetc(reference_tree);
        }
        else {
            tree_string[tree_len] = c;
            tree_len++;
            c = fgetc(reference_tree);
        }

    }
    fclose(reference_tree);
    // Now remove the last node
    while ( c != ')') {
        tree_string[tree_len + 1] = '\0';
        tree_string[tree_len] = ';';
        tree_len--;
        c = tree_string[tree_len];
    }

    return Py_BuildValue("s", tree_string);
//    char * reference_tree_elements = (char *) malloc ( tree_len * sizeof(char));
//    reference_tree_elements = split_tree_string(tree_string);
//    return Py_BuildValue("s", reference_tree_elements);


};


char * split_tree_string(char *tree_string) {
    int _MAX = 1000;
    char * tree_elements = (char *) malloc ( _MAX * sizeof(char));
    int i = 0; // tree_string counter
    int j = 0; // element counter
    int k = 0; // tree_element position counter
    int x;
    char *element = (char*) malloc(10*sizeof(char));
    while (tree_string[i]) {
        if (k <= _MAX && k >= _MAX - 100) {
            _MAX = _MAX + 1000;
            tree_elements = (char *) realloc(tree_elements, _MAX * sizeof(char));
        }
        sprintf(element, "%d", j);
        x = 0;
        while (element[x]) {
            tree_elements[k++] = element[x++];
        }
        j++;
        tree_elements[k++] = '.';
        if (is_char_substr(tree_string[i], real_number_chars) == 0) {
            tree_elements[k] = tree_string[i];
        }
        else {
            while (is_char_substr(tree_string[i], real_number_chars) == 1) {
                tree_elements[k++] = tree_string[i++];
            }
            i--; k--;
        }
        i++; k++;
        tree_elements[k++] = '.';
    }
    tree_elements[--k] = '\0';

    return tree_elements;
}


char * get_previous_node(char *&parsed_tree_string, int &end) {
    char *previous = (char*) malloc(10*sizeof(char));
    char *reversed = (char*) malloc(10*sizeof(char));
    char c = parsed_tree_string[end];
    int i = 0;
    // Skip through brackets and commas to the end of the previous node
    while (is_char_substr(c, real_number_chars) == 0) {
        parsed_tree_string[end--] = '\0';
        c = parsed_tree_string[end];
    }
    while (is_char_substr(c, real_number_chars) == 1) {
        reversed[i] = c;
        parsed_tree_string[end--] = '\0';
        c = parsed_tree_string[end];
        i++;
    }
    reversed[i] = '\0';
    reverse_char_array(reversed, previous, 0, i);
    return previous;
}

/**
 * Inserts a new Node (with key=newKey) at the head of the linked_list.
 */
void insert(Node*& head, int newKey) {
    Node * curr = new Node;
    curr->key  = newKey;
    curr->next = head;
//    printf("key = %d\n", newKey);

    head = curr;
}


void print(Node* head) {
    std::cout << std::endl;
    for (Node* curr = head; curr != NULL; curr = curr->next){
        std::cout << "key: " << curr->key;
        std::cout << "\tleft: " << curr->left;
        std::cout << "\tright: " << curr->right;
        if (curr->next != NULL) std::cout << "\n";
    }
    std::cout << std::endl;
}


void load_linked_list(char * tree_string, char * parsed_tree_string, Node *&head) {
    char c;
    int pos = 0;
    int i = 0;
    int newKey = -1;
    int retrace_pos = 0;
    char *curr = (char*) malloc(10*sizeof(char));
    char *right = (char*) malloc(10*sizeof(char));
    char *left = (char*) malloc(10*sizeof(char));

    while (tree_string[pos]) {
        c = tree_string[pos];
        parsed_tree_string[retrace_pos++] = c;
        if (c == ')') {
            // load the next node as curr
            c = tree_string[pos+1];
            i = 0;
            while (is_char_substr(c, real_number_chars) == 1) {
                curr[i] = c;
                pos++;
                c = tree_string[pos+1];
                i++;
            }
            curr[i] = '\0';
            newKey = atoi(curr);
            if (newKey == 0)
                newKey = -1;
            insert(head, newKey);

            // load the previous 2 nodes as children and remove these from the string
            right = get_previous_node(parsed_tree_string, retrace_pos);
            left = get_previous_node(parsed_tree_string, retrace_pos);
            head->right = atoi(right);
            head->left = atoi(left);


            // add the current node to the parsed_tree_string
            i = 0;
            while (curr[i])
                parsed_tree_string[retrace_pos++] = curr[i++];
        }
        pos++;
    }
}


int get_children_of_nodes(Node * head, char *&children) {
    char * buffer = (char*) malloc (10);
    int _MAX = 1000;
    children = (char *) malloc (_MAX * sizeof(char));
    int x = 0;
    for (Node* curr = head; curr != NULL; curr = curr->next){
        if (x <= _MAX && x >= _MAX - 20) {
            _MAX = _MAX + 100;
            children = (char *) realloc (children, _MAX * sizeof(char));
        }
        sprintf(buffer, "%d", curr->key);
        x = append_char_array(x, buffer, children);
        children[x++] = '=';
        sprintf(buffer, "%d", curr->left);
        x = append_char_array(x, buffer, children);
        children[x++] = ',';
        sprintf(buffer, "%d", curr->right);
        x = append_char_array(x, buffer, children);

        if (curr->next != NULL) children[x++] = ';';
    }
    children[x++] = '\n';
    children[x] = '\0';

    return x;
}


int get_parents_of_nodes(Node * head, char *&parents) {
    char * parent = (char*) malloc (10);
    char * child = (char*) malloc (10);
    int _MAX = 1000;
    parents = (char *) malloc (_MAX * sizeof(char));
    int x = 0;
    for (Node* curr = head; curr != NULL; curr = curr->next){
        if (x <= _MAX && x >= _MAX - 20) {
            _MAX = _MAX + 100;
            parents = (char *) realloc (parents, _MAX * sizeof(char));
        }
        sprintf(parent, "%d", curr->key);

        // Load the left child
        sprintf(child, "%d", curr->left);
        x = append_char_array(x, child, parents);
        parents[x++] = ':';
        x = append_char_array(x, parent, parents);
        parents[x++] = ',';

        // Now the right child
        sprintf(child, "%d", curr->right);
        x = append_char_array(x, child, parents);
        parents[x++] = ':';
        x = append_char_array(x, parent, parents);

        if (curr->next != NULL) parents[x++] = ',';
    }
    parents[x] = '\n';

    return x;
}


int get_node_relationships(char *tree_string, char *&children, char *&parents) {
    /*
    :param tree_string: tree_info['subtree_of_node'][node]
    Function loads the whole tree into a tree struct where each node has a child and a parent
    Then the tree is queried for its ONE PARENT and potentially MULTIPLE CHILDREN
    */

    // Step 1: Load the tree
    Node * linked_list = NULL;
    int tree_len = get_char_array_length(tree_string);
    char * parsed_tree_string = (char*) malloc(tree_len * sizeof(char));

    load_linked_list(tree_string, parsed_tree_string, linked_list);
//    print(linked_list);

    //Step 2: Traverse the tree to get parents and children strings for each node
    int len_children = get_children_of_nodes(linked_list, children);
//    printf("Children:\n%s\n", children);

    int len_parents = get_parents_of_nodes(linked_list, parents);
//    printf("Parents:\n%s\n", parents);

    return len_children + len_parents;
}

static PyObject *assign_parents_and_children(PyObject *self, PyObject *args) {
    //TODO: Dynamically allocate more space if needed
    char* tree_string;
    if (!PyArg_ParseTuple(args, "s", &tree_string)) {
        return NULL;
    }

//    printf("%s\n", tree_string);
    char* children;
    char* parents;

    int length = get_node_relationships(tree_string, children, parents);

//    printf("Children:\n%s\n", children);
//    printf("Parents:\n%s\n", parents);

    int c_pos = 0;
    while (children[c_pos]){
        c_pos++;
    }

    children = (char *) realloc(children, (length + 10));

    int p_pos = 0;
    while (parents[p_pos])
        children[c_pos++] = parents[p_pos++];

    free(parents);

    return Py_BuildValue("s", children);
}
