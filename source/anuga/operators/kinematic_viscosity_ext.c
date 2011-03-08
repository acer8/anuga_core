#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include "util_ext.h"

//Rough quicksort implementation (for build_operator_matrix)
int *quicksort(int *array, int n) {
	int *less, *more, pivot, i, num_less, num_more, *sorted;
	less = malloc(sizeof(int)*n);
	more = malloc(sizeof(int)*n);
	sorted = malloc(sizeof(int)*n);
	if (n<2) return array;
	pivot = array[n-1]; //choose the last element as a pivot
	num_less = 0; num_more = 0;
	for (i=0; i<n-1; i++) {
		if (array[i] <= pivot) {
			less[num_less] = array[i];
			num_less++;
		} else {
			more[num_more] = array[i];
			num_more++;
		}
	}
	
	less = quicksort(less,num_less);
	more = quicksort(more,num_more);
	
	for (i=0; i<n; i++) {
		if (i<num_less) {
			sorted[i] = less[i];
		} else if (i==num_less) {
			sorted[i] = pivot;
		} else {
			sorted[i] = more[i-num_less-1];
		}
	}
	return sorted;
}

int build_geo_structure(int n, 
                        int tot_len, 
                        double *centroids, 
                        long *neighbours, 
                        double *edgelengths, 
                        double *edge_midpoints, 
                        long *geo_indices, 
                        double *geo_values) {
    int i, edge, edge_counted, j, m;
	double dist, this_x, this_y, other_x, other_y, edge_length;
	edge_counted = 0;
	for (i=0; i<n; i++) {
		//The centroid coordinates of triangle i
		this_x = centroids[2*i];
		this_y = centroids[2*i+1];
		for (edge=0; edge<3; edge++) {
		
			j = neighbours[3*i+edge];
			
			//Get the index and the coordinates of the interacting point

                        // Edge
			if (j < 0 ) {
                                m = -j - 1;
				geo_indices[3*i+edge] = n + m;
				edge_counted++;
				
				other_x = edge_midpoints[2*(3*i+edge)];
				other_y = edge_midpoints[2*(3*i+edge)+1];
			} else {
				geo_indices[3*i+edge] = j;
				
				other_x = centroids[2*j];
				other_y = centroids[2*j+1];
			}
			
			//Compute the interaction
			edge_length = edgelengths[3*i+edge];
			dist = sqrt((this_x-other_x)*(this_x-other_x) + (this_y-other_y)*(this_y-other_y));
			geo_values[3*i+edge] = - edge_length / dist;
		}
	}
	return 0;
}

int build_operator_matrix(int n, 
                          int tot_len, 
                          long *geo_indices, 
                          double *geo_values, 
                          double *h, 
                          double *boundary_heights, 
                          double *data, 
                          long *colind) {
	int i, k, edge, j[4], *sorted_j, this_index;
	double h_j, v[3], v_i; //v[k] = value of the interaction of edge k in a given triangle, v_i = (i,i) entry
	for (i=0; i<n; i++) {
		v_i = 0.0;
		j[3] = i;
		//Get the values of each interaction, and the column index at which they occur
		for (edge=0; edge<3; edge++) {
			j[edge] = geo_indices[3*i+edge];
			if (j[edge]<n) { //interior
				h_j = h[j[edge]];
			} else { //boundary
				h_j = boundary_heights[j[edge]-n];
			}
			v[edge] = -0.5*(h[i] + h_j)*geo_values[3*i+edge]; //the negative of the individual interaction
			v_i += 0.5*(h[i] + h_j)*geo_values[3*i+edge]; //sum the three interactions
		}
		//Organise the set of 4 values/indices into the data and colind arrays
		sorted_j = quicksort((int *)j,4);
		for (k=0; k<4; k++) { //loop through the nonzero indices
			this_index = sorted_j[k];
			if (this_index == i) {
				data[4*i+k] = v_i;
				colind[4*i+k] = i;
			} else if (this_index == j[0]) {
				data[4*i+k] = v[0];
				colind[4*i+k] = j[0];
			} else if (this_index == j[1]) {
				data[4*i+k] = v[1];
				colind[4*i+k] = j[1];
			} else { //this_index == j[2]
				data[4*i+k] = v[2];
				colind[4*i+k] = j[2];
			}
		}
	}
	return 0;
}

/**
* Python wrapper methods
*/

static PyObject *py_build_geo_structure(PyObject *self, PyObject *args) {
    PyObject *kv_operator, *mesh;
    int n, tot_len, err;
    PyArrayObject
            *centroid_coordinates,
            *neighbours,
            *edgelengths,
            *edge_midpoint_coordinates,
			*geo_indices,
			*geo_values;
    
    //Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "Oii", &kv_operator, &n, &tot_len)) {
        PyErr_SetString(PyExc_RuntimeError, "build_geo_structure could not parse input");
        return NULL;
    }
    mesh = PyObject_GetAttrString(kv_operator, "mesh"); //kv_operator.mesh
    if (!mesh) {
        PyErr_SetString(PyExc_RuntimeError, "build_geo_structure could not obtain mesh object from kv_operator");
        return NULL;
    }

    //Extract parameters
    centroid_coordinates = get_consecutive_array(mesh,"centroid_coordinates");
    neighbours = get_consecutive_array(mesh,"neighbours");
    edgelengths = get_consecutive_array(mesh,"edgelengths");
    edge_midpoint_coordinates = get_consecutive_array(mesh,"edge_midpoint_coordinates");
    geo_indices = get_consecutive_array(kv_operator,"geo_structure_indices");
    geo_values = get_consecutive_array(kv_operator,"geo_structure_values");
    
    //Release
    Py_DECREF(mesh);
    
    err = build_geo_structure(n,tot_len, 
                             (double *)centroid_coordinates -> data, 
                             (long *)neighbours -> data, 
                             (double *)edgelengths->data, 
                             (double *)edge_midpoint_coordinates -> data, 
                             (long *)geo_indices -> data, 
                             (double *)geo_values -> data);
    if (err != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Could not build geo structure");
        return NULL;
    }
    
    //Release the arrays
    Py_DECREF(centroid_coordinates);
    Py_DECREF(neighbours);
    Py_DECREF(edgelengths);
    Py_DECREF(edge_midpoint_coordinates);
    Py_DECREF(geo_indices);
    Py_DECREF(geo_values);
    
    return Py_BuildValue("");
}

static PyObject *py_build_operator_matrix(PyObject *self, PyObject *args) {
	PyObject *kv_operator;
	int n, tot_len, err;
	PyArrayObject
		*h,
		*boundary_heights,
		*geo_indices,
		*geo_values,
		*_data,
		*colind;
	
	//Convert Python arguments to C
    if (!PyArg_ParseTuple(args, "OiiOO", &kv_operator, &n, &tot_len, &h, &boundary_heights)) {
        PyErr_SetString(PyExc_RuntimeError, "get_stage_height_interactions could not parse input");
        return NULL;
    }
	
	geo_indices = get_consecutive_array(kv_operator,"geo_structure_indices");
	geo_values = get_consecutive_array(kv_operator,"geo_structure_values");
	_data = get_consecutive_array(kv_operator,"operator_data");
	colind = get_consecutive_array(kv_operator,"operator_colind");
	
	err = build_operator_matrix(n,tot_len, 
								(long *)geo_indices -> data, 
								(double *)geo_values -> data, 
								(double *)h -> data, 
								(double *)boundary_heights -> data, 
								(double *)_data -> data, 
								(long *)colind -> data);
    if (err != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Could not get stage height interactions");
        return NULL;
    }
	
	Py_DECREF(geo_indices);
	Py_DECREF(geo_values);
	Py_DECREF(_data);
	Py_DECREF(colind);
    
    return Py_BuildValue("");
}

static struct PyMethodDef MethodTable[] = {
        {"build_geo_structure",py_build_geo_structure,METH_VARARGS,"Print out"},
		{"build_operator_matrix",py_build_operator_matrix,METH_VARARGS,"Print out"},
        {NULL,NULL,0,NULL} // sentinel
};
void initkinematic_viscosity_ext(){
    (void) Py_InitModule("kinematic_viscosity_ext", MethodTable);
    import_array();
}
