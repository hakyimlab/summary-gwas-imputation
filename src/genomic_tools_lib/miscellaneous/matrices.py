__author__ = "alvaro barbeira"
import numpy
import math
from .. import  Utilities

def _flatten_matrix_data(data):
    """data is expected to be a list of (name, id_labels, matrix) tuples"""
    results = []
    for name, id_labels, matrix in data:
        if len(id_labels) == 1:
            id = id_labels[0]
            results.append((name, id, id, float(matrix)))
            continue

        for i in range(0, len(id_labels)):
            for j in range(i, len(id_labels)):
                value = matrix[i][j]
                id1 = id_labels[i]
                id2 = id_labels[j]
                results.append((name, id1, id2, str(value)))
    return results

def _flatten_matrix_data_2(id_labels, matrix):
    """data is expected to be a list of (name, id_labels, matrix) tuples, but names will be unused"""
    id1=[]
    id2=[]
    value=[]
    for i in range(0, len(id_labels)):
        for j in range(i, len(id_labels)):
            value.append(float(matrix[i,j]))
            id1.append(id_labels[i])
            id2.append(id_labels[j])

    return {"id1":id1, "id2":id2, "value":value}

def _flatten_matrix_data_3(id_labels, id_key, matrix):
    """data is expected to be a list of (name, id_labels, matrix) tuples, but names will be unused"""
    id1=[]
    id2=[]
    value=[]
    for i in range(0, len(id_labels)):
        for j in range(i, len(id_labels)):
            value.append(float(matrix[i,j]))
            id1.append(id_key[id_labels[i]])
            id2.append(id_key[id_labels[j]])

    return {"id1":id1, "id2":id2, "value":value}

def matrix_data_to_dataframe(data, columns=["name", "id1", "id2", "value"]):
    f = _flatten_matrix_data([data])
    f = Utilities.to_dataframe(f, columns=columns)
    return f

def matrices_data_to_dataframe(data, columns=["name", "id1", "id2", "value"]):
    f = _flatten_matrix_data(data)
    f = Utilities.to_dataframe(f, columns=columns)
    return f

def matrix_data(name, labels, matrix):
    return (name, labels, matrix)