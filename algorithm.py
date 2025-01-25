#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import dna
import pandas as pd


def damerau_levenshtein(str1, str2):
    """
    Computes the Damerau-Levenshtein distance between two strings
    :param str1: first string
    :param str2: second string
    :return: Damerau-Levenshtein distance matrix
    """
    if not isinstance(str1, str) or not isinstance(str2, str):
        raise ValueError("Invalid input")
    if len(str1) == 0 or len(str2) == 0:
        raise ValueError("Empty strings")

    d1, d2 = len(str1), len(str2)
    nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    data_array = np.zeros(4, dtype=int)
    maximal_distance = d1 + d2
    distance_matrix = np.zeros((d1 + 2, d2 + 2), dtype=int)

    distance_matrix[0, 0] = maximal_distance
    for i in range(1, d1 + 2):
        distance_matrix[i, 0] = maximal_distance
        distance_matrix[i, 1] = i - 1
    for j in range(1, d2 + 2):
        distance_matrix[0, j] = maximal_distance
        distance_matrix[1, j] = j - 1

    for i in range(1, d1+1):
        db = 0
        for j in range(1, d2+1):
            k = data_array[nucleotides[str2[j-1]]]
            l = db
            if str1[i-1] == str2[j-1]:
                cost = 0
                db = j
            else:
                cost = 1

            distance_matrix[i+1, j+1] = min(
                distance_matrix[i, j] + cost,
                distance_matrix[i+1, j] + 1,
                distance_matrix[i, j+1] + 1,
                distance_matrix[k, l] + (i-k-1) + 1 + (j-l-1)
            )
        data_array[nucleotides[str1[i-1]]] = i
    df = pd.DataFrame(
        distance_matrix,
        index=[""] + list(str1) + [" "],  # Ensure the index matches the matrix size
        columns=[""] + list(str2) + [" "]  # Ensure the columns match the matrix size
    )

    print("Damerau-Levenshtein Distance Matrix:")
    print(df.to_string(index=True, header=True))

    return distance_matrix[d1+1, d2+1]


if __name__ == '__main__':
    dna1 = dna.DNA("gene_lib/ABL1_25.fasta")
    dna2 = dna.DNA("gene_lib/AR_231.fasta")
    damerau_levenshtein(dna1.get_sequence(), dna2.get_sequence())