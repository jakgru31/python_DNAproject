#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import dna
import matplotlib.pyplot as plt
import pandas as pd


class GeneAnalysis:
    def __init__(self, data1, data2):
        self.data1 = data1
        self.data2 = data2
        if not isinstance(data2, dna.DNA):
            raise ValueError("Invalid input")


    def damerau_levenshtein_dna(self):
        """
        A method to compute the Damerau-Levenshtein distance between two DNA sequences
        :return: df (DataFrame): Damerau-Levenshtein distance matrix
        """
        if not isinstance(self.data1, dna.DNA) or not isinstance(self.data2, dna.DNA):
            raise ValueError("Invalid input")

        seq1 = self.data1.get_sequence()
        seq2 = self.data2.get_sequence()

        if len(seq1) == 0 or len(seq2) == 0:
            raise ValueError("Empty strings")

        d1, d2 = len(seq1), len(seq2)
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
                k = data_array[nucleotides[seq2[j-1]]]
                l = db
                if seq1[i-1] == seq2[j-1]:
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
            data_array[nucleotides[seq1[i-1]]] = i
        df = pd.DataFrame(
            distance_matrix,
            index=[""] + list(seq1) + [" "],  # Ensure the index matches the matrix size
            columns=[""] + list(seq2) + [" "]  # Ensure the columns match the matrix size
        )
        return df

    def damerau_levenshtein_protein(self):
        """

        :return:
        """
        if not isinstance(self.data1, dna.DNA) or not isinstance(self.data2, dna.DNA):
            raise ValueError("Invalid input")

        seq1 = self.data1.get_protein_sequence()
        seq2 = self.data2.get_protein_sequence()

        if len(seq1) == 0 or len(seq2) == 0:
            raise ValueError("Empty strings")

        d1, d2 = len(seq1), len(seq2)
        amino_acids = {
            'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
            'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, '*': 20
        }
        data_array = np.zeros(21, dtype=int)
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
                k = data_array[amino_acids[seq2[j-1]]]
                l = db
                if seq1[i-1] == seq2[j-1]:
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
            data_array[amino_acids[seq1[i-1]]] = i
        df = pd.DataFrame(
            distance_matrix,
            index=[""] + list(seq1) + [" "],  # Ensure the index matches the matrix size
            columns=[""] + list(seq2) + [" "]  # Ensure the columns match the matrix size
        )
        return df

    def _shortest_path(self, distance_matrix):
        rows, cols = distance_matrix.shape
        path = [(rows - 1, cols - 1)]  # Start at the bottom-right corner
        r, c = rows - 1, cols - 1
        while r > 0 or c > 0:
            moves = []
        if r > 0:  # Move up
            moves.append((r - 1, c))
        if c > 0:  # Move left
            moves.append((r, c - 1))
        if r > 0 and c > 0:  # Move diagonally
            moves.append((r - 1, c - 1))

        # Choose the move with the minimum cost
        r, c = min(moves, key=lambda x: distance_matrix.iloc[x])
        path.append((r, c))

        return path[::-1]  # Reverse the path to start from the top-left

    def visualize_matrix(self, param):
        if param == 'dna':
            seq1 = self.data1.get_sequence()
            seq2 = self.data2.get_sequence()
            fig, ax = plt.subplots(figsize=(10, 10))
            heatmap = ax.imshow(self.damerau_levenshtein_dna(), cmap="coolwarm", aspect="equal")
            plt.colorbar(heatmap, label="Score")
            path = self._shortest_path(self.damerau_levenshtein_dna())
            path_x, path_y = zip(*path)
            ax.plot(path_y, path_x, color="red", linewidth=2, label="Optimal Path")
            ax.scatter(path_y, path_x, color="red", s=30)  # Mark path points
            for i in range(self.damerau_levenshtein_dna().shape[0]):
                for j in range(self.damerau_levenshtein_dna().shape[1]):
                    score = self.damerau_levenshtein_dna().iloc[i, j]
                    ax.text(j, i, f"{score}", ha="center", va="center", color="black", fontsize=8)
                    if (i, j) in path:
                        if (i-1, j-1) in path:
                            ax.arrow(j-0.3, i-0.3, -0.4, -0.4, head_width=0.1, head_length=0.1, fc='black', ec='black')
                        elif (i-1, j) in path:
                            ax.arrow(j, i-0.4, 0, -0.4, head_width=0.1, head_length=0.1, fc='black', ec='black')
                        elif (i, j-1) in path:
                            ax.arrow(j-0.4, i, -0.4, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')

            ax.set_xticks(range(1, len(seq2) + 1))
            ax.set_xticklabels(list(seq2), fontsize=10)
            ax.set_yticks(range(1, len(seq1) + 1))
            ax.set_yticklabels(list(seq1), fontsize=10)
            ax.set_xlabel("Query DNA Sequence", fontsize=12)
            ax.set_ylabel("Reference DNA Sequence", fontsize=12)
            ax.set_title("Alignment Matrix", fontsize=14)
            plt.legend(loc="upper left")
            plt.show()




        elif param == 'protein':
            seq1 = self.data1.get_protein_sequence()
            seq2 = self.data2.get_protein_sequence()
            fig, ax = plt.subplots(figsize=(10, 10))
            heatmap = ax.imshow(self.damerau_levenshtein_protein(), cmap="coolwarm", aspect="equal")
            plt.colorbar(heatmap, label="Score")
            path = self._shortest_path(self.damerau_levenshtein_protein())
            path_x, path_y = zip(*path)
            ax.plot(path_y, path_x, color="red", linewidth=2, label="Optimal Path")
            ax.scatter(path_y, path_x, color="red", s=30)  # Mark path points
            for i in range(self.damerau_levenshtein_protein().shape[0]):
                for j in range(self.damerau_levenshtein_protein().shape[1]):
                    score = self.damerau_levenshtein_protein().iloc[i, j]
                    ax.text(j, i, f"{score}", ha="center", va="center", color="black", fontsize=8)
                    if (i, j) in path:
                        if (i-1, j-1) in path:
                            ax.arrow(j-0.3, i-0.3, -0.4, -0.4, head_width=0.1, head_length=0.1, fc='black', ec='black')
                        elif (i-1, j) in path:
                            ax.arrow(j, i-0.4, 0, -0.4, head_width=0.1, head_length=0.1, fc='black', ec='black')
                        elif (i, j-1) in path:
                            ax.arrow(j-0.4, i, -0.4, 0, head_width=0.1, head_length=0.1, fc='black', ec='black')

            ax.set_xticks(range(1, len(seq2) + 1))
            ax.set_xticklabels(list(seq2), fontsize=10)
            ax.set_yticks(range(1, len(seq1) + 1))
            ax.set_yticklabels(list(seq1), fontsize=10)
            ax.set_xlabel("Query Protein Sequence", fontsize=12)
            ax.set_ylabel("Reference Protein Sequence", fontsize=12)
            ax.set_title("Alignment Matrix", fontsize=14)
            plt.legend(loc="upper left")
            plt.show()

        else:
            raise ValueError(f"{param} parameter is not valid")




if __name__ == '__main__':
    dna1 = dna.DNA("gene_lib/ATR_545.fasta")
    dna2 = dna.DNA("gene_lib/ATR_390502.fasta")
    analysis1 = GeneAnalysis(dna1, dna2)
    analysis1.damerau_levenshtein_protein()
    print(analysis1.damerau_levenshtein_protein())
