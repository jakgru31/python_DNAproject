#!/usr/bin/env python3
# -*- coding: utf-8 -*-


class FastaDNA:
    def __init__(self, file):
        self.file = file

    def _get_identifier(self):
        with open(self.file, 'r') as f:
            return f.readline().strip()

    def _get_sequence(self):
        with open(self.file, 'r') as f:
            return f.read().splitlines()[1:]

    def _get_sequence_string(self):
        return ''.join(self._get_sequence())


class DNA(FastaDNA):
    def __init__(self, file):
        super().__init__(file)
        self.identifier = self._get_identifier()
        self.sequence = self._get_sequence_string()

    def __str__(self):
        return f">{self.identifier}\n{self.sequence}"

    def __len__(self):
        return len(self.sequence)

    def get_sequence(self):
        return self.sequence

    def get_protein_sequence(self):
        protein_dict = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }

        protein_sequence = []
        for i in range(0, len(self.sequence), 3):
            codon = self.sequence[i:i + 3]
            if len(codon) == 3:
                protein_sequence.append(protein_dict.get(codon, ''))
        return ''.join(protein_sequence)





