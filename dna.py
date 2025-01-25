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









