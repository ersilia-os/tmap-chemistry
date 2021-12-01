from rdkit.Chem.Descriptors import MolWt
import pandas as pd


class DefaultDescriptors(object):

	def __init__(self, input_file, output_file):
		self.input_file = input_file
		self.output_file = output_file

	def run(self):
		pass