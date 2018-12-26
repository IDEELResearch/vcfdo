#! /usr/bin/env python3

from setuptools import setup

setup(
	name = "vcfdo",
	version = "0.1",
	description = "Utlity functions for processing and annotating VCF files, focused on Plasmodium spp.",
	url = "http://github.com/IDEELResearch/vcfdo",
	author = "Andrew Parker Morgan",
	author_email='andrew.parker.morgan@gmail.com',
	license= "MIT",
	packages = ["vcfdo"],
	scripts=["bin/vcfdo"],
	install_requires = [
		"cyvcf2",
		"scikit-allel",
		"pyfaidx"
	]
	zip_safe = False
)
