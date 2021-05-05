#! /usr/bin/env python3

import os
import sys

from vcfdo.core import VcfIterator

print(sys.argv[1])
vcf = VcfIterator(sys.argv[1], nsites = 10000)
vcf.add_counter("mycount","A custom counter")
for site in vcf:
	print(site.CHROM, site.POS, vcf.ii, site.get_wsaf())
	#print(site)
