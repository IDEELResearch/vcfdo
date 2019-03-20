#! /usr/bin/env python3

import os
import sys
import numpy as np
import logging

from re import split as resplit
from collections import defaultdict
from datetime import datetime as dt
from cyvcf2 import VCF, Variant

class SiteCounter:

	def __init__(self, name = "counter", description = "I'm a counter", value = 0):
		self.name = name
		self.description = description
		self.count = value

	def __iadd__(self, value):
		self.count += value

	def reset(self):
		self.count = 0

	def increment(self, by = 1):
		self.count += by


class VcfIterator(VCF):

	def __init__(self, fname = "/dev/stdin", threads = None, log = None, nsites = None, log_name = "vcfdo", breadcrumb = None, **kwargs):

		if isinstance(log, logging.Logger):
			self.log = log
		else:
			logging.basicConfig(level = logging.INFO)
			logging.StreamHandler(stream = sys.stderr)
			self.log = logging.getLogger(log_name)
		self._fname = fname
		self.ploidy = None
		self.assume_ref_ancestral = False

		self._ns = nsites if nsites else -1
		self.ii = 0
		self.counters = {}
		self.counter_names = []
		self._has_standard_counters = False

		self.log.info("Connecting to VCF file <{}>".format(fname))
		super().__init__(fname, gts012 = True, threads = threads, **kwargs)
		self.nsamples = len(self.samples)
		if breadcrumb is not None:
			timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")
			self.add_to_header("##{}; timestamp={}".format(breadcrumb, timestamp))

	def set_ploidy(self, site):
		if self.ploidy is None:
			self.ploidy = int(site.ploidy)
		elif self.ploidy != site.ploidy:
			self.log.warning("Warning: possible ploidy conflict: was {}, now {} ".format(self.ploidy, site.ploidy))
			self.ploidy = int(site.ploidy)

	def add_counter(self, name, description, value = 0):
		self.counter_names.append(name)
		self.counters[name] = SiteCounter(name, description, value)

	def add_standard_counter(self, name):
		self._has_standard_counters = True
		if name == "multiallelic":
			self.add_counter("nmulti","sites had >2 alleles")
		elif name == "nowsaf":
			self.add_counter("nnowsaf","sites had no WSAF annotations")
		elif name == "noancestral":
			self.add_counter("nnoanc","sites had undefined ancestral allele")
		elif name == "nonpolymorphic":
			self.add_counter("nnonpoly","sites not polymorphic in called genotypes")
		else:
			self._has_standard_counters = False

	def __next__(self):

		self.ii += 1
		if (self.ii-1) and not (self.ii-1) % 1000:
			self.log.info("\t... processed {} sites ...".format(self.ii-1))

		if self.ii > self._ns and self._ns > -1:
			self.ii -= 1
			#self.show_count_summary()
			raise StopIteration
		else:
			try:
				next_site = VcfSite( super().__next__(), self.assume_ref_ancestral )
				self.set_ploidy(next_site)
				if self._has_standard_counters:
					if "multiallelic" in self.counter_names:
						if next_site.is_multiallelic:
							self.counters["nmulti"].increment()
						else:
							pass
					if "nowsaf" in self.counter_names:
						if not next_site.has_wsaf:
							self.counters["nnowsaf"].increment()
						else:
							pass
					if "noancestral" in self.counter_names:
						if not next_site.has_ancestral:
							self.counters["nnoanc"].increment()
						else:
							pass
					if "nonpolymorphic" in self.counter_names:
						if not next_site.is_polymorphic:
							self.counters["nnonpoly"].increment()
						else:
							pass
				return next_site
			except StopIteration:
				#self.show_count_summary()
				self.ii -= 1
				raise StopIteration

	def fetch_region(self, region = None):

		site_iter = self.__call__(region)
		while True:

			self.ii += 1
			if (self.ii-1) and not (self.ii-1) % 1000:
				self.log.info("\t... processed {} sites ...".format(self.ii-1))

			if self.ii > self._ns and self._ns > -1:
				self.ii -= 1
				#self.show_count_summary()
				raise StopIteration

			else:
				try:
					next_site = VcfSite( next(site_iter), self.assume_ref_ancestral )
					self.set_ploidy(next_site)
					if self._has_standard_counters:
						if "multiallelic" in self.counter_names:
							if next_site.is_multiallelic:
								self.counters["nmulti"].increment()
							else:
								pass
						if "nowsaf" in self.counter_names:
							if not next_site.has_wsaf:
								self.counters["nnowsaf"].increment()
							else:
								pass
						if "noancestral" in self.counter_names:
							if not next_site.has_ancestral and not self.assume_ref_ancestral:
								self.counters["nnoanc"].increment()
							else:
								pass
						if "nonpolymorphic" in self.counter_names:
							if not next_site.is_polymorphic:
								self.counters["nnonpoly"].increment()
							else:
								pass
					yield next_site
				except StopIteration:
					#self.show_count_summary()
					self.ii -= 1
					raise StopIteration

	def write_header(self):
		sys.stdout.write( str(self.raw_header) )

	def write_site(self, site):
		sys.stdout.write( str(site) )

	def add_chunk(self, template, size, dtype = np.float, nsamples = None):
		if nsamples is None:
			nsamples = self.nsamples
		if template is None:
			return np.zeros( (size, nsamples), dtype = dtype )
		else:
			return np.append( template, self.add_chunk(None, size, template.dtype, template.shape[1]), axis = 0 )

	def show_count_summary(self):

		self.log.info("Done; site summary:".format(self.ii))
		for count_what in self.counter_names:
			counter = self.counters[count_what]
			abs_count = counter.count
			rel_count = 100.0 * abs_count/self.ii
			self.log.info("\t{: >12d} {} ({:0.3f}%)".format(abs_count, counter.description, rel_count))
		self.log.info("\t------------")
		self.log.info("\t{: >12d} total sites".format(self.ii))

	def reconcile_samples(self, samples, samples_file, allow_all = False):

		if samples is None:
			samples = []

		if hasattr(samples_file, "read"):
			keep_samples = []
			for line in samples_file:
				if line.startswith("#"):
					continue
				else:
					iid = line.strip().split().pop(0)
					if not iid in self.samples:
						self.log.warning("Sample '{}' in populations file but not in VCF; skpping it.".format(iid))
						continue
					keep_samples.append(iid)
			self.log.info("Read {} samples from file <{}>".format(len(keep_samples), samples_file.name))
		elif len(samples):
			keep_samples = []
			for iid in samples:
				keep_samples.extend(resplit("[,;|]", iid))
			self.log.info("Read {} samples from command line".format(len(keep_samples)))
		else:
			if not allow_all:
				self.log.error("Must specify at least one focal sample.")
				raise ValueError
			else:
				self.log.info("No samples specified, so keeping all of them.")
				keep_samples = self.samples

		keep_samples = set(self.samples) & set(keep_samples)
		keep_samples = list(keep_samples)
		nsamples = len(keep_samples)
		self.log.info("Retained {} distinct samples that were verified in VCF header.".format(nsamples))

		self._working_samples = keep_samples
		return keep_samples

	def reconcile_populations(self, pops_file = None):

		if not pops_file:
			self.log.info("Treating all samples in VCF as coming from one population.")
			pops = defaultdict(list)
			pops["pop1"].extend( list(range(0, self.nsamples)) )
		else:
			self.log.info("Reading population assignments from <{}>".format(pops_file.name))
			## decide which samples to use
			pops = defaultdict(list)
			seen = defaultdict(int)
			for line in pops_file:
				if line.startswith("#"):
					continue
				pieces = line.strip().split()
				# not an iid,pop tuple? skip
				if len(pieces) < 2:
					continue
				iid, pop = pieces[:2]
				# sample not in VCF? skip
				if not iid in self.samples:
					self.log.warning("Sample '{}' in populations file but not in VCF; skpping it.".format(iid))
					continue
				# already seen this sample? assume first assignment was right
				if iid in seen:
					continue
				else:
					seen[iid] += 1
				# finally, add sample to population
				pops[pop].append( self.samples.index(iid) )

		pop_order = sorted(pops, key = lambda k: -1*len(pops[k]))
		pop_sizes = [ len(pops[_]) for _ in pop_order ]
		self.log.info("Populations ({}) and their sizes:".format(len(pop_order)))
		for ii,pop in enumerate(pop_order):
			self.log.info("\t  {: >4d} {}".format(len(pops[pop]), pop))
		self.log.info("\t------")
		self.log.info("\t{: >6d} total samples in {} populations".format(sum(pop_sizes), len(pop_order)))
		self._pop_order = pop_order
		self._pops = pops
		return pop_order, pops

class VcfSite:

	def __init__(self, site, assume_ref_ancestral = False):
		self._site = site
		self._cached = None

		## expose properties of cyvcf2.Variant; apparently can't inherit in the usual way
		## this is bad form, but not clear how to avoid it ...
		self.CHROM = self._site.CHROM
		self.POS = self._site.POS
		self.REF = self._site.REF
		self.ALT = self._site.ALT
		self.INFO = self._site.INFO
		self.FORMAT = self._site.FORMAT
		self.ploidy = self._site.ploidy
		self.nsamples = len(self._site.gt_types)
		self.gt_types = self._site.gt_types
		self.gt_bases = self._site.gt_bases
		self.is_snp = self._site.is_snp
		self.is_indel = self._site.is_indel

		## some useful site properties
		self.is_multiallelic = len(self._site.ALT) > 1
		self.is_polymorphic = len(self._site.ALT) >= 1
		self.has_wsaf = self.has_format("WSAF")

		## annotations of ancestral allele
		self.has_ancestral = False
		self.ref_is_ancestral = False
		self.alt_is_ancestral = False
		self.ancestral = None
		self.derived = None
		anc = self.get_info("AA")
		self.set_ancestral(anc, assume_ref_ancestral)

		## add numeric genotypes with no-calls masked out
		self.gt_masked = np.ma.masked_values(self._site.gt_types, 3)

	def set_ancestral(self, anc, assume_ref_ancestral):
		if assume_ref_ancestral:
			anc = self._site.REF
		if anc is not None and anc != "X":
			self.ancestral = anc
			self.has_ancestral = True
			if self._site.REF == anc:
				self.ref_is_ancestral = True
				self.derived = self._site.ALT[0]
			elif self._site.ALT[0] == anc:
				self.derived = self._site.REF
				self.alt_is_ancestral = True

	def get_wsaf(self, samples = None):
		if samples is None:
			samples = list(range(0, self.nsamples))
		if self.has_format("WSAF"):
			wsaf = self._site.format("WSAF")[samples,0]
		else:
			wsaf = np.tile(-1, len(samples))
		wsaf_masked = np.ma.masked_values(wsaf, -1)
		return wsaf_masked

	def get_imputed_wsaf(self, samples = None):
		wsaf = self.get_wsaf()
		plaf = np.average(wsaf)
		if plaf is np.ma.masked:
			plaf = -1
		wsaf[ wsaf.mask ] = plaf
		return wsaf, plaf

	def recalc_plmaf(self, samples = None):
		wsaf, plaf = self.get_imputed_wsaf(samples)
		if wsaf is np.ma.masked:
			return -1
		else:
			plmaf = min(plaf, 1-plaf)
			return plmaf

	def to_n_alt(self):
		return self.gt_masked/2 * self.ploidy

	def to_n_der(self):
		nalt = self.to_n_alt()
		if self.ref_is_ancestral:
			return nalt
		elif self.alt_is_ancestral:
			return self.ploidy - nalt
		else:
			return None

	def format(self, field):
		return self._site.format(field)

	def set_format(self, name, data):
		self._site.set_format(name, data)

	def get_info(self, field):
		## safely retrieve stuff from INFO field; give None instead of KeyError if not present
		try:
			return self._site.INFO[field]
		except KeyError as e:
			return None

	def has_format(self, field):
		return field in self._site.FORMAT

	def __str__(self):
		if self._cached:
			return self._cached
		else:
			return self._site.__str__()

	def cache(self):
		self._cached = self._site.__str__()

	def write_out(self):
		sys.stdout.write( self.__str__() )
