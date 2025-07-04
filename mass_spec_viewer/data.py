#!/usr/bin/env python3
#
#  data.py
"""
Data gathering.
"""
#
#  Copyright © 2024 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#  OR OTHER DEALINGS IN THE SOFTWARE.
#

# stdlib
import itertools
from collections import defaultdict
from typing import Dict, Iterator, List, MutableMapping, NamedTuple, Optional, Sequence, Tuple, Type, Union, cast

# 3rd party
import numpy
from chemistry_tools.spectrum_similarity import SpectrumSimilarity
from gunshotmatch_reports.alignment import CSVRow, _max_peak_area, get_csv_data
from libgunshotmatch.consolidate import ConsolidatedPeak, combine_spectra
from libgunshotmatch.project import Project

# this package
from mass_spec_viewer.types import (
		CSVReportRow,
		IntensityList,
		JSONData,
		MassList,
		PaddedPeakList,
		PeakData,
		TopMassesData
		)

if "TYPE_CHECKING":
	# 3rd party
	from pyms_nist_search import ReferenceData  # nodep

__all__ = [
		"PeakInfo",
		"SimilarityScores",
		"csv_reports",
		"get_max_mass",
		"get_similarity",
		"get_spectra_data",
		"get_top_masses",
		"get_top_masses_data",
		"get_within_similarity",
		"normalize_intensities"
		]


class PeakInfo(NamedTuple):
	"""
	Information about a GC-MS peak.
	"""

	peak_no: int
	name: str

	#: Retention time in minutes
	rt: float
	area: float
	area_percentage: float  # as decimal
	match_factor: float
	reference_data: "ReferenceData"

	def to_dict(self) -> PeakData:
		"""
		Return a dictionary representation of the :class:`~.PeakInfo` object.
		"""

		return cast(PeakData, self._asdict())

	@classmethod
	def for_peak(cls: Type["PeakInfo"], peak: ConsolidatedPeak, peak_no: int, max_area: float) -> "PeakInfo":
		"""
		Create :class:`~.PeakInfo` for the given peak.

		:param project:
		:param peak:
		:param max_area: The area of the largest peak.
		"""

		first_hit = peak.hits[0]
		area = peak.area
		area_percentage = area / max_area
		assert first_hit.reference_data is not None
		return PeakInfo(
				peak_no,
				first_hit.name,
				peak.rt / 60,
				area,
				area_percentage,
				first_hit.match_factor,
				first_hit.reference_data,
				)


def get_max_mass(mass_list: MassList, intensity_list: IntensityList, cutoff: float = 0.01) -> int:
	"""
	Returns the maximum mass in the mass spectrum, excluding masses below the given intensity threshold.

	:param mass_list:
	:param intensity_list:
	:param cutoff: Masses less intense than this fraction of the maxiumum intensity will be ignored.
	"""

	max_intensity = max(intensity_list)
	max_mass = max(m for m, i in zip(mass_list, intensity_list) if i >= cutoff * max_intensity)
	return max_mass


_ExpRefSpectra = Dict[str, Optional[numpy.ndarray]]


class SimilarityScores(NamedTuple):
	"""
	Matrix of mass spectrum similarity scores between experimental and reference spectra.
	"""

	forward_scores: Dict[str, Dict[str, float]]
	reverse_scores: Dict[str, Dict[str, float]]

	@staticmethod
	def _get_experimental_spectra(json_data: JSONData) -> _ExpRefSpectra:
		experimental_spectra = {k: numpy.array(v).transpose() for k, v in json_data["ms"].items()}
		return {k: v if v.size else None for k, v in experimental_spectra.items()}

	@staticmethod
	def _get_reference_spectra(json_data: JSONData) -> _ExpRefSpectra:

		reference_spectra: _ExpRefSpectra = {}

		for sample, sample_peak_info in json_data["peak"].items():
			if sample_peak_info is None:
				reference_spectra[sample] = None
			else:
				ref_spec_array = numpy.array([
						sample_peak_info["reference_data"]["mass_spec"]["mass_list"],
						sample_peak_info["reference_data"]["mass_spec"]["intensity_list"]
						])
				reference_spectra[sample] = ref_spec_array.transpose()

		return reference_spectra

	@staticmethod
	def _get_scores(
			experimental_spectra: _ExpRefSpectra,
			reference_spectra: _ExpRefSpectra,
			) -> Iterator[Tuple[str, str, float, float]]:

		# fillvalue = None
		fillvalue = -1
		for sample_name, sample_ms in experimental_spectra.items():
			if sample_ms is None:
				for reference_sample_name, reference_ms in reference_spectra.items():
					yield sample_name, reference_sample_name, fillvalue, fillvalue
				continue

			for reference_sample_name, reference_ms in reference_spectra.items():
				if reference_ms is None:
					yield sample_name, reference_sample_name, fillvalue, fillvalue
					continue

				# print(sample_name, "/", reference_sample_name)
				ss = SpectrumSimilarity(sample_ms, reference_ms)
				mf, rmf = ss.score()
				yield sample_name, reference_sample_name, mf * 1000, rmf * 1000


def get_similarity(json_data: JSONData) -> SimilarityScores:
	"""
	Calculate matrix of mass spectrum similarity scores between experimental and reference spectra.

	:param json_data:
	"""

	experimental_spectra = SimilarityScores._get_experimental_spectra(json_data)
	reference_spectra = SimilarityScores._get_reference_spectra(json_data)

	forward_scores: Dict[str, Dict[str, float]] = defaultdict(dict)
	reverse_scores: Dict[str, Dict[str, float]] = defaultdict(dict)

	for sample_name, reference_sample_name, mf, rmf in SimilarityScores._get_scores(experimental_spectra, reference_spectra):

		# def fill(sample_name: str, reference_sample_name: str) -> None:
		# 	forward_scores[sample_name][reference_sample_name] = fillvalue
		# 	reverse_scores[sample_name][reference_sample_name] = fillvalue

		# # fillvalue = None
		# fillvalue = -1
		# for sample_name, sample_ms in experimental_spectra.items():
		# 	if sample_ms is None:
		# 		for reference_sample_name, reference_ms in reference_spectra.items():
		# 			fill(sample_name, reference_sample_name)
		# 		continue

		# 	for reference_sample_name, reference_ms in reference_spectra.items():
		# 		if reference_ms is None:
		# 			fill(sample_name, reference_sample_name)
		# 			continue

		# 		# print(sample_name, "/", reference_sample_name)
		# 		ss = SpectrumSimilarity(sample_ms, reference_ms)
		# 		mf, rmf = ss.score()
		# 		mf *= 1000
		# 		rmf *= 1000

		forward_scores[sample_name][reference_sample_name] = mf
		reverse_scores[sample_name][reference_sample_name] = rmf

	# print({k: dict(v) for k, v in forward_scores.items()})

	return SimilarityScores(forward_scores, reverse_scores)


def get_within_similarity(json_data: JSONData) -> SimilarityScores:
	"""
	Calculate matrix of mass spectrum similarity scores between experimental spectra.

	:param json_data:
	"""

	experimental_spectra = SimilarityScores._get_experimental_spectra(json_data)

	forward_scores: Dict[str, Dict[str, float]] = defaultdict(dict)
	reverse_scores: Dict[str, Dict[str, float]] = defaultdict(dict)

	for sample_name, reference_sample_name, mf, rmf in SimilarityScores._get_scores(experimental_spectra, experimental_spectra):

		# def fill(sample_name: str, reference_sample_name: str) -> None:
		# 	forward_scores[sample_name][reference_sample_name] = fillvalue
		# 	reverse_scores[sample_name][reference_sample_name] = fillvalue

		# # fillvalue = None
		# fillvalue = -1
		# for sample_name, sample_ms in experimental_spectra.items():
		# 	if sample_ms is None:
		# 		for reference_sample_name, reference_ms in experimental_spectra.items():
		# 			fill(sample_name, reference_sample_name)
		# 		continue

		# 	for reference_sample_name, reference_ms in experimental_spectra.items():
		# 		if reference_ms is None:
		# 			fill(sample_name, reference_sample_name)
		# 			continue

		# 		# print(sample_name, "/", reference_sample_name)
		# 		ss = SpectrumSimilarity(sample_ms, reference_ms)
		# 		mf, rmf = ss.score()
		# 		mf *= 1000
		# 		rmf *= 1000

		if sample_name == reference_sample_name:
			if mf != -1 and rmf != -1:
				assert numpy.isclose(mf, 1000, atol=0.01), mf
				assert numpy.isclose(rmf, 1000, atol=0.01), rmf
				forward_scores[sample_name][reference_sample_name] = -2
				reverse_scores[sample_name][reference_sample_name] = -2
				continue

		forward_scores[sample_name][reference_sample_name] = mf
		reverse_scores[sample_name][reference_sample_name] = rmf

	# print({k: dict(v) for k, v in forward_scores.items()})

	return SimilarityScores(forward_scores, reverse_scores)


def normalize_intensities(intensity_list: List[float]) -> List[float]:
	"""
	Normalise intensities in the intensity list (0-100).

	:param intensity_list:
	"""

	max_intensity = max(intensity_list) if intensity_list else 0
	return [(intensity / max_intensity) * 100 for intensity in intensity_list]


def get_top_masses(mass_list: List[int], intensity_list: List[float]) -> Iterator[Tuple[int, int]]:
	"""
	Calculate the top 10 most intense masses in the mass spectrum.

	:param mass_list:
	:param intensity_list:
	"""

	intensity_list = normalize_intensities(intensity_list)
	ms_dict: Dict[int, float] = dict(zip(mass_list, intensity_list))
	result = sorted(
			ms_dict,
			key=ms_dict.get,  # type: ignore[arg-type]
			reverse=True
			)

	for mass in result[:10]:
		yield mass, int(ms_dict[mass] * 9.99)


def get_top_masses_data(json_data: JSONData) -> TopMassesData:
	"""
	Calculate the most intense mass fragments for each project in the given data.

	:param json_data:
	"""

	top_masses_data = []

	for project_name, spectrum in json_data["ms"].items():
		peak_info = cast(PeakData, json_data["peak"][project_name] or {})

		project_top_masses_data: List[Tuple[str, str, str, str]] = [
				(project_name, '', '', ''),
				(str(peak_info.get("name")), '', '', ''),
				("Experimental", '', "Reference", ''),
				("Mass", "Intensity", "Mass", "Intensity"),
				]

		experimental_top_masses = get_top_masses(*spectrum)

		if peak_info:
			reference_ms = peak_info["reference_data"]["mass_spec"]
		else:
			reference_ms = {"mass_list": [], "intensity_list": []}

		reference_mass_list, reference_intensity_list = [], []
		for ref_mass, ref_i in zip(reference_ms["mass_list"], reference_ms["intensity_list"]):
			if ref_mass >= 50:
				reference_mass_list.append(ref_mass)
				reference_intensity_list.append(ref_i)

		reference_top_masses = get_top_masses(reference_mass_list, reference_intensity_list)

		exp: Union[Tuple[int, int], Tuple[str, str]]
		ref: Union[Tuple[int, int], Tuple[str, str]]
		for exp, ref in itertools.zip_longest(experimental_top_masses, reference_top_masses, fillvalue=('', '')):
			project_top_masses_data.append((str(exp[0]), str(exp[1]), str(ref[0]), str(ref[1])))

		top_masses_data.append(project_top_masses_data)

	return top_masses_data


_Spectrum = Tuple[List[int], List[float]]


def get_spectra_data(
		p1: Project,
		padded_p1_cp: PaddedPeakList,
		p2: Project,
		padded_p2_cp: PaddedPeakList,
		u: Optional[Project] = None,
		padded_unkn_cp: Optional[PaddedPeakList] = None,
		) -> Iterator[JSONData]:
	"""
	Get spectra data for each peak in the aligned projects and unknown.

	:param p1: The first project
	:param padded_p1_cp: List of consolidated peaks for the first project,
		padded to align with the other project's and the unknown's.
	:param p2: The second project
	:param padded_p2_cp: List of consolidated peaks for the second project,
		padded to align with the other project's and the unknown's.
	:param u: The unknown sample
	:param padded_unkn_cp: List of consolidated peaks for the unknown,
		padded to align with the projects'.
	"""

	p1_max_pa = _max_peak_area(p1)
	p2_max_pa = _max_peak_area(p2)
	ms_data: Dict[str, Tuple[MassList, IntensityList]]
	peak_info: MutableMapping[str, Optional[PeakData]]

	if u:
		assert padded_unkn_cp is not None

		unkn_max_pa = _max_peak_area(u)

		for idx, (p1_peak, p2_peak, unkn_peak) in enumerate(zip(padded_p1_cp, padded_p2_cp, padded_unkn_cp)):
			if sum([peak is not None for peak in [p1_peak, p2_peak, unkn_peak]]) < 2:
				continue

			p1_combined_spectrum, unkn_combined_spectrum, p2_combined_spectrum = _get_combined_spectra((p1_peak, unkn_peak, p2_peak))

			# for ax in axes:
			# 	ax.set_xlim(45, max_mass+5)

			ms_data = {}
			peak_info = {}

			for project, peak, spec, max_area in [
				(p1, p1_peak, p1_combined_spectrum, p1_max_pa),
				(u, unkn_peak, unkn_combined_spectrum, unkn_max_pa),
				(p2, p2_peak, p2_combined_spectrum, p2_max_pa),
				]:
				assert project.consolidated_peaks is not None

				ms_data[project.name] = spec
				peak_info[project.name] = _peak_data_for_peak(project, peak, max_area)

			yield {"row": idx + 3, "peak": peak_info, "ms": ms_data}

	else:

		for idx, (p1_peak, p2_peak) in enumerate(zip(padded_p1_cp, padded_p2_cp)):
			# fig, axes = plt.subplots(3, 1, figsize=(8, 12), layout="constrained")

			p1_combined_spectrum, p2_combined_spectrum = _get_combined_spectra((p1_peak, p2_peak))

			# for ax in axes:
			# 	ax.set_xlim(45, max_mass+5)

			ms_data = {}
			peak_info = {}

			for project, peak, spec, max_area in [
				(p1, p1_peak, p1_combined_spectrum, p1_max_pa),
				(p2, p2_peak, p2_combined_spectrum, p2_max_pa),
				]:
				assert project.consolidated_peaks is not None

				ms_data[project.name] = spec
				peak_info[project.name] = _peak_data_for_peak(project, peak, max_area)

			yield {"row": idx + 3, "peak": peak_info, "ms": ms_data}


def _peak_data_for_peak(
		project: Project,
		peak: Optional[ConsolidatedPeak],
		max_area: float,
		) -> Optional[PeakData]:
	assert project.consolidated_peaks is not None

	if peak is None:
		return None
	else:

		project_peak_info = PeakInfo.for_peak(
				peak,
				project.consolidated_peaks.index(peak) + 1,
				max_area,
				)
		project_peak_data = project_peak_info.to_dict()
		return project_peak_data


def _get_csv_header(p1: Project, p2: Project, u: Optional[Project]) -> Tuple[CSVReportRow, CSVReportRow]:
	top_row_pad = ('', ) * (5)

	if u:
		return (
				['', p1.name, *top_row_pad, u.name, *top_row_pad, p2.name, *top_row_pad],
				['', *(CSVRow.header() * 3)],
				)
	else:
		return (
				['', p1.name, *top_row_pad, p2.name, *top_row_pad],
				['', *(CSVRow.header() * 2)],
				)


def csv_reports(
		p1: Project,
		padded_p1_cp: PaddedPeakList,
		p2: Project,
		padded_p2_cp: PaddedPeakList,
		u: Optional[Project] = None,
		padded_unkn_cp: Optional[PaddedPeakList] = None,
		) -> Tuple[Tuple[CSVReportRow, CSVReportRow], List[Tuple[CSVReportRow, bool]]]:
	"""
	Generate CSV reports showing the alignment between the two reference profiles and the unknown.

	:param p1: The first project.
	:param padded_p1_cp: The peak list for the first project, padded to ensure alignment with the other peak lists.
	:param p2: The second project.
	:param padded_p2_cp: The peak list for the second project, padded to ensure alignment with the other peak lists.
	:param u: The unknown sample.
	:param padded_unkn_cp: The peak list for the unknown sample, padded to ensure alignment with the other peak lists.
	"""

	p1_max_pa = _max_peak_area(p1)
	p2_max_pa = _max_peak_area(p2)

	if u:
		unkn_max_pa = _max_peak_area(u)
	else:
		padded_unkn_cp = [None] * len(padded_p1_cp)

	csv_header = _get_csv_header(p1, p2, u)
	csv_data = []

	assert padded_unkn_cp is not None

	for idx, (cp1, cp2, cpu) in enumerate(zip(padded_p1_cp, padded_p2_cp, padded_unkn_cp)):
		idx += 3

		cp1_data = get_csv_data(p1, cp1, p1_max_pa)
		cp2_data = get_csv_data(p2, cp2, p2_max_pa)

		row: CSVReportRow
		if u:
			unknown_data = get_csv_data(u, cpu, unkn_max_pa)
			row = [str(idx), *cp1_data, *unknown_data, *cp2_data]
			csv_data.append((row, sum(tuple(map(all, [cp1_data, cp2_data, unknown_data]))) >= 2))
		else:
			row = [str(idx), *cp1_data, *cp2_data]
			csv_data.append((row, True))

	return csv_header, csv_data


def _get_combined_spectra(peaks: Sequence[Optional[ConsolidatedPeak]]) -> List[_Spectrum]:
	combined_spectra = []
	max_mass = 0
	# for ax, peak in zip(axes, (p1_peak, unkn_peak, p2_peak)):

	peak: Optional[ConsolidatedPeak]
	for peak in peaks:

		combined_spectrum: _Spectrum
		if peak is None:
			combined_spectrum = [], []
			# ax.set_xlabel("m/z")
			# ax.set_ylabel("Intensity")
		else:
			combined_spectrum = combine_spectra(peak)
			# draw_ms(fig, ax, *combined_spectrum)
			max_mass = max(max_mass, get_max_mass(*combined_spectrum))

		combined_spectra.append(combined_spectrum)

	return combined_spectra
