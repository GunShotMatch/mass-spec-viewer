#!/usr/bin/env python3
#
#  confidence.py
"""
Confidence metric for similarity of profiles.
"""
#
#  Copyright © 2025 Dominic Davis-Foster <dominic@davis-foster.co.uk>
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
from typing import Dict, Iterator, MutableMapping, Optional, Tuple

# 3rd party
from domdf_python_tools.utils import redirect_output
from gunshotmatch_reports.alignment import _max_peak_area
from libgunshotmatch.comparison import align_projects, get_padded_peak_lists
from libgunshotmatch.project import Project

# this package
from mass_spec_viewer.data import _get_combined_spectra, _peak_data_for_peak, get_within_similarity
from mass_spec_viewer.types import IntensityList, JSONData, MassList, PaddedPeakList, PeakData

__all__ = ["calculate_confidence", "get_spectra_data"]


def get_spectra_data(
		reference_profile: Project,
		padded_ref_cp: PaddedPeakList,
		unknown: Project,
		padded_unkn_cp: PaddedPeakList,
		) -> Iterator[JSONData]:
	"""
	Returns an iterator over spectra and related data for peaks in the reference profile and unknown data.

	:param reference_profile:
	:param padded_ref_cp:
	:param unknown:
	:param padded_unkn_cp:
	"""

	reference_profile_max_pa = _max_peak_area(reference_profile)
	ms_data: Dict[str, Tuple[MassList, IntensityList]]
	peak_info: MutableMapping[str, Optional[PeakData]]

	assert padded_unkn_cp is not None

	unkn_max_pa = _max_peak_area(unknown)

	for idx, (reference_profile_peak, unkn_peak) in enumerate(zip(padded_ref_cp, padded_unkn_cp)):

		reference_profile_combined_spectrum, unkn_combined_spectrum, = _get_combined_spectra((reference_profile_peak, unkn_peak))

		ms_data = {}
		peak_info = {}

		for project, peak, spec, max_area in [
			(reference_profile, reference_profile_peak, reference_profile_combined_spectrum, reference_profile_max_pa),
			(unknown, unkn_peak, unkn_combined_spectrum, unkn_max_pa),
			]:
			assert project.consolidated_peaks is not None

			ms_data[project.name] = spec
			peak_info[project.name] = _peak_data_for_peak(project, peak, max_area)

		yield {"row": idx + 3, "peak": peak_info, "ms": ms_data}


def calculate_confidence(reference_profile: Project, unknown: Project) -> float:
	"""
	Calculate the confidence metric for similarity of profiles.

	:param referece_profile:
	:param unknown:
	"""

	with redirect_output() as (stdout, stderr):
		A1 = align_projects((reference_profile, ), (unknown, ), D=5)

	(padded_p1_cp, ), (padded_unkn_cp, ) = get_padded_peak_lists(
			A1,
			(reference_profile, ),
			(unknown, ),
			)

	max_ref_peak_area = _max_peak_area(reference_profile)
	max_unkn_peak_area = _max_peak_area(unknown)

	theoretical_confidence: float = 0
	actual_confidence: float = 0

	MAX_CONFIDENCE: float = 10

	for ref_peak, unkn_peak, spectra in zip(padded_p1_cp, padded_unkn_cp, get_spectra_data(reference_profile, padded_p1_cp, unknown, padded_unkn_cp)):
		if ref_peak:
			theoretical_confidence += MAX_CONFIDENCE
		else:
			continue

		if unkn_peak:
			confidence = MAX_CONFIDENCE
			# Confidence comprises peak area and mass spec similarity
			# Lose a point for every 100 MF away from 1000
			# Lose a point for every 10 percentage points difference
			# Overall confidence is percentage of what could be obtained (10 points for every peak in reference spectrum)
			# TODO: What to do about peaks in unknown only? Forward and reverse scores?
			ref_area_percentage = ref_peak.area / max_ref_peak_area
			unkn_area_percentage = unkn_peak.area / max_unkn_peak_area
			within_similarity = get_within_similarity(spectra)

			percentage_difference = abs(ref_area_percentage - unkn_area_percentage)
			confidence_penalty_area = min(percentage_difference / 0.1, 5)
			assert confidence_penalty_area <= 5
			confidence -= confidence_penalty_area

			match_factor = within_similarity.forward_scores[reference_profile.name][unknown.name]
			confidence_penalty_mf = min((1000 - match_factor) / 100, 5)
			assert confidence_penalty_mf <= 5
			confidence -= confidence_penalty_mf

			assert confidence >= 0
		else:
			confidence = 0

		actual_confidence += confidence

	return actual_confidence / theoretical_confidence
