#!/usr/bin/env python3
#
#  types.py
"""
Classes for type annotations.
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
from typing import Any, Dict, List, MutableMapping, MutableSequence, Optional, Tuple, TypedDict

# 3rd party
from libgunshotmatch.consolidate import ConsolidatedPeak

__all__ = [
		"PaddedPeakList",
		"MassList",
		"IntensityList",
		"ReferenceData",
		"PeakData",
		"JSONData",
		"TopMassesData",
		"CSVReportRow",
		]

PaddedPeakList = MutableSequence[Optional[ConsolidatedPeak]]
MassList = List[int]
IntensityList = List[float]


class ReferenceData(TypedDict):
	"""
	Represents reference mass spectrometry data, such as from the NIST Library.
	"""

	#: The compound name.
	name: str

	#: The compound's Chemical Abstracts Service Registry Number.
	cas: str

	#: The compound's chemical formula.
	formula: str

	#: The contributor of the reference data.
	contributor: str

	#: The NIST number of the compound.
	nist_no: int

	#: Identification string (e.g. library primary key) for the compound.
	id: str

	#: The compound's molecular weight.
	mw: int

	#: The compound's exact mass.
	exact_mass: float

	#: List of synonyms of the compound's name.
	synonyms: List[str]

	#: The compound's mass spectrum. TODO: Structure?
	mass_spec: Dict[str, List[int]]


class PeakData(TypedDict):
	"""
	Represents data for a peak.
	"""

	peak_no: int
	name: str

	#: The retention time in minutes.
	rt: float

	#: The peak area.
	area: float

	#: The peak area as a percentage.
	area_percentage: float

	match_factor: float
	reference_data: Dict[str, Any]


class JSONData(TypedDict):
	row: int
	peak: MutableMapping[str, Optional[PeakData]]
	ms: Dict[str, Tuple[MassList, IntensityList]]


TopMassesData = List[List[Tuple[str, str, str, str]]]

CSVReportRow = Tuple[
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str,
	str
	]
