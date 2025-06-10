#!/usr/bin/env python3
#
#  _batch_process.py
"""
Internal CLI utils.
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
import concurrent.futures
import csv
import io
from typing import List, Optional, TypedDict

# 3rd party
import mpld3.urls  # type: ignore[import-untyped]
import sdjson
from domdf_python_tools.compat import importlib_resources
from domdf_python_tools.paths import PathPlus
from libgunshotmatch.comparison import align_projects, get_padded_peak_lists
from libgunshotmatch.project import Project
from matplotlib import pyplot as plt
from natsort import natsorted
from tqdm import tqdm

# this package
from mass_spec_viewer import data
from mass_spec_viewer.index import BOOTSTRAP_CSS_URL, create_index_page
from mass_spec_viewer.plot import draw_mass_spectrum, plot_scores
from mass_spec_viewer.types import PaddedPeakList

__all__ = [
		"Comparison",
		"copy_assets",
		"draw_mass_spectra",
		"get_mass_spectra",
		"get_spectra_data",
		"make_csv_reports",
		"make_index_html",
		"process_comparison"
		]


def get_spectra_data(
		p1: Project,
		padded_p1_cp: PaddedPeakList,
		p2: Project,
		padded_p2_cp: PaddedPeakList,
		u: Optional[Project] = None,
		padded_unkn_cp: Optional[PaddedPeakList] = None,
		*,
		output_dir: PathPlus,
		) -> None:
	"""
	Write out spectra data for each peak in the aligned projects and unknown.

	:param p1: The first project
	:param padded_p1_cp: List of consolidated peaks for the first project,
		padded to align with the other project's and the unknown's.
	:param p2: The second project
	:param padded_p2_cp: List of consolidated peaks for the second project,
		padded to align with the other project's and the unknown's.
	:param u: The unknown sample
	:param padded_unkn_cp: List of consolidated peaks for the unknown,
		padded to align with the projects'.
	:param output_dir: The directory to write the JSON files into.
	"""

	for json_data in data.get_spectra_data(
			p1=p1,
			padded_p1_cp=padded_p1_cp,
			p2=p2,
			padded_p2_cp=padded_p2_cp,
			u=u,
			padded_unkn_cp=padded_unkn_cp,
			):
		idx = json_data["row"]
		json_filename = output_dir / f"{idx}_spectra.json"
		json_filename.write_clean(sdjson.dumps(json_data, indent=2))


def make_csv_reports(
		p1: Project,
		padded_p1_cp: PaddedPeakList,
		p2: Project,
		padded_p2_cp: PaddedPeakList,
		u: Optional[Project] = None,
		padded_unkn_cp: Optional[PaddedPeakList] = None,
		*,
		output_dir: PathPlus,
		) -> None:
	"""
	Generate CSV reports showing the alignment between the two reference profiles and the unknown.

	:param p1: The first project.
	:param padded_p1_cp: The peak list for the first project, padded to ensure alignment with the other peak lists.
	:param p2: The second project.
	:param padded_p2_cp: The peak list for the second project, padded to ensure alignment with the other peak lists.
	:param u: The unknown sample.
	:param padded_unkn_cp: The peak list for the unknown sample, padded to ensure alignment with the other peak lists.
	:param output_dir: The directory to write the CSV files into.
	"""

	csv_header, csv_data = data.csv_reports(p1, padded_p1_cp, p2, padded_p2_cp, u, padded_unkn_cp)

	full_csv_fp = io.StringIO()
	pair_only_csv_fp = io.StringIO()

	csvwriter_kwargs = dict(quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
	full_csvwriter = csv.writer(
			full_csv_fp,
			**csvwriter_kwargs,  # type: ignore[arg-type]
			)
	pair_only_csvwriter = csv.writer(
			pair_only_csv_fp,
			**csvwriter_kwargs,  # type: ignore[arg-type]
			)

	full_csvwriter.writerows(csv_header)
	pair_only_csvwriter.writerows(csv_header)

	for row, include_in_pair_only in csv_data:
		full_csvwriter.writerow(row)
		if include_in_pair_only:
			pair_only_csvwriter.writerow(row)

	(output_dir / "alignment.csv").write_clean(full_csv_fp.getvalue())
	(output_dir / "alignment_pair_only.csv").write_clean(pair_only_csv_fp.getvalue())


def get_mass_spectra(
		project1: Project,
		project2: Project,
		unknown: Optional[Project],
		output_dir: PathPlus,
		) -> None:
	"""
	Write CSV reports of the alignment between the projects and unknown, and JSON files giving the spectra for the individual peaks.

	:param project1: The first project
	:param project2: The second project
	:param unknown: The unknown sample
	:param output_dir: The directory to write the CSV files into.
	"""

	if unknown:
		A1 = align_projects([project1, project2], (unknown, ), D=5)
		(padded_p1_cp, padded_p2_cp), (padded_unkn_cp, ) = get_padded_peak_lists(
				A1,
				[project1, project2],
				(unknown, ),
				)
	else:
		A1 = align_projects([project1, project2], D=5)
		(padded_p1_cp, padded_p2_cp), () = get_padded_peak_lists(A1, [project1, project2])
		padded_unkn_cp = None

	make_csv_reports(
			project1,
			padded_p1_cp,
			project2,
			padded_p2_cp,
			unknown,
			padded_unkn_cp,
			output_dir=output_dir,
			)
	get_spectra_data(
			project1,
			padded_p1_cp,
			project2,
			padded_p2_cp,
			unknown,
			padded_unkn_cp,
			output_dir=output_dir,
			)


def _json_filename_to_html(filename: PathPlus) -> str:
	"""
	Return the HTML filename corresponding to the given spectrum JSON filename.
	"""

	return f"./{filename.name.strip('*_spectra.json')}.html"


def _draw_ms(
		idx: int,
		filename: PathPlus,
		filenames: List[PathPlus],
		output_dir: PathPlus,
		d3_url: str,
		mpld3_url: str,
		) -> None:
	json_data = filename.load_json()
	spectrum_basename = output_dir / str(json_data["row"])
	previous_file = _json_filename_to_html(filenames[idx - 1]) if idx > 0 else '#'
	next_file = _json_filename_to_html(filenames[idx + 1]) if idx < len(filenames) - 1 else '#'

	row_num = json_data["row"]
	html_output = draw_mass_spectrum(
			json_data,
			row_num,
			previous_file,
			next_file,
			d3_url=d3_url,
			mpld3_url=mpld3_url,
			)
	spectrum_basename.with_suffix(".html").write_clean(html_output)

	fig = plot_scores(json_data, row_num)
	plt.savefig(output_dir / f"{row_num}.matrix.png")
	plt.close(fig)


def draw_mass_spectra(
		output_dir: PathPlus,
		jobs: int = 1,
		d3_url: str = mpld3.urls.D3_URL,
		mpld3_url: str = mpld3.urls.MPLD3MIN_URL,
		) -> None:
	"""
	Plot mass spectra to PNG files.

	:param output_dir: The directory in which the :file:`*_spectra.json` files are to be found, and in which the PNG files are to be written.
	:param jobs: The number of jobs to use, to allow for parallel execution.
	:param d3_url: URL to use for the ``d3.js`` library.
	:param mpld3_url: URL to use for the mpld3 library.
	"""

	filenames = natsorted(output_dir.glob("*_spectra.json"))

	if jobs == 1:
		for idx, filename in enumerate(tqdm(filenames)):
			_draw_ms(idx, filename, filenames, output_dir, d3_url, mpld3_url)
	else:
		with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
			futures = {
					executor.submit(_draw_ms, idx, filename, filenames, output_dir, d3_url, mpld3_url)
					for idx,
					filename in enumerate(filenames)
					}
			for future in tqdm(concurrent.futures.as_completed(futures), total=len(filenames)):
				future.result()


def make_index_html(
		output_dir: PathPlus,
		bootstrap_css_url: str = BOOTSTRAP_CSS_URL,
		) -> None:
	r"""
	Write the ``index.html`` page with an overview table of compounds and retention times.

	:param output_dir: Directory conaining \*.spectra.json files and in which to place the ``index.html`` file.
	:param table_data: For each tuple the first element is the row number,
		the second is list of tuples of (retention time, compound name).
	:param bootstrap_css_url: Optional URL for boostrap css. May include additional HTML attributes for the ``link`` tag.
	"""

	index_table_data = []
	project_names = []

	for filename in natsorted(output_dir.glob("*_spectra.json")):
		json_data = filename.load_json()
		table_row = []
		project_names = list(json_data["peak"])
		for project_data in json_data["peak"].values():
			if project_data:
				table_row.append((f'{project_data["rt"]:0.3f}', project_data["name"]))
			else:
				table_row.append(('', ''))

		index_table_data.append((json_data["row"], table_row))

	index_html_output = create_index_page(project_names, index_table_data, bootstrap_css_url=bootstrap_css_url)
	(output_dir / "index.html").write_clean(index_html_output)


class Comparison(TypedDict):
	"""
	Represents a comparison between the projects and unknown.
	"""

	#: The True class.
	project1: str

	#: For incorrect predictions, the predicted class. For correct predictions, the 2nd prediction.
	project2: str

	#: The unknown sample for comparison
	unknown: Optional[str]

	#: Output directory for this comparison.
	output_dir: str

	# Optional name for project1. Useful if two projects with the same name are being compared.
	project1_name_override: Optional[str]

	# Optional name for project2. Useful if two projects with the same name are being compared.
	project2_name_override: Optional[str]

	# Optional name for unknown. Useful if two projects with the same name are being compared.
	unknown_name_override: Optional[str]


# def download_assets(output_dir: PathPlus):
# 	import httpx, hashlib

# 	bootstrap_url = BOOTSTRAP_CSS_URL.split('" rel="')[0]
# 	bootstrap_hash_algo, bootstrap_hash = BOOTSTRAP_CSS_URL.split('integrity="')[1].split('" crossorigin="')[0].split("-", 1)

# 	response = httpx.get(bootstrap_url)
# 	response.raise_for_status()

# 	hash_algo_function = hashlib.new(bootstrap_hash_algo)
# 	hash_algo_function.update(response.text.encode())
# 	actual_hash = base64.b64encode(hash_algo_function.digest()).decode()
# 	if actual_hash != bootstrap_hash:
# 		raise ValueError(f"Hash mismatch for bootstrap.min.css:\n\tExpected {bootstrap_hash!r}\n\tGot {actual_hash!r}")

# 	(output_dir / "bootstrap.min.css").write_text(response.text)
_BS_MIN_CSS = "bootstrap.min.css"


def copy_assets(output_dir: PathPlus) -> List[str]:
	"""
	Copy CSS and JS assets for future offline use.

	:param output_dir:

	:returns: List of Boostrap CSS, d3 js, and mpld3 js filenames (in that order).
	"""

	bootstrap_css = importlib_resources.read_text("mass_spec_viewer.css", _BS_MIN_CSS)
	(output_dir / _BS_MIN_CSS).write_clean(bootstrap_css)
	return_names = [_BS_MIN_CSS]

	for filename in [mpld3.urls.D3_LOCAL, mpld3.urls.MPLD3MIN_LOCAL]:
		path = PathPlus(filename)
		(output_dir / path.name).write_clean(path.read_text())
		return_names.append(path.name)

	return return_names  # boostrap css, d3 js, mpld3 js


def process_comparison(
		comparison: Comparison,
		jobs: int = 1,
		should_copy_assets: bool = False,
		) -> None:
	"""
	Produce HTML, PNG and JSON files for the given comparison.

	:param comparison:
	:param jobs: The number of jobs to use, to allow for parallel execution.
	:param should_copy_assets: Whether to copy the Bootstrap CSS, d3.js, and mpl3 Javascript files into the output directory,
		to allow for offline viewing.
	"""

	output_dir = PathPlus(comparison["output_dir"])
	output_dir.maybe_make(parents=True)

	if should_copy_assets:
		bootstrap_css, d3js, mpl3js = copy_assets(output_dir)
	else:
		bootstrap_css = BOOTSTRAP_CSS_URL
		d3js, mpl3js = mpld3.urls.D3_URL, mpld3.urls.MPLD3MIN_URL

	project1 = Project.from_file(comparison["project1"])
	project1_name_override = comparison.get("project1_name_override")
	if project1_name_override:
		project1.name = project1_name_override

	project2 = Project.from_file(comparison["project2"])
	project2_name_override = comparison.get("project2_name_override")
	if project2_name_override:
		project2.name = project2_name_override

	if project2.name == project1.name:
		project2.name = f"{project2.name} (1)"

	if "unknown" in comparison:
		assert comparison["unknown"] is not None
		unknown = Project.from_file(comparison["unknown"])
		unknown_name_override = comparison.get("unknown_name_override")
		if unknown_name_override:
			unknown.name = unknown_name_override

		if unknown.name == project1.name:
			unknown.name = f"{unknown.name} (1)"

	else:
		unknown = None

	get_mass_spectra(project1, project2, unknown, output_dir)
	del project1
	del project2
	del unknown

	make_index_html(output_dir, bootstrap_css_url=bootstrap_css)
	draw_mass_spectra(output_dir, jobs=jobs, d3_url=d3js, mpld3_url=mpl3js)
