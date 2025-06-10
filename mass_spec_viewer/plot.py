#!/usr/bin/env python3
#
#  plot.py
"""
Plotting functions.
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
import json
import re
from typing import Any, Dict, List, cast

# 3rd party
import matplotlib
import mpld3  # type: ignore[import-untyped]
import numpy
from domdf_python_tools.stringlist import StringList
from domdf_python_tools.words import truncate_string
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.container import BarContainer
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.ticker import MultipleLocator
from mpld3 import plugins
from mpld3.mpld3renderer import MPLD3Renderer  # type: ignore[import-untyped]
from mpld3.mplexporter import Exporter  # type: ignore[import-untyped]

# this package
from mass_spec_viewer.data import get_max_mass, get_similarity, get_top_masses_data, get_within_similarity
from mass_spec_viewer.types import JSONData, PeakData, TopMassesData

__all__ = [
		"add_html_tooltips",
		"draw_mass_spectrum",
		"draw_ms",
		"format_top_masses_html_table",
		"make_similarity_grid",
		"plot_scores",
		"plot_spectra"
		]

html_tooltip_css = """
table {
  border-collapse: collapse;
}

th {
  color: #ffffff;
  background-color: #000000;
}

td {
  background-color: #cccccc;
}

table, th, td {
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
"""


def add_html_tooltips(fig: Figure, bars: BarContainer, mass_list: List[int]) -> None:
	"""
	Add mass tooltips for the HTML output.

	:param fig:
	:param bars:
	:param mass_list:
	"""

	for box, mass in zip(bars.get_children(), mass_list):
		tooltip = mpld3.plugins.LineLabelTooltip(box, label=str(mass))
		mpld3.plugins.connect(fig, tooltip)

	tooltip = plugins.PointHTMLTooltip(
			bars[0],
			["hello world"] * len(mass_list),
			voffset=100,
			hoffset=100,
			css=html_tooltip_css,
			)
	plugins.connect(fig, tooltip)


def draw_ms(ax: Axes, mass_list: List[int], intensity_list: List[float]) -> BarContainer:
	"""
	Plot a mass spectrum.

	:param ax:
	:param mass_list:
	:param intensity_list:
	"""

	max_intensity = max(intensity_list)
	chart_intensity_list = [(intensity / max_intensity) * 100 for intensity in intensity_list]

	bars = ax.bar(
			mass_list,
			chart_intensity_list,
			width=0.5,
			)
	ax.set_xlabel("m/z")
	ax.set_ylabel("Intensity")

	ax.xaxis.set_major_locator(MultipleLocator(5))
	ax.xaxis.set_minor_locator(MultipleLocator(1))
	return bars


def plot_spectra(json_data: JSONData) -> Figure:
	"""
	Plot mass spectra from the given data.

	:param json_data:

	:returns: The figure the spectrum was plotted on. Figure size 16.5" x 9".
	"""

	# fig, axes = plt.subplots(3, 2, figsize=(16.5, 11.7), layout="constrained")
	fig, axes = plt.subplots(len(json_data["peak"]), 2, figsize=(16.5, 9), layout="constrained")
	experimental_axes: List[Axes] = axes[:, 0]  # Left column
	reference_axes: List[Axes] = axes[:, 1]  # Right column

	max_mass = 0

	for ax, (project_name, spectrum) in zip(experimental_axes, json_data["ms"].items()):
		peak_info = json_data["peak"][project_name]
		if peak_info:

			add_html_tooltips(fig, draw_ms(ax, *spectrum), spectrum[0])
			max_mass = max(max_mass, get_max_mass(*spectrum))
			ax.set_title(f"{project_name} – {truncate_string(peak_info['name'], 70)} @ {peak_info['rt']:0.3f} min")
		else:
			ax.set_xlabel("m/z")
			ax.set_ylabel("Intensity")
			ax.set_title(f"{project_name} – No Peak")

	for ax, (project_name, peak_info) in zip(reference_axes, json_data["peak"].items()):

		if peak_info:
			reference_ms = peak_info["reference_data"]["mass_spec"]

			add_html_tooltips(
					fig,
					draw_ms(ax, reference_ms["mass_list"], reference_ms["intensity_list"]),
					reference_ms["mass_list"],
					)
			max_mass = max(max_mass, get_max_mass(reference_ms["mass_list"], reference_ms["intensity_list"]))
			ax.set_title(f"Reference – {truncate_string(peak_info['reference_data']['name'], 70)}")
		else:
			ax.set_xlabel("m/z")
			ax.set_ylabel("Intensity")
			ax.set_title(f"No Peak")

	for ax in [*experimental_axes, *reference_axes]:
		ax.set_xlim(45, max_mass + 5)

	return fig


def make_similarity_grid(ax: Axes, scores: Dict[str, Dict[str, float]], labels: List[str]) -> AxesImage:
	"""
	Create a grid of mass spectral similarity scores.

	:param ax:
	:param scores: Mapping of y names to a mapping of x names to scores.
	:param labels: Labels for the samples (used on both x and y axes).
	"""

	grid_data = [list(v.values()) for v in scores.values()]

	im = ax.imshow(grid_data, cmap=matplotlib.colormaps["Blues"], vmin=0, vmax=1000)

	# Show all ticks and label them with the respective list entries
	ax.set_xticks(numpy.arange(len(labels)), labels=labels)
	ax.set_yticks(numpy.arange(len(labels)), labels=labels)

	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

	thresh = 500
	assert im.cmap is not None
	assert im.cmap is not None
	cmap_min, cmap_max = im.cmap(0), im.cmap(1.0)

	# Loop over data dimensions and create text annotations.
	for i in range(len(labels)):
		for j in range(len(labels)):
			color = cmap_max if grid_data[i][j] < thresh else cmap_min

			# Ref: https://github.com/python/mypy/issues/18481
			kwds: Dict[str, Any] = dict(ha="center", va="center", color=color)
			value = grid_data[i][j]

			if value == -1:
				ax.text(j, i, '', **kwds)
			elif value == -2:
				ax.text(j, i, "1000.0", **kwds)
			else:
				ax.text(j, i, f"{value:0.1f}", **kwds)

	return im


def plot_scores(json_data: JSONData, row_num: int) -> Figure:
	"""
	Plot mass spectral similarity scores calculated from the given data.

	Four grids are plotted. The top two grids show forward and reverse match factors for the experimental spectra against reference spectra for the top hits.
	The bottom two grids show forward and reverse match factors for the experimental spectra against experimental spectra.
	Comparing a spectrum to itself gives are score of 1000, so these are de-emphasised for ease of reading.

	:param json_data:

	:returns: The figure the spectrum was plotted on. Figure size 10" x 11".
	"""

	# fig, axes = plt.subplots(2, 2, figsize=(8, 9), layout="constrained", sharey=True)
	fig, axes = plt.subplots(2, 2, figsize=(10, 11), layout="constrained", sharey=True)

	forward_scores, reverse_scores = get_similarity(json_data)
	grid_labels = []

	for sample_name in forward_scores:
		peak_data = cast(PeakData, json_data["peak"][sample_name] or {})
		compound = truncate_string(peak_data.get("name", '').strip(), 30)

		label = f"{sample_name}\n{compound}".strip()
		grid_labels.append(label)

	make_similarity_grid(axes[0, 0], forward_scores, grid_labels)
	make_similarity_grid(axes[0, 1], reverse_scores, grid_labels)

	fig.suptitle(f"Mass Spectrum Similarity Scores – Row {row_num}")

	axes[0, 0].set_ylabel("Experimental Spectra")
	axes[0, 0].set_xlabel("Reference Spectra")
	axes[0, 0].set_title("Forward Match")

	axes[0, 1].set_xlabel("Reference Spectra")
	axes[0, 1].set_title("Reverse Match")

	forward_scores, reverse_scores = get_within_similarity(json_data)
	grid_labels = []

	for sample_name in forward_scores:
		peak_data = cast(PeakData, json_data["peak"][sample_name] or {})
		compound = truncate_string(peak_data.get("name", '').strip(), 30)
		label = f"{sample_name}\n{compound}".strip()
		grid_labels.append(label)

	make_similarity_grid(axes[1, 0], forward_scores, grid_labels)
	make_similarity_grid(axes[1, 1], reverse_scores, grid_labels)

	# fig.suptitle(f"Mass Spectrum Similarity Scores – Row {row_num}")

	axes[1, 0].set_ylabel("Experimental Spectra")
	axes[1, 0].set_xlabel("Experimental Spectra")
	axes[1, 0].set_title("Forward Match")

	axes[1, 1].set_xlabel("Reference Spectra")
	axes[1, 1].set_title("Reverse Match")

	# fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1, wspace=0, hspace=0)
	fig.get_layout_engine().set(w_pad=0.1, h_pad=0.1, wspace=0, hspace=0)  # type: ignore[call-arg,union-attr]

	return fig


def format_top_masses_html_table(top_masses_data: TopMassesData) -> str:
	"""
	Format data on the most intense fragment masses in the mass spectra as an HTML table.

	:param top_masses_data:
	"""

	output = []
	output.append("<table>\n  <tbody>")

	html_table_rows = list(itertools.zip_longest(*top_masses_data, fillvalue=[''] * 4))

	header_rows = html_table_rows[:3]

	# First header row, sample names
	output.append("<tr>")
	for sample_subrow in header_rows[0]:
		sample_name, _, _, _ = sample_subrow
		output.append(f'<td colspan="4">{sample_name}</td>')
	output.append("    </tr>")

	# Second header row, compound names
	output.append("    <tr>")
	for sample_subrow in header_rows[1]:
		compound, _, _, _ = sample_subrow
		output.append(f'      <td colspan="4">{compound or ""}</td>')
	output.append("    </tr>")

	# Third header row, "Experimental" and "Reference"
	output.append("    <tr>")
	output.append('      <td colspan="2">Experimental</td><td colspan="2">Reference</td>' * len(header_rows[2]))
	output.append("    </tr>")

	# Body
	for row in html_table_rows[3:]:
		output.append("    <tr>")
		for cell in itertools.chain.from_iterable(row):
			output.append(f"      <td>{cell}</td>")
		output.append("    </tr>")

	output.append("</tbody>\n  </table>")

	return '\n'.join(output)


def _fig_to_html(fig: Figure, d3_url: str, mpld3_url: str) -> str:
	# 3rd party
	from mpld3._display import GENERAL_HTML, NumpyEncoder  # type: ignore[import-untyped]

	renderer = MPLD3Renderer()
	Exporter(renderer, close_mpl=False).run(fig)

	fig, figure_json, extra_css, extra_js = renderer.finished_figures[0]

	figure_json = json.dumps(figure_json, cls=NumpyEncoder, indent=2)
	current_id = 0
	id_map = {}
	for m in re.finditer('"id": "(el.*)"', figure_json):
		element_id = m[1]
		if element_id not in id_map:
			id_map[element_id] = current_id
			current_id += 1

	for old_id, new_id in id_map.items():
		figure_json = figure_json.replace(old_id, f"el{new_id}")

	return GENERAL_HTML.render(
			figid='"mass_spectra_d3"',
			d3_url=d3_url,
			mpld3_url=mpld3_url,
			figure_json=figure_json,
			extra_css=extra_css,
			extra_js=extra_js,
			include_libraries=True
			)


def draw_mass_spectrum(
		json_data: JSONData,
		peak_number: int,
		previous_file: str,
		next_file: str,
		*,
		d3_url: str = mpld3.urls.D3_URL,
		mpld3_url: str = mpld3.urls.MPLD3MIN_URL,
		) -> str:
	"""
	Generate HTML page with mass spectra, top peaks table and similarity matrices.

	:param json_data:
	:param peak_number:
	:param previous_file:
	:param next_file:
	:param d3_url: The URL of the d3 library.  If not specified, a standard web path will be used.
	:param mpld3_url: The URL of the mpld3 library.  If not specified, a standard web path will be used.
	"""

	fig = plot_spectra(json_data)

	# fig.savefig(spectrum_basename.with_suffix(".png"))
	project_names = list(json_data["peak"])
	page_title = f"Row {json_data['row']} – " + " / ".join(project_names)

	html_output = StringList()
	html_output.set_indent_type("  ")
	html_output.extend([
			"<!DOCTYPE html>",
			"<html lang='en'>",
			"<head>",
			"  <meta charset='utf-8'>",
			f"  <title>{page_title}</title>",
			"  <style>",
			"    body {",
			'      font-family: "Segoe UI", Roboto, "Helvetica Neue", "Noto Sans", "Liberation Sans", Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";',
			"    }",
			"    .column {",
			"      float: left;",
			"      width: 50%;",
			"    }",
			"    .row:after {",
			'      content: "";',
			"      display: table;",
			"      clear: both;",
			"    }",
			"  </style>",
			"</head>",
			"<body>",
			])

	with html_output.with_indent_size(1):
		html_output.extend([
				'  <a href="index.html">Home</a>',
				f'  <a href="{previous_file}">Previous Peak</a>',
				f'  <a href="{next_file}">Next Peak</a>',
				])

		html_output.extend(_fig_to_html(fig, d3_url=d3_url, mpld3_url=mpld3_url).splitlines())

	html_output.extend([
			"  <p></p>",
			'  <div class="row">',
			'    <div class="column">',
			f'      <img src="{peak_number}.matrix.png" width="800px" height="880px">',
			"    </div>",
			'    <div class="column">',
			])

	with html_output.with_indent_size(2):
		html_output.extend(format_top_masses_html_table(get_top_masses_data(json_data)).splitlines())

	html_output.extend([
			"    </div>",
			"  </div>",
			"  <script>",
			"    document.addEventListener('keydown', (event) => {",
			"      switch (event.key) {",
			'        case "ArrowLeft":',
			f"          window.location.href = '{previous_file}';",
			"          break;",
			'        case "ArrowRight":',
			f"          window.location.href = '{next_file}';",
			"          break;",
			"      }",
			"    });",
			"  </script>",
			"</body>",
			"</html>",
			])
	plt.close(fig)
	return str(html_output)
