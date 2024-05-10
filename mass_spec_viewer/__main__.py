#!/usr/bin/env python3
#
#  __main__.py
"""
Combined spectrum viewer and comparison tool.
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

# 3rd party
import click
import consolekit
from consolekit.options import flag_option

__all__ = ["main"]


@click.argument("comparison")
@click.option("-c", "--config", help="Configuration TOML file", default="viewer.toml")
@click.option("-j", "--jobs", help="Number of parallel workers.", default=1)
# @flag_option("-D", "--download-assets", "should_download_assets", help="Download CSS and JS assets for future offline use.", default=False)
@flag_option(
		"-C",
		"--copy-assets",
		"should_copy_assets",
		help="Copy CSS and JS assets for future offline use.",
		default=False
		)
@consolekit.click_command()
def main(
		comparison: str,
		config: str = "viewer.toml",
		jobs: int = 1,
		should_copy_assets: bool = False,
		) -> None:
	"""
	Generate HTML combined mass spectrum comparison.
	"""

	# 3rd party
	import matplotlib  # type: ignore[import-untyped]
	matplotlib.use("Agg")

	# 3rd party
	from domdf_python_tools.paths import PathPlus
	from gunshotmatch_pipeline.utils import tomllib

	# this package
	from mass_spec_viewer._batch_process import process_comparison

	config_file = PathPlus(config)
	viewer_config = tomllib.loads(config_file.read_text())

	comparison = comparison.strip()
	if comparison == ":all:":
		tables = viewer_config.keys()
	else:
		tables = [t.strip() for t in comparison.split(',')]

	for table_name in tables:

		print(table_name)
		process_comparison(
				viewer_config[table_name],
				jobs=jobs,
				should_copy_assets=should_copy_assets,
				)


if __name__ == "__main__":
	main()
