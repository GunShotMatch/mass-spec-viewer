#!/usr/bin/env python3
#
#  index.py
"""
Overview index.html page.
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
from typing import List, Tuple

# 3rd party
from domdf_python_tools.stringlist import StringList

__all__ = ["create_index_page"]


def create_index_page(project_names: List[str], table_data: List[Tuple[int, List[Tuple[str, str]]]]):
	# Table data: first row element is row number, 2nd is list of tuples of retention time and compound name

	page_title = " / ".join(project_names)
	index_html_output = StringList()
	index_html_output.set_indent_type("  ")
	index_html_output.extend([
			"<!DOCTYPE html>",
			"<html lang='en'>",
			"<head>",
			"  <meta charset='utf-8'>",
			f"  <title>{page_title}</title>",
			'  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-QWTKZyjpPEjISv5WaRU9OFeRpok6YctnYmDr5pNlyT2bRjXh0JMhjY6hW+ALEwIH" crossorigin="anonymous">',
			"  <style>",
			"    table.table tbody tr td.row-idx, table.table tbody tr td.compound {",
			"      border-right: 1px solid rgb(222, 226, 230);",
			"    }",
			"    tr:hover {",
			"      background: #f2f3ff;",
			"      outline: none;",
			"      cursor: pointer;",
			"    }",
			"    td {",
			"      background-color:inherit !important;",
			"    }",
			"  </style>",
			"</head>",
			"<body>",
			'  <table class="table">',
			"    <tbody>",
			"      <tr>",
			"        <th></th>",
			])

	index_html_output.indent_size = 4
	index_html_output.extend(f'<th colspan="2">{project_name}</th>' for project_name in project_names)
	index_html_output.indent_size = 3
	index_html_output.append("</tr>")
	index_html_output.indent_size = 4
	index_html_output.append("<th></th>")
	index_html_output.extend("<th>Rt</th><th>Compound</th>" for _ in project_names)
	index_html_output.indent_size = 3
	index_html_output.append("</tr>")

	for row in table_data:
		index_html_output.append(f'<tr data-href="{row[0]}.html">')
		with index_html_output.with_indent_size(4):
			index_html_output.append(f'<td class="row-idx">{row[0]}</td>')
			for (rt, compound) in row[1]:
				index_html_output.append(f'<td class="rt">{rt}</td>')
				index_html_output.append(f'<td class="compound">{compound}</td>')

		index_html_output.append("</tr>")

	index_html_output.indent_size = 0
	index_html_output.extend([
			"    </tbody>",
			"  </table>",
			"  <script>",
			"    var rows = document.getElementsByTagName('table')[0].rows;",
			"    Array.from(rows).forEach(row => {",
			"      row.addEventListener('click', function() {",
			"        window.location.href = this.getAttribute('data-href'); ",
			"      });",
			"    });",
			"  </script>",
			"</body>",
			"</html>",
			])

	return str(index_html_output)
