---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Conformance Status

Which requirements are met, which are not, and how each claim is checked.

```{code-cell} python
import itacart
report = itacart.conformance()
report["summary"]
```

Twenty-four requirements. Nineteen met in full, four partially, one not at all.
**Core conformant; EAERS not.**

```{figure} ../../../_static/f9/f9_17_conformance_board.png
:alt: The conformance board
:width: 100%

`docs/_static/f9/f9_17_conformance_board.png`
```

## The five divergences

```{code-cell} python
[(r["requirement"], r["coverage"]) for r in report["requirements"]
 if r["coverage"] != "full"]
```

Each carries its own justification, and the justifications are the substance:

```{code-cell} python
import textwrap
for r in report["requirements"]:
    if r["coverage"] != "full":
        print(f"--- requirement {r['requirement']} ({r['coverage']}) ---")
        print(textwrap.fill(r["justification"], 76))
        print()
```

## Formally non-conformant, substantively close

Requirement 28 is the clearest case of the distinction
[the overview](ogc_iso_19170_overview.md) draws. The budget is one per cent or
less against the theoretical average cell area — at resolution 1 that average
is 100 008 533 square metres over 5 100 221 cells.

Formal conformance is **no**, because a budget covering the whole domain cannot
be stated. Substantive coverage is **partial**, and the exception is small and
named: 3 941 antemeridian trapezoids and the 2 polar caps, which is 3 943 cells
or 0.077 per cent of the grid.

Every other cell — parallelogram and meridian triangle alike — sits at 0.999915
of the theoretical average, which is 0.0085 per cent off and well inside the
budget. The trapezoid family runs from 0.032 to 2.067 and the caps sit at
0.121; those figures describe **that family**, not the grid.

```{figure} ../../../_static/f9/f9_20_equal_area_budget.png
:alt: The equal-area budget against the theoretical average
:width: 100%

`docs/_static/f9/f9_20_equal_area_budget.png`
```

Requirement 27 is the same shape of argument. `cell_to_centroid` lies on the
surface and within the cell everywhere. Over the ordinary families it agrees
with the geodesic centre of surface area to within ten metres — on a cell ten
kilometres across, a thousandth of the cell — and over the parallelogram
families to 0.7 metres, stable under refinement of the search. It is
nonetheless the area centroid of the equal-area projection plane inverted to
geodetic, and the requirement names the geodesic centre. The two part company
at the polar row and the two caps, where the cap's geodesic centre is the pole
itself and the plane answer sits about 1.3 kilometres short: three orders of
magnitude, in three families out of ten.

**A *shall* is not met by seven.** That sentence is the whole reason these are
reported as non-conformant rather than as rounding.

```{figure} ../../../_static/f9/f9_18_centroid_against_geodesic_centre.png
:alt: The centroid against the geodesic centre of surface area
:width: 100%

`docs/_static/f9/f9_18_centroid_against_geodesic_centre.png`
```

## How the claims are checked

The table above is generated, not written. `itacart.conformance()` reads the
package's own record, and a conformance test suite exercises the requirements
in continuous integration — so a claim that stops being true fails a build
rather than surviving in a document.

That is also why this page shows the report instead of restating it. There is
one copy of the conformance position, and you are reading it.

## Reading the numbering

The requirement numbers here are the standard's, as the implementation records
them, and they do not line up item for item with the table printed in the
paper. Where the two differ, this report is the one that is checked.
