ITACaRT
=======

**ITA Cadastral Ellipsoidal Reference Tessellation** — an equal-area
parallelogram Discrete Global Grid System for terrestrial cadastral
mapping, tessellated directly on the WGS84 ellipsoid.

.. code-block:: python

   import itacart

   cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
   lon, lat = itacart.cell_to_centroid(cell)

.. Every document below appears in exactly one toctree. A document listed
   twice is reported at "checking consistency" and produces a duplicated
   sidebar entry, so a new page goes in one place and is cross-referenced
   from the others.

.. toctree::
   :maxdepth: 3
   :caption: Getting Started

   getting_started/index

.. toctree::
   :maxdepth: 3
   :caption: Concepts & Guides

   concepts/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

.. The phase galleries are organised by the development phase that produced
   each figure, which is how the work was done rather than how the grid is
   learned. They stay built and linkable — the concept pages reference the
   figures they need — but they are not a section a reader navigates by, so
   they are not in the sidebar.

.. toctree::
   :hidden:

   _generated/figures/index

Citing
------

Silva, I. N., Dietzsch, G., & Shiguemori, E. H. (2025). ITACaRT: An
Equal-Area Parallelogram Discrete Global Grid System for Terrestrial
Cadastral Mapping — Designed for Usability and Blockchain Integration.
*Revista Brasileira de Cartografia*, 77.
https://doi.org/10.14393/rbcv77n0a-79281

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
