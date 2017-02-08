
.. image:: https://github.com/3DGenomes/tadbit/raw/master/doc/source/pictures/TADbit_logo.png
   :height: 50
   :width: 240

+-------------------------------------+---------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|                                     | .. image:: https://travis-ci.org/3DGenomes/TADbit.png?branch=master       | .. image:: https://coveralls.io/repos/github/3DGenomes/tadbit/badge.svg?branch=master :target: https://coveralls.io/github/3DGenomes/tadbit?branch=master |
| Current version: v0.2.0.4           |   :target: https://travis-ci.org/3DGenomes/TADbit                         |   :target: https://coveralls.io/github/3DGenomes/tadbit?branch=master                                                                                     |
|                                     |                                                                           |                                                                                                                                                           |
+-------------------------------------+---------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+


TADbit is a complete Python library to deal with all steps to analyze,
model and explore 3C-based data. With TADbit the user can map FASTQ
files to obtain raw interaction binned matrices (Hi-C like matrices),
normalize and correct interaction matrices, identify adn compare the
so-called Topologically Associating Domains (TADs), build 3D models
from the interaction matrices, and finally, extract structural
properties from the models. TADbit is complemented by TADkit for
visualizing 3D models.

Hi-C experiments generate genomic interaction between loci located in
the same or in different chromosomes. TADbit is built around the
concept of a chromosome, and uses it as a central item to store and
compare different Hi-C experiments. The library has been designed to
be used by researchers with no expertise in computer
science. All-in-one scripts provided in TADbit allow to run the full
analysis using one single command line; advanced users may produce
their own programs using TADbit as a complementary library.


Documentation
*************

* `TADbit documentation <http://3dgenomes.github.io/TADbit/>`_
* `Summary of TADbit classes and functions <https://github.com/3DGenomes/tadbit/blob/master/doc/summary.rst>`_

Citation
--------
Serra, F., Baù, D., Filion, G., & Marti-Renom, M. A. (2016).
**Structural features of the fly chromatin colors revealed by automatic three-dimensional modeling.**
*bioRxiv*. `doi:10.1101/036764 <http://biorxiv.org/cgi/content/short/036764>`_

Feedback
--------
If you have any question remaining, we would be happy to answer informally:

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/3DGenomes/tadbit
   :target: https://gitter.im/3DGenomes/tadbit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

TADbit training
***************

Next editions
-------------

* April 3rd to April 7th 2016: `Chromosomal Conformation course
  <http://www.crg.eu/en/event/coursescrg-chromosomal-conformation-0>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Bracelona (Spain)


Past editions
-------------

* October 10th to October 14th 2016: `3DAROC16 3C-based data analysis and 3D reconstruction of chromatin folding
  <http://gtpb.igc.gulbenkian.pt/bicourses/3DAROC16/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)
* September 28th to October 2nd 2015: `Chromosomal Conformation course
  <http://gtpb.igc.gulbenkian.pt/bicourses/2014/CSDM14/>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Bracelona (Spain)
* November 25th to November 28th 2014: `CSDM 2014
  <http://gtpb.igc.gulbenkian.pt/bicourses/2014/CSDM14/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)
* September 6th 2014: `TADbit: Automated Analysis and
  Three-Dimensional Modeling of Genomic Domains
  <http://www.eccb14.org/program/tutorials/tadbit>`_ at `ECCB14
  <http://www.eccb14.org/>`_ Strasbourg (France)
* November 27th to November 29th 2013: `CSDM 2013
  <http://gtpb.igc.gulbenkian.pt/bicourses/2013/CSDM13/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)