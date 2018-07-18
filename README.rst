
.. image:: https://github.com/3DGenomes/tadbit/raw/master/doc/source/pictures/TADbit_logo.png
   :height: 50
   :width: 240

+-------------------------------------+---------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+---------------------------------------------------------------+
|                                     | .. image:: https://travis-ci.org/3DGenomes/TADbit.png?branch=master       | .. image:: https://coveralls.io/repos/github/3DGenomes/TADbit/badge.svg?branch=master       | .. image:: https://img.shields.io/badge/license-GPL-green.svg |
| Current version: v0.2.0.472         |   :target: https://travis-ci.org/3DGenomes/TADbit                         |   :target: https://coveralls.io/github/3DGenomes/TADbit?branch=master                       |                                                               |
|                                     |                                                                           |                                                                                             |                                                               |
+-------------------------------------+---------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+---------------------------------------------------------------+


TADbit is a complete Python library to deal with all steps to analyze,
model and explore 3C-based data. With TADbit the user can map FASTsQ
files to obtain raw interaction binned matrices (Hi-C like matrices),
normalize and correct interaction matrices, identify and compare the
Topologically Associating Domains (TADs), build 3D models
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

Feedback
--------
If you have any question remaining, we would be happy to answer informally:

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/3DGenomes/tadbit
   :target: https://gitter.im/3DGenomes/tadbit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge


Docker
------

This minimal `Dockerfile <https://docs.docker.com/engine/reference/builder/>`_ (resulting in a 3 Gb docker image) can be used to run TADbit in any computer::

    FROM debian:8

    ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

    RUN apt-get update --fix-missing && \
        apt-get install -y wget bzip2 --no-install-recommends && \
        rm -rf /var/lib/apt/lists/*

    RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
        wget --quiet --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh \
            -O ~/miniconda.sh && \
        /bin/bash ~/miniconda.sh -b -p /opt/conda && \
        rm ~/miniconda.sh

    ENV PATH /opt/conda/bin:$PATH

    RUN conda config --add channels salilab && conda config --add channels bioconda && \
        conda install -y -q imp scipy matplotlib jupyter mcl samtools sra-tools pysam && \
        conda clean -y --all  && rm -rf /opt/conda/pkgs/*

    RUN wget --quiet --no-check-certificate https://newcontinuum.dl.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%202/GEM-binaries-Linux-x86_64-core_i3-20121106-022124.tbz2 \
            -O GEM.tbz2 && \
        tar xvf GEM.tbz2 && cd GEM-*/ && \
        mv * /usr/local/bin/ && cd .. && rm -rf GEM*

    RUN apt-get update --fix-missing && \
        apt-get install -y unzip build-essential --no-install-recommends && \
        wget --quiet --no-check-certificate https://github.com/fransua/TADbit/archive/dev.zip && unzip dev.zip && \
        cd TADbit-dev && python setup.py install && cd .. && rm -rf TADbit-dev dev.zip && \
        apt-get remove -y --purge unzip build-essential && \
        apt-get autoremove -y && \
        apt-get autoclean -y && \
        rm -rf /var/lib/apt/lists/*

    CMD [ "/bin/bash" ]

Build the image by saving this file as :code:`Dockerfile` into an empty folder
and build the image from inside this empty folder with :code:`docker build -t tadbit .` (~20 minutes)

Once built, run it as :code:`docker run tadbit tadbit map -h`

This image contains all dependencies for TADbit and also `jupyter <http://jupyter.org/>`_.

To run a notebook from inside the docker container run :code:`tadbit` docker image as::

    docker run -it -p 8888:8888 -v /LOCAL_PATH:/mnt tadbit

:code:`LOCAL_PATH` *would be for example a local folder with data*
*(e.g. FASTQs or reference genomes). And* :code:`/mnt` *a directory*
*inside the Docker container where the* :code:`LOCAL_PATH` *would be mounted.*

From inside docker run::

  jupyter notebook --ip 0.0.0.0 --allow-root --NotebookApp.token=''

And finally write the url :code:`http://localhost:8888` in your browser.

*Note: this can also be done in a single line and running in the background:*
::

  docker run -d -p 8888:8888 -v /LOCAL_PATH:/mnt tadbit jupyter notebook --ip 0.0.0.0 --allow-root --NotebookApp.token='' > /dev/null &


Citation
********
Please, cite this article if you use TADbit.

Serra, F., Baù, D., Goodstadt, M., Castillo, D. Filion, G., & Marti-Renom, M.A. (2017).
**Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors.**
*PLOS Comp Bio* 13(7) e1005665. `doi:10.1371/journal.pcbi.1005665 <https://doi.org/10.1371/journal.pcbi.1005665>`_

Methods implemented in TADbit
-----------------------------
In addition to the general citation for the TADbit library, please cite these articles if you used TADbit for:

- Mapping and read filtering [Imakaev2012]_ [Ay2015]_
- Hi-C normalization [Imakaev2012]_ [Rao2014]_
- A/B compartment calling [Lieberman-Aiden2009]_
- Model assessement [Trussart2015]_
- Chromatin 3D Model Building [BaùMarti-Renom2012]_

Applications
------------
TADbit has been previously used for modeling genomes and genomic domains. Here is the list of published articles:

- Alpha-globin domain [Baù2011]_
- *Caulobacter crescentus* genome [Umbarger2011]_
- TADs as regulons [Le_Dily2014]_
- Yeast chromosome III [Belton2015]_
- *Mycoplasma pneumoniae* genome [Trussart2017]_


TADbit training
***************

Next editions
-------------

* To be announced.

Past editions
-------------

* April 10th to April 11th 2017: `MuG
  <http://www.multiscalegenomics.eu/MuGVRE/>`_ workshop: `Multi-scale study of 3D Chromatin structure
  <http://www.multiscalegenomics.eu/MuGVRE/multi-scale-study-of-3d-chromatin-structure/>`_ at the
  `European Bioinformatics Institute (EMBL-EBI) <https://www.embl.de/training/cco/>`_,
  Hinxton, Cambridge, (United Kingdom)
* April 3rd to April 7th 2017: `Chromosomal Conformation course
  <http://www.crg.eu/en/event/coursescrg-chromosomal-conformation-0>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Barcelona (Spain)
* October 10th to October 14th 2016: `3DAROC16 3C-based data analysis and 3D reconstruction of chromatin folding
  <http://gtpb.igc.gulbenkian.pt/bicourses/3DAROC16/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)
* September 28th to October 2nd 2015: `Chromosomal Conformation course
  <http://gtpb.igc.gulbenkian.pt/bicourses/2014/CSDM14/>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Barcelona (Spain)
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


Bibliography
************

.. [Ay2015] Ay, F., Vu, T.H., Zeitz, M.J., Varoquaux, N., Carette, J.E., Vert, J.-P., Hoffman, A.R. and Noble, W.S. 2015. Identifying multi-locus chromatin contacts in human cells using tethered multiple 3C. BMC Genomics 16, p. 121.

.. [BaùMarti-Renom2012] Baù, D. and Marti-Renom, M.A. 2012. Genome structure determination via 3C-based data integration by the Integrative Modeling Platform. Methods 58(3), pp. 300–306.

.. [Baù2011] Baù, D., Sanyal, A., Lajoie, B.R., Capriotti, E., Byron, M., Lawrence, J.B., Dekker, J. and Marti-Renom, M.A. 2011. The three-dimensional folding of the α-globin gene domain reveals formation of chromatin globules. Nature Structural & Molecular Biology 18(1), pp. 107–114.

.. [Belton2015] Belton, J.-M., Lajoie, B.R., Audibert, S., Cantaloube, S., Lassadi, I., Goiffon, I., Baù, D., Marti-Renom, M.A., Bystricky, K. and Dekker, J. 2015. The conformation of yeast chromosome III is mating type dependent and controlled by the recombination enhancer. Cell reports 13(9), pp. 1855–1867.

.. [Imakaev2012] Imakaev, M., Fudenberg, G., McCord, R.P., Naumova, N., Goloborodko, A., Lajoie, B.R., Dekker, J. and Mirny, L.A. 2012. Iterative correction of Hi-C data reveals hallmarks of chromosome organization. Nature Methods 9(10), pp. 999–1003.

.. [Le_Dily2014] Le Dily, F., Baù, D., Pohl, A., Vicent, G.P., Serra, F., Soronellas, D., Castellano, G., Wright, R.H.G., Ballare, C., Filion, G., Marti-Renom, M.A. and Beato, M. 2014. Distinct structural transitions of chromatin topological domains correlate with coordinated hormone-induced gene regulation. Genes & Development 28(19), pp. 2151–2162.

.. [Lieberman-Aiden2009] Lieberman-Aiden, E., van Berkum, N.L., Williams, L., Imakaev, M., Ragoczy, T., Telling, A., Amit, I., Lajoie, B.R., Sabo, P.J., Dorschner, M.O., Sandstrom, R., Bernstein, B., Bender, M.A., Groudine, M., Gnirke, A., Stamatoyannopoulos, J., Mirny, L.A., Lander, E.S. and Dekker, J. 2009. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 326(5950), pp. 289–293.

.. [Rao2014] Rao, S.S.P., Huntley, M.H., Durand, N.C., Stamenova, E.K., Bochkov, I.D., Robinson, J.T., Sanborn, A.L., Machol, I., Omer, A.D., Lander, E.S. and Aiden, E.L. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159(7), pp. 1665–1680.

.. [Trussart2015] Trussart, M., Serra, F., Baù, D., Junier, I., Serrano, L. and Marti-Renom, M.A. 2015. Assessing the limits of restraint-based 3D modeling of genomes and genomic domains. Nucleic Acids Research 43(7), pp. 3465–3477.

.. [Trussart2017] Trussart, M., Yus, E., Martinez, S., Baù, D., Tahara, Y.O., Pengo, T., Widjaja, M., Kretschmer, S., Swoger, J., Djordjevic, S., Turnbull, L., Whitchurch, C., Miyata, M., Marti-Renom, M.A., Lluch-Senar, M. and Serrano, L. 2017. Defined chromosome structure in the genome-reduced bacterium Mycoplasma pneumoniae. Nature Communications 8, p. 14665.

.. [Umbarger2011] Umbarger, M.A., Toro, E., Wright, M.A., Porreca, G.J., Baù, D., Hong, S.-H., Fero, M.J., Zhu, L.J., Marti-Renom, M.A., McAdams, H.H., Shapiro, L., Dekker, J. and Church, G.M. 2011. The three-dimensional architecture of a bacterial genome and its alteration by genetic perturbation. Molecular Cell 44(2), pp. 252–264.